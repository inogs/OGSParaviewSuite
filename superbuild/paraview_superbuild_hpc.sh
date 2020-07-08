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
# Arnau Miro, OGS 2020

PV_VERS=${1}
MACHINE=${2}
QT5_VERS=${3}
PROJ_VERS=${4}
PROJ_DATV=${5}
GEOS_VERS=${6}
CCOMPILER=${7}
CFLAGS=${8}
CXXCOMPILER=${9}
CXXFLAGS=${10}

NPROCS=8
INSTALL_PREFIX="${PWD}/../paraview-${PV_VERS}"
SUITEDIR="${PWD}"
SUPERBUILD_DIR="${PWD}/../paraview-superbuild"
BUILD_DIR="${PWD}/../paraview-build"

VERSIONSDIR="${SUITEDIR}/superbuild/versions"
MAPPLOTLIBDIR="${SUITEDIR}/superbuild/matplotlib_py27"
PLUGINDIR="${SUITEDIR}/OGSPlugins/"
BINDIR="${SUITEDIR}/bin"
BITSEADIR="${SUITEDIR}/bit.sea"

# Load modules - provide a basic building environment
module purge
if [ "$MACHINE" = "GALILEO" ]; then
   module load profile/advanced
   module load gnu/7.3.0
   module load openmpi/3.1.4--gnu--7.3.0
   module load blas/3.8.0--gnu--7.3.0
   module load cmake git
fi
if [ "$MACHINE" = "MARCONI" ]; then
   module load profile/advanced
   module load gnu/6.1.0 
   module load openmpi/3.0.0--gnu--6.1.0
   module load blas/3.6.0--gnu--6.1.0
   module load cmake git
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
# FIX: matplotlib
cp $MAPPLOTLIBDIR/matplotlib.cmake superbuild/projects
cp $MAPPLOTLIBDIR/matplotlib-kiwisolver.patch superbuild/projects/patches
# FIX: libpng
cp $SUITEDIR/superbuild/projects/png.cmake superbuild/projects
# FIX: freetype marconi
if [ "$MACHINE" = "MARCONI" ]; then
   cp $SUITEDIR/superbuild/projects/freetype.cmake superbuild/projects/apple-unix
fi
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
   -DENABLE_mpi=ON \
   -DUSE_SYSTEM_mpi=OFF \
   -DENABLE_cuda=OFF \
   -DENABLE_boost=OFF \
   -DENABLE_lapack=OFF \
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
   -DENABLE_vrpn=OFF \
   -DENABLE_pvNVIDIAIndeX=ON \
   -DENABLE_vortexfinder2=OFF \
   -DENABLE_paraview=ON \
   -DENABLE_paraviewsdk=OFF \
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
printf "Deploying ParaView $PV_VERS in $INSTALL_PREFIX.\n"
#read -s -p "Please check install configuration and press [enter]..."

# Make
make -j $NPROCS
# Do some fixing...
ln -s $BUILD_DIR/install/lib/python2.7/site-packages/matplotlib-*/matplotlib/ $BUILD_DIR/install/lib/python2.7/site-packages
# Install
make -j $NPROCS install
# FIX: matplotlib depends on backports and kiwisolver
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/backports.*/backports $INSTALL_PREFIX/lib/python2.7/site-packages/
cp $BUILD_DIR/install/lib/python2.7/site-packages/kiwisolver-*/kiwisolver.so $INSTALL_PREFIX/lib/python2.7/site-packages/

# Load environment
export PATH=$PATH:$INSTALL_PREFIX/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_PREFIX/lib:$INSTALL_PREFIX/lib/paraview-5.6
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$INSTALL_PREFIX/include
export PYTHONPATH=$PYTHONPATH:$INSTALL_PREFIX/lib/python2.7/site-packages

# Deploy GEOS library
bash $PLUGINDIR/_utils/geos/install_geos.sh "${GEOS_VERS}" "${INSTALL_PREFIX}" "${CCOMPILER}" "${CFLAGS}" "${CXXCOMPILER}" "${CXXFLAGS}"
bash $PLUGINDIR/_utils/geos/install_geos.sh "${GEOS_VERS}" "${INSTALL_PREFIX}" "${BUILD_DIR}/install" "${CFLAGS}" "${CXXCOMPILER}" "${CXXFLAGS}"

# Deploy PROJ library
bash $PLUGINDIR/_utils/proj/install_proj.sh "${PROJ_VERS}" "${PROJ_DATV}" "SHARED" "${INSTALL_PREFIX}" "${CCOMPILER}" "${CFLAGS}" "${CXXCOMPILER}" "${CXXFLAGS}"
bash $PLUGINDIR/_utils/proj/install_proj.sh "${PROJ_VERS}" "${PROJ_DATV}" "SHARED" "${BUILD_DIR}/install" "${CCOMPILER}" "${CFLAGS}" "${CXXCOMPILER}" "${CXXFLAGS}"

# Install netCDF4, configparser, cython
$BUILD_DIR/install/bin/pip install --upgrade pip
$BUILD_DIR/install/bin/pip install netcdf4 configparser cython cartopy pyepsg
cp -r $BUILD_DIR/install/include/python* $INSTALL_PREFIX/include
cp -r $BUILD_DIR/install/include/paraview-* $INSTALL_PREFIX/include
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/cftime $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/netCDF4 $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/configparser.* $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/backports $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/backports.functools_lru_cache-1.6.1-py2.7.egg/backports/* $INSTALL_PREFIX/lib/python2.7/site-packages/backports/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/cython.* $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/Cython $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/lib-dynload $INSTALL_PREFIX/lib/python2.7/
cp -r $BUILD_DIR/install/lib/python2.7/lib-tk $INSTALL_PREFIX/lib/python2.7/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/cartopy $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/shape* $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/pyepsg.* $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/requests $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/urllib3 $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/chardet $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/certifi $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/idna $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/numpy-*/numpy/core/include $INSTALL_PREFIX/lib/python2.7/site-packages/numpy/core

# Copy ffmpeg
cp $BUILD_DIR/install/bin/ffmpeg $INSTALL_PREFIX/bin
cd ..
printf "Install done successfully!\n"

# Deploy img2video
printf "Deploying img2video... "
cp $BINDIR/img2video $INSTALL_PREFIX/bin
printf "OK\n"

# Deploy bit.sea inside the ParaView installation
printf "Deploying bit.sea... "
cp -r $BITSEADIR/commons    $INSTALL_PREFIX/lib/python*/site-packages
cp -r $BITSEADIR/basins     $INSTALL_PREFIX/lib/python*/site-packages
cp -r $BITSEADIR/MapPlotter $INSTALL_PREFIX/lib/python*/site-packages
printf "OK\n"

# Deploy OGSMesh and OGS2Paraview inside the installation
printf "Deploying OGSMesh and OGS2Paraview... "
mv $SUITEDIR/libOGS.so $INSTALL_PREFIX/lib
cp $PLUGINDIR/_utils/python/OGSmesh.py $INSTALL_PREFIX/lib/python*/site-packages
cp $PLUGINDIR/_utils/python/OGSlonlat2m.py $INSTALL_PREFIX/bin
cp $PLUGINDIR/_utils/python/OGS2Paraview.py $INSTALL_PREFIX/bin
cp $PLUGINDIR/_utils/python/default.ini $INSTALL_PREFIX/bin
printf "OK\n"

if [ "$MACHINE" = "GALILEO" ]; then
   cp $SUITEDIR/superbuild/env-galileo.sh $INSTALL_PREFIX/env.sh
fi
if [ "$MACHINE" = "MARCONI" ]; then
   cp $SUITEDIR/superbuild/env-marconi.sh $INSTALL_PREFIX/env.sh
fi

# Clean-up
cd $SUITEDIR
rm -rf $BUILD_DIR $SUPERBUILD_DIR
tar cvzf "${INSTALL_PREFIX}.tar.gz" $INSTALL_PREFIX
