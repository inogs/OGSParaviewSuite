# OGS ParaView Suite

The **OGS ParaView Suite** contains ParaView filters to work with data from [OGSTM-BFM](http://marine.copernicus.eu/services-portfolio/access-to-products/). The suite also links with *bit.sea* and other python libraries for pre- and post-processing of the data. Each filter has its own readme, containing a brief explanation of itself. The code can be obtained by either compiling the filters in a working development build of ParaView (as instructed later) or downloading the latest [Release](https://github.com/inogs/OGSParaviewSuite/releases).

The project can be cloned from github by doing
```bash
git clone --recurse-submodules https://github.com/inogs/OGSParaviewSuite.git
```
or 
```bash
git clone https://github.com/inogs/OGSParaviewSuite.git
cd OGSParaviewSuite
git submodule update --init 
```

## Plugin List

A list of the available plugins is as follows:
1. [OGS Reader](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSReader)
2. [OGS Annotate Date Time](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSAnnotateDateTime)
3. [OGS Select Tools](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSSelectTools)
4. [OGS Variable Aggregator](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSVariableAggregator)
5. [OGS Spatial Statistics](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSSpatialStatistics)
6. [OGS Time Statistics](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSTimeStatistics)
7. [OGS Derivatives](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSDerivatives)
8. [OGS Vortex Identification](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSVortexIdentification)
9. [OGS Vortex Detection](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSVortexDetection)
10. [OGS Depth Profile](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSDepthProfile)
11. [OGS Hovmoeller](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSHovmoeller)
12. [OGS Spaghetti](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSSpaghetti)
13. [OGS Utilities](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSUtils)
14. [OGS Plot Views](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSPlotViews)
15. [OGS Writer](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSWriter)
16. [OGS Compare Variables](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSCompareVariables)
17. [OGS Water Density](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSDensity)
18. [OGS Mixing Layer Depth](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSMixingLayerDepth)
19. [OGS Rossby Radius](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/OGSRossbyRadius)
20. [Map Plotter View](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/MapPlotterView)

Details on the Field and VTK interface and the conversion from UVW grid to T grid can be found [here](https://github.com/inogs/OGSParaviewSuite/tree/master/OGSPlugins/_utils).

## Compilation

The provided **Makefile** enables easy compilation in UNIX environments. Inside the **Makefile** there are a number of options for the user to configure. These refer to the compilation of the external libraries needed for the **OGS ParaView Suite**. The compilation options for the plugins are inside *OGSPlugins/OGSPlugins.cmake*.

_NOTE_: Activating the **VECTORIZATION** option can generate code that is machine dependent. Use this option with caution.

_NOTE_: MacOS X Clang does not support **OPENMP** parallelization. 

### Compiling the OGS ParaView Suite

In a proper environment with [ParaView](https://www.paraview.org/download/) and [CMake](https://cmake.org/download/), the user only needs to use the command
```bash
make
```
which triggers the compilation of the whole suite. For a first time installation, the command 
```bash
make firstime
```
copies also the launchers into the ParaView directories. Finally, the command 
```bash
make plugins
```
compiles only the OGS ParaView plugins. 

### Superbuild compilation

The [superbuild](https://gitlab.kitware.com/paraview/paraview-superbuild/) installation provides with an autonomous, standalone build of ParaView and its dependencies along with the OGS ParaView plugins. Superbuild can be built in Linux or in MacOS X. In both cases, it is necessary to install the dependencies of gfortran, CMake and OpenSSL.

Moreover, before launching any superbuild compilation, make sure that **VECTORIZATION** and **OPENMP** are deactivated inside the Makefile and the OGSPlugins.cmake file.

#### Linux

The dependencies can be easily installed as
```bash
sudo apt install libssl-dev gfortran
```
Compilation of the superbuild can then be achieved through the make command
```bash
make superbuild-linux
```
_NOTE_: to compile in **GALILEO** or **MARCONI** use
```bash
make superbuild-galileo
```
or
```bash
make superbuild-marconi
```

#### MacOS X

Compilation in MacOS X requires of a clean environment with [gfortran](https://github.com/fxcoudert/gfortran-for-macOS/releases) and [CMake](https://cmake.org/download/) installed. They can be obtained from their respective websites. Regarding the CMake app, it is necessary to export the _PATH_ environmental variable as
```bash
export PATH=$PATH:/Applications/CMake.app/Contents/bin
```
Then, the prerequisites for MacOS X (OpenSSL and pkg-config) can be installed using the make tool as
```bash
make prereq-osx
```
Finally, the compilation of the ParaView superbuild can be started as
```bash
make superbuild-osx
```
