# Makefile to compile OGSPlugins
#
# This makefile helps compile the OGS ParaView suite and 
# deploy bit.sea into ParaView.
#
# To add a new plugin just copy paste the compiling rule 
# and rename it according to your plugin folder. Don't 
# forget to add the plugin into the plugins rule for the
# general compilation.
#
# Arnau Miro - OGS (2020)

# TODO: Tests

# Optimization, host and CPU type
#
OPTL = 2
HOST = Host
TUNE = skylake

# Options
#
VECTORIZATION  = ON
OPENMP_PARALL  = ON
FORCE_GCC      = OFF
DEBUGGING      = OFF

# Paths to the plugins
#
BIN_PATH  = bin
BPS_PATH  = bit.sea
PVPL_PATH = OGSPlugins
SPB_PATH  = superbuild

# Versions of the libraries
#
PARAVIEW_VERS = 5.6.3
LAPACK_VERS   = 3.9.0
PROJ_VERS     = 5.2.0
PROJ_DATV     = 1.8
GEOS_VERS     = 3.7.0
QT5_VERS      = 5.10
OSX_SDK       = macosx10.14

# Path to ParaView
#
PV_BIN_PATH = $$(dirname $$(which paraview-config))
PV_LIB_PATH = $(PV_BIN_PATH)/../lib
PV_PYT_PATH = $$(echo $(PV_LIB_PATH)/python*/site-packages)
PV_PLG_PATH = $(PV_LIB_PATH)/plugins

# Compilers
#
ifeq ($(FORCE_GCC),ON) 
	# Forcing the use of GCC
	# C Compiler
	CC = gcc
	# C++ Compiler
	CXX = g++
	# Fortran Compiler
	FC = gfortran
else
	ifeq (,$(shell which icc))
		# C Compiler
		CC = gcc
		# C++ Compiler
		CXX = g++
		# Fortran Compiler
		FC = gfortran
	else
		# C Compiler
		CC = icc
		# C++ Compiler
		CXX = icpc
		# Fortran Compiler
		FC = ifort
	endif
endif

# Compiler flags
#
ifeq ($(CC),gcc)
	# Using GCC as a compiler
	ifeq ($(DEBUGGING),ON)
		# Debugging flags
		CFLAGS   += -O0 -g -rdynamic -fPIC
		CXXFLAGS += -O0 -g -rdynamic -fPIC
		FFLAGS   += -O0 -g -rdynamic -fPIC
	else
		CFLAGS   += -O$(OPTL) -fPIC
		CXXFLAGS += -O$(OPTL) -fPIC
		FFLAGS   += -O$(OPTL) -fPIC
	endif
	# Vectorization flags
	ifeq ($(VECTORIZATION),ON)
		CFLAGS   += -ffast-math -march=native -ftree-vectorize
		CXXFLAGS += -ffast-math -march=native -ftree-vectorize
		FFLAGS   += -ffast-math -march=native -ftree-vectorize
	endif
	# OpenMP flag
	ifeq ($(OPENMP_PARALL),ON)
		CFLAGS   += -fopenmp
		CXXFLAGS += -fopenmp
	endif
else
	# Using INTEL as a compiler
	ifeq ($(DEBUGGING),ON)
		# Debugging flags
		CFLAGS   += -O0 -g -traceback -fPIC
		CXXFLAGS += -O0 -g -traceback -fPIC
		FFLAGS   += -O0 -g -traceback -fPIC
	else
		CFLAGS   += -O$(OPTL) -fPIC
		CXXFLAGS += -O$(OPTL) -fPIC
		FFLAGS   += -O$(OPTL) -fPIC
	endif
	# Vectorization flags
	ifeq ($(VECTORIZATION),ON)
		CFLAGS   += -x$(HOST) -mtune=$(TUNE)
		CXXFLAGS += -x$(HOST) -mtune=$(TUNE)
		FFLAGS   += -x$(HOST) -mtune=$(TUNE)
	endif
	# OpenMP flag
	ifeq ($(OPENMP_PARALL),ON)
		CFLAGS   += -qopenmp
		CXXFLAGS += -qopenmp
	endif
endif
# C++ standard
CXXFLAGS += -std=c++11

# Generic launchers
#
default: tools bit.sea plugins

firstime: launchers tools bit.sea plugins

# OGS Launchers
#
launchers: img2video paraview-launch pvpython-egl pvserver-egl
	@echo ""
	@echo "   Thanks for waiting! OGS ParaView launchers have been "
	@echo "   successfully installed."
	@echo ""

img2video: $(BIN_PATH)/img2video
	@echo "$< -> $(PV_BIN_PATH)/$@"
	@rm -f $(PV_BIN_PATH)/$@
	@cp $(PWD)/$< $(PV_BIN_PATH)/$@
	@chmod +x $(PV_BIN_PATH)/$@

paraview-launch: $(BIN_PATH)/paraview-launch
	@echo "$< -> $(PV_BIN_PATH)/$@"
	@rm -f $(PV_BIN_PATH)/$@
	@cp $(PWD)/$< $(PV_BIN_PATH)/$@
	@chmod +x $(PV_BIN_PATH)/$@

pvpython-egl: $(BIN_PATH)/pvpython-egl
	@echo "$< -> $(PV_BIN_PATH)/$@"
	@rm -f $(PV_BIN_PATH)/$@
	@cp $(PWD)/$< $(PV_BIN_PATH)/$@
	@chmod +x $(PV_BIN_PATH)/$@

pvserver-egl: $(BIN_PATH)/pvserver-egl
	@echo "$< -> $(PV_BIN_PATH)/$@"
	@rm -f $(PV_BIN_PATH)/$@
	@cp $(PWD)/$< $(PV_BIN_PATH)/$@
	@chmod +x $(PV_BIN_PATH)/$@

# OGS Tools
#
tools: OGS2Paraview OGSlonlat2m OGSmesh ParaViewBlender
	@echo ""
	@echo "   Thanks for waiting! OGS ParaView tools have been "
	@echo "   successfully installed."
	@echo ""	

libOGS.so: $(PVPL_PATH)/_utils/OGS.cpp
	$(CXX) -shared -Wl,-soname,$@ $(CXXFLAGS) $< -o $@ -DOGS_NO_NETCDF

libOGS.dylib: $(PVPL_PATH)/_utils/OGS.cpp
	$(CXX) -shared $(CXXFLAGS) $< -o $@ -DOGS_NO_NETCDF

libOGS_install: libOGS.so
	@mv $< $(PV_LIB_PATH)

OGS2Paraview: $(PVPL_PATH)/_utils/python/OGS2Paraview.py
	@echo "$< -> $(PV_PYT_PATH)/$@.py"
	@rm -f $(PV_PYT_PATH)/$@.py
	@rm -f $(PV_BIN_PATH)/$@.py
	@cp $(PWD)/$< $(PV_PYT_PATH)/$@.py
	@cp $(PWD)/$< $(PV_PYT_PATH)/default.ini
	@ln -s $(PV_PYT_PATH)/$@.py $(PV_BIN_PATH)/$@.py
	@chmod +x $(PV_BIN_PATH)/$@.py

OGSlonlat2m: $(PVPL_PATH)/_utils/python/OGSlonlat2m.py
	@echo "$< -> $(PV_PYT_PATH)/$@.py"
	@rm -f $(PV_PYT_PATH)/$@.py
	@rm -f $(PV_BIN_PATH)/$@.py
	@cp $(PWD)/$< $(PV_PYT_PATH)/$@.py
	@ln -s $(PV_PYT_PATH)/$@.py $(PV_BIN_PATH)/$@.py
	@chmod +x $(PV_BIN_PATH)/$@.py

OGSmesh: $(PVPL_PATH)/_utils/python/OGSmesh.py libOGS_install
	@echo "$< -> $(PV_PYT_PATH)/$@.py"
	@rm -f $(PV_PYT_PATH)/$@.py
	@cp $(PWD)/$< $(PV_PYT_PATH)/$@.py

ParaViewBlender: $(PVPL_PATH)/_utils/python/ParaViewBlender.py
	@echo "$< -> $(PV_PYT_PATH)/$@.py"
	@rm -f $(PV_PYT_PATH)/$@.py
	@cp $(PWD)/$< $(PV_PYT_PATH)/$@.py

# OGS bit.sea
#
bit.sea: commons basins MapPlotter
	@echo ""
	@echo "   Thanks for waiting! bit.sea has been "
	@echo "   successfully installed."
	@echo ""

commons: $(BPS_PATH)/commons
	@echo "$< -> $(PV_PYT_PATH)/$@"
	@rm -rf $(PV_PYT_PATH)/$@
	@cp -r $(PWD)/$< $(PV_PYT_PATH)/$@

basins: $(BPS_PATH)/basins
	@echo "$< -> $(PV_PYT_PATH)/$@"
	@rm -rf $(PV_PYT_PATH)/$@
	@cp -r $(PWD)/$< $(PV_PYT_PATH)/$@

MapPlotter: $(BPS_PATH)/MapPlotter
	@echo "$< -> $(PV_PYT_PATH)/$@"
	@rm -rf $(PV_PYT_PATH)/$@
	@cp -r $(PWD)/$< $(PV_PYT_PATH)/$@

# Prerequisites
#
prereq: lapack proj
	@echo ""
	@echo "   Thanks for waiting! Prerequisites have been "
	@echo "   successfully installed."
	@echo ""
lapack: $(PVPL_PATH)/_utils/lapack/
	@bash $</install_lapack.sh "${LAPACK_VERS}" "${PWD}/$<" "${CC}" "${CFLAGS}" "${FC}" "${FFLAGS}"
proj: $(PVPL_PATH)/_utils/proj/
	@bash $</install_proj.sh "${PROJ_VERS}" "${PROJ_DATV}" "STATIC" "${PWD}/$<" "${CC}" "${CFLAGS}" "${CXX}" "${CXXFLAGS}"
geos: $(PVPL_PATH)/_utils/geos/
	@bash $</install_geos.sh "${GEOS_VERS}" "${PWD}/$<" "${CC}" "${CFLAGS}" "${CXX}" "${CXXFLAGS}"
prereq-osx: /usr/local/bin/pkg-config /usr/local/bin/openssl
	@echo ""
	@echo "   Thanks for waiting! Prerequisites for MacOS X"
	@echo "   have been successfully installed."
	@echo ""
/usr/local/bin/pkg-config: $(SPB_PATH)/fix_pkgconfig_osx.sh
	@bash $<
/usr/local/bin/openssl: $(SPB_PATH)/fix_openssl_osx.sh
	@bash $<

# ParaView Plugins
#
ALL_PLUGINS = $(notdir $(shell find $(PVPL_PATH)/* -maxdepth 0 -type d))
DEV_PLUGINS = $(notdir $(shell find $(PVPL_PATH)/_* -maxdepth 0 -type d))
PLUGINS     = $(filter-out $(DEV_PLUGINS),$(ALL_PLUGINS))

plugins: plugindir prereq $(PLUGINS)
	@echo ""
	@echo "   Thanks for waiting! OGS ParaView plugins have been "
	@echo "   successfully installed."
	@echo ""
	@echo "   Have a nice day!"
plugindir:
	@echo ""
	@echo "   Installation of OGS ParaView Plugins"
	@echo "   ____________________________________"
	@echo ""
	@echo "   The process will take several minutes."
	@echo "   You can go grab a coffee meanwhile."
	@echo ""
	@mkdir -p $(PV_PLG_PATH)

# The following rule finds all plugins in the OGSPlugins folder and 
# automatically links them to the Makefile.
#
# Plugins starting with "_" are considered under development and thus
# are not considered for compilation in this rule.

$(PLUGINS): $(addprefix $(PVPL_PATH)/,$(PLUGINS))
	@echo "$(PVPL_PATH)/$@ -> $(PV_PLG_PATH)"
	-@(cd $(PVPL_PATH)/$@ && mkdir -p build)
	@(cd $(PVPL_PATH)/$@/build && rm -rf * && cmake .. )
	@(cd $(PVPL_PATH)/$@/build && make )
	-@(cd $(PVPL_PATH)/$@/build && cp *.so $(PV_PLG_PATH) )
	-@(cd $(PVPL_PATH)/$@ && rm -rf build)

# Superbuild compilations
#
superbuild-linux: $(SPB_PATH)/paraview_superbuild_linux_gen.sh prereq libOGS.so 
	@echo ""
	@echo "   OGS ParaView v${PARAVIEW_VERS} Superbuild"
	@echo "   _________________________________________"
	@echo ""
	@echo "   The process will take some time."
	@echo "   You can go grab a coffee or two meanwhile..."
	@echo ""

	@bash $< "${PARAVIEW_VERS}" "${QT5_VERS}" "${PROJ_VERS}" "${PROJ_DATV}" "${GEOS_VERS}" "${CC}" "${CFLAGS}" "${CXX}" "${CXXFLAGS}"
	
	@echo ""
	@echo "   Thanks for waiting! OGS ParaView superbuild has been "
	@echo "   successfully compiled with the OGS plugins."
	@echo ""

superbuild-osx: $(SPB_PATH)/paraview_superbuild_osx.sh prereq-osx prereq libOGS.dylib
	@echo ""
	@echo "   OGS ParaView v${PARAVIEW_VERS} Superbuild"
	@echo "   _________________________________________"
	@echo ""
	@echo "   The process will take some time."
	@echo "   You can go grab a coffee or two meanwhile..."
	@echo ""

	@bash $< "${PARAVIEW_VERS}" "${OSX_SDK}" "${QT5_VERS}" "${PROJ_VERS}" "${PROJ_DATV}" "${GEOS_VERS}" "${CC}" "${CFLAGS}" "${CXX}" "${CXXFLAGS}" "${FC}" "${FFLAGS}"
	
	@echo ""
	@echo "   Thanks for waiting! OGS ParaView superbuild has been "
	@echo "   successfully compiled with the OGS plugins."
	@echo ""

superbuild-galileo: $(SPB_PATH)/paraview_superbuild_hpc.sh prereq libOGS.so 
	@echo ""
	@echo "   OGS ParaView v${PARAVIEW_VERS} Superbuild"
	@echo "   _________________________________________"
	@echo ""
	@echo "   The process will take some time."
	@echo "   You can go grab a coffee or two meanwhile..."
	@echo ""

	@bash $< "${PARAVIEW_VERS}" "GALILEO" "${QT5_VERS}" "${PROJ_VERS}" "${PROJ_DATV}" "${GEOS_VERS}" "${CC}" "${CFLAGS}" "${CXX}" "${CXXFLAGS}" "${FC}" "${FFLAGS}"
	
	@echo ""
	@echo "   Thanks for waiting! OGS ParaView superbuild has been "
	@echo "   successfully compiled with the OGS plugins."
	@echo ""

superbuild-marconi: $(SPB_PATH)/paraview_superbuild_hpc.sh prereq libOGS.so 
	@echo ""
	@echo "   OGS ParaView v${PARAVIEW_VERS} Superbuild"
	@echo "   _________________________________________"
	@echo ""
	@echo "   The process will take some time."
	@echo "   You can go grab a coffee or two meanwhile..."
	@echo ""

	@bash $< "${PARAVIEW_VERS}" "MARCONI" "${QT5_VERS}" "${PROJ_VERS}" "${PROJ_DATV}" "${GEOS_VERS}" "${CC}" "${CFLAGS}" "${CXX}" "${CXXFLAGS}" "${FC}" "${FFLAGS}"
	
	@echo ""
	@echo "   Thanks for waiting! OGS ParaView superbuild has been "
	@echo "   successfully compiled with the OGS plugins."
	@echo ""


# Generic object makers
#
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $< $(DFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(DFLAGS)

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $< $(DFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(DFLAGS)

# Cleaning rules 
#
clean_launchers:
	-@rm -rf $(PV_BIN_PATH)/img2video
	-@rm -rf $(PV_BIN_PATH)/*-launch
	-@rm -rf $(PV_BIN_PATH)/*-egl
clean_libogs.so:
	-@rm -rf $(PV_LIB_PATH)/libOGS.so
clean_tools: clean_libogs.so
	-@rm -f $(PV_PYT_PATH)/OGS2Paraview.py
	-@rm -f $(PV_BIN_PATH)/OGS2Paraview.py
	-@rm -f $(PV_PYT_PATH)/OGSlonlat2m.py
	-@rm -f $(PV_BIN_PATH)/OGSlonlat2m.py
	-@rm -f $(PV_PYT_PATH)/OGSmesh.py
	-@rm -f $(PV_PYT_PATH)/ParaViewBlender.py
clean_bit.sea:
	-@rm -rf $(PV_PYT_PATH)/commons
	-@rm -rf $(PV_PYT_PATH)/basins
	-@rm -rf $(PV_PYT_PATH)/MapPlotter
clean_plugins:
	-@rm -rf $(PV_PLG_PATH)/libMapPlotterView.so
	-@rm -rf $(PV_PLG_PATH)/libOGS*

.PHONY: clean		
clean: clean_tools clean_bit.sea clean_plugins
	@echo "OGS ParaView Suite successfully cleaned!"

.PHONY: uninstall		
uninstall: clean_launchers clean_tools clean_bit.sea clean_plugins
	@echo "OGS ParaView Suite successfully uninstalled!"
