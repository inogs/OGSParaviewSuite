if (BUILD_SHARED_LIBS)
  set(python_shared_args --enable-shared)
else ()
  set(python_shared_args --disable-shared --enable-static)
endif ()

if (NOT _python_unicode_default)
  set(_python_unicode_default "UCS2")
endif ()

set(python_USE_UNICODE "${_python_unicode_default}" CACHE STRING "Enable Unicode support for python")
set_property(CACHE python_USE_UNICODE PROPERTY STRINGS "OFF;UCS2;UCS4")
mark_as_advanced(python_USE_UNICODE)

if (python_USE_UNICODE STREQUAL "UCS2")
  set(python_unicode_args "--enable-unicode=ucs2")
elseif (python_USE_UNICODE STREQUAL "UCS4")
  set(python_unicode_args "--enable-unicode=ucs4")
else ()
  set(python_unicode_args "--disable-unicode")
endif ()

#superbuild_apply_patch(python ssl
#    "SSL with macports")

superbuild_add_project(python
  CAN_USE_SYSTEM
  DEPENDS bzip2 zlib png
  CONFIGURE_COMMAND
    <SOURCE_DIR>/configure
      --prefix=<INSTALL_DIR>
      --with-ensurepip=install
      ${python_unicode_args}
      ${python_shared_args}
  BUILD_COMMAND
    $(MAKE)
  INSTALL_COMMAND
    make install)

if (NOT CMAKE_CROSSCOMPILING)
  # Pass the -rpath flag when building python to avoid issues running the
  # executable we built.
  superbuild_append_flags(
    ld_flags "-Wl,-rpath,${superbuild_install_location}/lib"
    PROJECT_ONLY)
endif ()

if (python_enabled)
  set(superbuild_python_executable "${superbuild_install_location}/bin/python"
    CACHE INTERNAL "")
  set(superbuild_python_pip "${superbuild_install_location}/bin/pip"
    CACHE INTERNAL "")
else ()
  set(superbuild_python_executable ""
    CACHE INTERNAL "")
  set(superbuild_python_pip ""
    CACHE INTERNAL "")
endif ()

superbuild_add_extra_cmake_args(
  -DPYTHON_EXECUTABLE:FILEPATH=<INSTALL_DIR>/bin/python2.7
  -DPYTHON_INCLUDE_DIR:PATH=<INSTALL_DIR>/include/python2.7
  -DPYTHON_LIBRARY:FILEPATH=<INSTALL_DIR>/lib/libpython2.7.dylib)
