if ((WIN32) OR (APPLE))
  list(APPEND qt5_options
    -qt-zlib)
else ()
  list(APPEND qt5_depends
    zlib)
  list(APPEND qt5_options
    -system-zlib)
endif ()

set(qt5_EXTRA_CONFIGURATION_OPTIONS ""
    CACHE STRING "Extra arguments to be passed to Qt when configuring.")
mark_as_advanced(qt5_EXTRA_CONFIGURATION_OPTIONS)

set(qt5_configure_ext)
if (WIN32)
  set(qt5_configure_ext ".bat")
endif ()

set(qt5_build_commands
  BUILD_COMMAND   $(MAKE)
  INSTALL_COMMAND make install)
if (WIN32)
  if ((NOT CMAKE_GENERATOR MATCHES "^NMake.*$") OR
      (NOT CMAKE_GENERATOR MATCHES "^Visual Studio.*$"))
    find_program(NMAKE_PATH nmake)
  endif ()

  set(qt5_build_commands
    BUILD_COMMAND   ${NMAKE_PATH}
    INSTALL_COMMAND ${NMAKE_PATH} install)
endif ()

# If not using system qt5, add qt5_ENABLE_OPENSSL option
option(qt5_ENABLE_OPENSSL
  "Build with OpenSSL support. Requires system-installed OpenSSL at runtime." OFF)
mark_as_advanced(qt5_ENABLE_OPENSSL)
if (qt5_ENABLE_OPENSSL)
  # Require build machines to have OpenSSL
  find_package(OpenSSL)
  if (NOT OpenSSL_FOUND)
    message(FATAL_ERROR "Cannot build with qt5_ENABLE_OPENSSL option because OpenSSL not found")
  endif ()
  list(APPEND qt5_options "-openssl-linked")
else ()
  list(APPEND qt5_options "-no-openssl")
endif ()

# Add option to build qtsvg
option(qt5_ENABLE_SVG "Build Qt5 SVG library." OFF)
mark_as_advanced(qt5_ENABLE_SVG)
if (NOT qt5_ENABLE_SVG)
  list(APPEND qt5_options
    -skip qtsvg)
endif()

foreach(module IN LISTS qt5_skip_modules)
  list(APPEND qt5_skip_args -skip ${module})
endforeach()

superbuild_add_project(qt5
  CAN_USE_SYSTEM
  DEPENDS ${qt5_depends} ${qt5_extra_depends} cxx11
  CONFIGURE_COMMAND
    <SOURCE_DIR>/configure${qt5_configure_ext}
      -opensource
      -confirm-license

      -release

      -prefix <INSTALL_DIR>
      -I <INSTALL_DIR>/include
      -L <INSTALL_DIR>/lib

      ${qt5_skip_args}

      -nomake examples
      -nomake tests

      -no-dbus

      -qt-libjpeg
      -qt-pcre

      ${qt5_options}
      ${qt5_extra_options}
      ${qt5_EXTRA_CONFIGURATION_OPTIONS}
  ${qt5_build_commands}
  ${qt5_process_environment})

superbuild_add_extra_cmake_args(
  -DPARAVIEW_QT_VERSION:STRING=5
  -DQt5_DIR:PATH=<INSTALL_DIR>/lib/cmake/Qt5)
