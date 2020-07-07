if (BUILD_SHARED_LIBS)
  set(freetype_shared_args --enable-shared --disable-static)
else ()
  set(freetype_shared_args --disable-shared --enable-static)
endif ()

superbuild_add_project(freetype
  CAN_USE_SYSTEM
  DEPENDS zlib
  CONFIGURE_COMMAND
    <SOURCE_DIR>/configure
      --prefix=<INSTALL_DIR>
      --with-harfbuzz=no
      ${shared_args}
      --with-sysroot=<INSTALL_DIR>
  BUILD_COMMAND
    $(MAKE)
  INSTALL_COMMAND
    $(MAKE) install)

if (APPLE AND __BUILDBOT_INSTALL_LOCATION)
  superbuild_project_add_step(clean-build
    COMMAND   make
              clean
    DEPENDEES configure
    DEPENDERS build
    COMMENT   "Cleaning the build tree for install name fixes"
    WORKING_DIRECTORY <BINARY_DIR>)
endif ()
