# This maintains the links for all sources used by this superbuild.
# Simply update this file to change the revision.
# One can use different revision on different platforms.
# e.g.
# if (UNIX)
#   ..
# else (APPLE)
#   ..
# endif()

include(CMakeDependentOption)

# NOTE: if updating bzip2 version, fix patch in bzip2.cmake
superbuild_set_revision(bzip2
  URL     "https://www.paraview.org/files/dependencies/bzip2-1.0.6.tar.gz"
  URL_MD5 00b516f4704d4a7cb50a1d97e6e8e15b)

superbuild_set_revision(zlib
#  URL     "https://www.paraview.org/files/dependencies/zlib-1.2.11.tar.xz"
#  URL_MD5 85adef240c5f370b308da8c938951a68)
  URL     "https://www.paraview.org/files/dependencies/zlib-1.2.7.tar.gz"
  URL_MD5 60df6a37c56e7c1366cca812414f7b85)

superbuild_set_revision(ffmpeg
  URL     "https://www.paraview.org/files/dependencies/ffmpeg-2.3.3.tar.bz2"
  URL_MD5 72361d3b8717b6db3ad2b9da8df7af5e)

superbuild_set_revision(szip
  URL     "https://www.paraview.org/files/dependencies/szip-2.1.1.tar.gz"
  URL_MD5 dd579cf0f26d44afd10a0ad7291fc282)

superbuild_set_revision(hdf5
  URL     "https://www.paraview.org/files/dependencies/hdf5-1.10.5.tar.bz2"
  URL_MD5 7c19d6b81ee2a3ba7d36f6922b2f90d3)

superbuild_set_revision(boost
  URL     "https://www.paraview.org/files/dependencies/boost_1_59_0.tar.bz2"
  URL_MD5 6aa9a5c6a4ca1016edd0ed1178e3cb87)

superbuild_set_revision(png
  URL     "https://www.paraview.org/files/dependencies/libpng-1.6.37.tar.xz"
  URL_MD5 015e8e15db1eecde5f2eb9eb5b6e59e9)

if (WIN32 AND (NOT superbuild_building_prebuilt_python OR superbuild_use_prebuilt_python))
  if (superbuild_is_64bit)
    superbuild_set_revision(python
      URL     "https://www.paraview.org/files/dependencies/python-2.7.15-win64-20180905.tar.gz"
      URL_MD5 6cfab07945bf75474d4ed2d2ea799c57)
  else ()
    message(FATAL_ERROR
      "Prebuilt Python binaries for Windows 32 bit are not provided.")
  endif ()
else()
  superbuild_set_revision(python
    URL     "https://www.paraview.org/files/dependencies/Python-2.7.15.tar.xz"
    URL_MD5 a80ae3cc478460b922242f43a1b4094d)
endif()

if (WIN32)
  superbuild_set_revision(python3
    URL "https://www.paraview.org/files/dependencies/Python-3.7.4-win64-no-ssl.tar.xz"
    URL_MD5 c73455c755f2c79031fc92036040a90b)
else()
  superbuild_set_revision(python3
    URL     "https://www.paraview.org/files/dependencies/Python-3.7.4.tar.xz"
    URL_MD5 d33e4aae66097051c2eca45ee3604803)
endif()

superbuild_set_revision(ftjam
  URL     "https://www.paraview.org/files/dependencies/ftjam-2.5.2-win32.tar.bz2"
  URL_MD5 ee52f3faff6d31ffb89a2fedb3b0caf6)

superbuild_set_revision(freetype
#  URL     "https://www.paraview.org/files/dependencies/freetype-2.10.0.tar.bz2"
#  URL_MD5 a717e6925b61b9dda946322ecd278a42)
  URL     "https://www.paraview.org/files/dependencies/freetype-2.8.1.tar.bz2"
  URL_MD5 bf0a210b6fe781228fa0e4a80691a521)

superbuild_set_revision(gperf
  URL     "https://www.paraview.org/files/dependencies/gperf-3.1.tar.gz"
  URL_MD5 9e251c0a618ad0824b51117d5d9db87e)

superbuild_set_revision(fontconfig
#  URL     "https://www.paraview.org/files/dependencies/fontconfig-2.13.1.tar.bz2"
#  URL_MD5 36cdea1058ef13cbbfdabe6cb019dc1c)
  URL     "https://www.paraview.org/files/dependencies/fontconfig-2.12.6.tar.bz2"
  URL_MD5 733f5e2371ca77b69707bd7b30cc2163)

superbuild_set_revision(libxml2
#  URL     "https://www.paraview.org/files/dependencies/libxml2-2.9.9.tar.gz"
#  URL_MD5 c04a5a0a042eaa157e8e8c9eabe76bd6)
  URL     "https://www.paraview.org/files/dependencies/libxml2-2.7.8.tar.gz"
  URL_MD5 8127a65e8c3b08856093099b52599c86)

superbuild_set_revision(nlohmannjson
  URL     "https://www.paraview.org/files/dependencies/nlohmannjson-3.6.1.tar.gz"
  URL_MD5 c53592d55e7fec787cf0a406d36098a3)

if (WIN32)
  set(qt5_ext "zip")
else ()
  set(qt5_ext "tar.xz")
endif ()
set(qt5_8_ver "5.8.0")
if (WIN32)
  set(qt5_8_md5 "1e372fabc9d97a32877cb4adb377e7c8")
else ()
  set(qt5_8_md5 "66660cd3d9e1a6fed36e88adcb72e9fe")
endif ()
set(qt5_9_ver "5.9.2")
if (WIN32)
  set(qt5_9_md5 "d5239e19f6b80dcf44f4dd2de04c7d3d")
else ()
  set(qt5_9_md5 "738d1b98106e1bd39f00cc228beb522a")
endif ()
set(qt5_10_ver "5.10.1")
if (WIN32)
  set(qt5_10_md5 "60c4ea41950857c65015fb6cffcb2497")
else ()
  set(qt5_10_md5 "7e167b9617e7bd64012daaacb85477af")
endif ()
set(qt5_12_ver "5.12.3")
if (WIN32)
  set(qt5_12_md5 "5da2e14a9f5db620c578fafd219592d7")
else ()
  set(qt5_12_md5 "38017e0ed88b9baba063bd723d9ca8b2")
endif ()
set(qt5_13_ver "5.13.2")
set(qt5_13_md5 "7c04c678d4ecd9e9c06747e7c17e0bb9")
set(qt5_14_ver "5.14.2")
set(qt5_14_md5 "b3d2b6d00e6ca8a8ede6d1c9bdc74daf")
set(qt5_15_ver "5.15.0")
set(qt5_15_md5 "610a228dba6ef469d14d145b71ab3b88")
superbuild_set_selectable_source(qt5
  SELECT 5.8
    URL     "https://www.paraview.org/files/dependencies/qt-everywhere-opensource-src-${qt5_8_ver}.${qt5_ext}"
    URL_MD5 "${qt5_8_md5}"
  SELECT 5.9
    URL     "https://www.paraview.org/files/dependencies/qt-everywhere-opensource-src-${qt5_9_ver}.${qt5_ext}"
    URL_MD5 "${qt5_9_md5}"
  SELECT 5.10
    URL     "https://www.paraview.org/files/dependencies/qt-everywhere-src-${qt5_10_ver}.${qt5_ext}"
    URL_MD5 "${qt5_10_md5}"
  SELECT 5.12
    URL     "https://www.paraview.org/files/dependencies/qt-everywhere-src-${qt5_12_ver}.${qt5_ext}"
    URL_MD5 "${qt5_12_md5}"
  SELECT 5.13
    URL     "http://download.qt.io/archive/qt/5.13/5.13.0/single/qt-everywhere-src-${qt5_13_ver}.${qt5_ext}" 
    URL_MD5 "${qt5_13_md5}"
  SELECT 5.14
    URL     "https://download.qt.io/archive/qt/5.14/5.14.2/single/qt-everywhere-src-${qt5_14_ver}.${qt5_ext}" 
    URL_MD5 "${qt5_14_md5}"
  SELECT 5.15 DEFAULT
    URL     "https://download.qt.io/archive/qt/5.15/5.15.0/single/qt-everywhere-src-${qt5_15_ver}.${qt5_ext}" 
    URL_MD5 "${qt5_15_md5}")

if (WIN32 AND NOT superbuild_building_prebuilt_python)
  if (superbuild_is_64bit)
    superbuild_set_revision(numpy
      URL     "https://www.paraview.org/files/dependencies/numpy-1.15.1-win64-20180906.tar.gz"
      URL_MD5 d75f1c5c111de3fed8556174fe353f0c)
  else ()
    message(FATAL_ERROR
      "Prebuilt Python binaries for Windows 32 bit are not provided.")
  endif ()
else ()
  superbuild_set_revision(numpy
    URL     "https://www.paraview.org/files/dependencies/numpy-1.16.4.tar.gz"
    URL_MD5 6edf7334d04d8e8849ad058ccd3b3803)
  superbuild_set_revision(scipy
    URL     "https://www.paraview.org/files/dependencies/scipy-1.2.2.tar.xz"
    URL_MD5 136c5ee1bc4b259a12a7efe331b15d64)
endif ()

if (WIN32 AND NOT superbuild_building_prebuilt_python AND NOT ENABLE_python3 AND NOT python3_enabled)
  if (superbuild_is_64bit)
    superbuild_set_revision(matplotlib
      URL     "https://www.paraview.org/files/dependencies/matplotlib-1.1.1-win64-20180905.tar.gz"
      URL_MD5 0c96b84e87b4db50cdc4d18869ae74ed)
  else ()
    message(FATAL_ERROR
      "Prebuilt Python binaries for Windows 32 bit are not provided.")
  endif ()
else ()
  if (ENABLE_python3 OR python3_enabled)
    superbuild_set_revision(matplotlib
      URL "https://www.paraview.org/files/dependencies/matplotlib-3.1.1-no-jquery.tar.xz"
      URL_MD5 20f4c3a9ba1df1d6f7b2a509f617ddb1)
  else()
    superbuild_set_revision(matplotlib
      URL     "https://files.pythonhosted.org/packages/eb/a0/31b6ba00bc4dcbc06f0b80d1ad6119a9cc3081ecb04a00117f6c1ca3a084/matplotlib-2.2.3.tar.gz"
      URL_MD5 403b0bddd751d71187416f20d4cff100)
  endif()
endif ()

if (WIN32 AND NOT superbuild_building_prebuilt_python AND NOT ENABLE_python3 AND NOT python3_enabled)
  if (superbuild_is_64bit)
    superbuild_set_revision(pywin32
      URL     "https://www.paraview.org/files/dependencies/pywin32-220-win64-20180905.tar.gz"
      URL_MD5 08a6ab778e459e6752d54083c29dbb13)
  else ()
    message(FATAL_ERROR
      "Prebuilt Python binaries for Windows 32 bit are not provided.")
  endif ()
elseif (ENABLE_python3 OR python3_enabled)
  superbuild_set_revision(pywin32
    URL "https://www.paraview.org/files/dependencies/pywin32-224-amd64.tar.xz"
    URL_MD5 1a79e2571d276eeb26f54d944daf39e4)
else ()
  superbuild_set_revision(pywin32
    URL     "https://www.paraview.org/files/dependencies/pywin32-220.zip"
    URL_MD5 9c386839c1485b2047c03fab66e69b9e)
endif ()

superbuild_set_revision(mpi
  URL     "https://www.paraview.org/files/dependencies/mpich-3.3.tar.gz"
  URL_MD5 574af413dc0dc7fbb929a761822beb06)

superbuild_set_revision(lapack
  URL     "https://github.com/Reference-LAPACK/lapack/archive/v3.9.0.tar.gz"
  URL_MD5 0b251e2a8d5f949f99b50dd5e2200ee2)

superbuild_set_revision(netcdf
  URL     "https://www.paraview.org/files/dependencies/netcdf-c-4.7.0.tar.gz"
  URL_MD5 37134a12a49e80c45fb58777aa3e9e3b)

# Using Intel Threading Building Blocks 2018 Update 2
set(tbb_ver "2018_20171205oss")
if (WIN32)
  set(tbb_file "tbb${tbb_ver}_win.zip")
  set(tbb_md5 "e37abf02d74a638f7a6629c992f23918")
elseif (APPLE)
  # using older version on Mac since newer version has issues
  # with changing install name.
  set(tbb_ver "44_20150728oss")
  set(tbb_file "tbb${tbb_ver}_osx.tgz")
  set(tbb_md5 "a767d7a8b375e6b054e44e2317d806b8")
else ()
  set(tbb_file "tbb${tbb_ver}_lin.tgz")
  set(tbb_md5 "d637d29f59ee31fe5830a0366e2e973a")
endif ()

superbuild_set_revision(tbb
  URL     "https://www.paraview.org/files/dependencies/${tbb_file}"
  URL_MD5 "${tbb_md5}")

superbuild_set_revision(pytz
  URL     "https://www.paraview.org/files/dependencies/pytz-2016.10.tar.bz2"
  URL_MD5 88b1d6c50c764579292edce3493c8a3a)

superbuild_set_revision(pythondateutil
  URL     "https://www.paraview.org/files/dependencies/python-dateutil-2.6.0.tar.gz"
  URL_MD5 6e38f91e8c94c15a79ce22768dfeca87)

superbuild_set_revision(pythonpyparsing
  URL     "https://www.paraview.org/files/dependencies/pyparsing-2.2.0.tar.gz"
  URL_MD5 0214e42d63af850256962b6744c948d9)

superbuild_set_revision(pythoncycler
  URL     "https://www.paraview.org/files/dependencies/cycler-0.10.0.tar.gz"
  URL_MD5 4cb42917ac5007d1cdff6cccfe2d016b)

superbuild_set_revision(pythonsetuptools
  URL     "https://www.paraview.org/files/dependencies/setuptools-23.0.0.tar.gz"
  URL_MD5 100a90664040f8ff232fbac02a4c5652)

superbuild_set_revision(pythonautobahn
  URL     "https://www.paraview.org/files/dependencies/autobahn-17.10.1.tar.gz"
  URL_MD5 f8c8d74bf73644719b751e6fb11dc4a3)

superbuild_set_revision(pythonconstantly
  URL     "https://www.paraview.org/files/dependencies/constantly-15.1.0.tar.gz"
  URL_MD5 f0762f083d83039758e53f8cf0086eef)

superbuild_set_revision(pythonhyperlink
  URL     "https://www.paraview.org/files/dependencies/hyperlink-17.3.1.tar.gz"
  URL_MD5 eaccb9845b559817e838846669cbc68a)

superbuild_set_revision(pythonincremental
  URL     "https://www.paraview.org/files/dependencies/incremental-17.5.0.tar.gz"
  URL_MD5 602746e0d438e075a5a9e0678140bba2)

superbuild_set_revision(pythontwisted
  URL     "https://www.paraview.org/files/dependencies/Twisted-19.2.1.tar.bz2"
  URL_MD5 528b7856938edc2a752c244aebd94981)

superbuild_set_revision(pythontxaio
  URL     "https://www.paraview.org/files/dependencies/txaio-2.8.2.tar.gz"
  URL_MD5 a20e3431c95896a49fa3ffa84f77cde1)

superbuild_set_revision(pythonwslink
  URL     "https://www.paraview.org/files/dependencies/wslink-0.1.11.tar.gz"
  URL_MD5 35e6285c2a74304da0557f1402c400e5)

superbuild_set_revision(pythonzope
  URL     "https://www.paraview.org/files/dependencies/Zope-4.0b3.tar.gz"
  URL_MD5 9a63e8c8b614dc6d6944fcbd9c105f45)

superbuild_set_revision(pythonzopeinterface
  URL     "https://www.paraview.org/files/dependencies/zope.interface-4.4.3.tar.gz"
  URL_MD5 8700a4f527c1203b34b10c2b4e7a6912)

superbuild_set_revision(pythonsix
  URL     "https://www.paraview.org/files/dependencies/six-1.11.0.tar.gz"
  URL_MD5 d12789f9baf7e9fb2524c0c64f1773f8)

superbuild_set_revision(pythonpygments
  URL     "https://www.paraview.org/files/dependencies/Pygments-2.2.0.tar.gz"
  URL_MD5 13037baca42f16917cbd5ad2fab50844)

superbuild_set_revision(pythonmako
  URL     "https://www.paraview.org/files/dependencies/Mako-1.0.7.tar.gz"
  URL_MD5 5836cc997b1b773ef389bf6629c30e65)

superbuild_set_revision(pythonkiwisolver
  URL     "https://www.paraview.org/files/dependencies/kiwisolver-1.1.0.tar.gz"
  URL_MD5 fc8a614367f7ba0d34a02fd08c535afc)

superbuild_set_revision(pythonattrs
  URL     "https://www.paraview.org/files/dependencies/attrs-19.1.0.tar.gz"
  URL_MD5 2be7bce157988928f5ff2bb50a0b510d)

superbuild_set_revision(ffi
  URL     "https://www.paraview.org/files/dependencies/libffi-3.3.tar.gz"
  URL_MD5 6313289e32f1d38a9df4770b014a2ca7)

superbuild_set_revision(utillinux
  URL     "https://www.paraview.org/files/dependencies/util-linux-2.34.tar.xz"
  URL_MD5 a78cbeaed9c39094b96a48ba8f891d50)