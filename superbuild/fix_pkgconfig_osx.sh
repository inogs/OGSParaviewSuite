#!/bin/bash
#
# Prerequisites to compile superbuild in a 
# MAC OS X environment.
#
# Arnau Miro, OGS 2020

PKG_VERS=0.29.2

PKG_DIR="pkg-config-${PKG_VERS}"
PKG_TAR="${PKG_DIR}.tar.gz"
PKG_URL="https://pkg-config.freedesktop.org/releases/${PKG_TAR}"

curl -o $PKG_TAR $PKG_URL
tar xvzf $PKG_TAR
cd $PKG_DIR

LDFLAGS="-framework CoreFoundation -framework Carbon" ./configure --with-internal-glib
make
sudo make install

cd ..
rm -rf $PKG_DIR $PKG_TAR