#!/bin/bash
#
# If cartopy doesn't work, this fix might be able to
# make it work again by recompiling from source
#
# Arnau Miro, OGS 2020

VERS=v0.18.0

source env.sh

# Clone git repository
git clone -b $VERS https://github.com/SciTools/cartopy.git cartopy
cd cartopy
# Compile
LDFLAGS=-L$PWD/../lib pvpython setup.py build
# Copy it to the destination
rm -rf ../lib/python2.7/site-packages/cartopy
cp -r build/lib.*/cartopy ../lib/python2.7/site-packages/
# Cleanup and exit
cd ..
rm -rf cartopy