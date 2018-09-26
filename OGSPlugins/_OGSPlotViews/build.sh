#!/bin/bash
#

rm -rf build
mkdir build && cd build
cmake .. -DPARAVIEW_QT_VERSION=5
make
cd ..
