#!/bin/bash
#
# Load OGS ParaView Suite into your environment
#
DIR=$PWD

export PATH=$PATH:$DIR/bin:$DIR/MacOS
export DYLD_FALLBACK_LIBRARY_PATH=$DYLD_FALLBACK_LIBRARY_PATH:$DIR/Libraries
#export C_INCLUDE_PATH=$C_INCLUDE_PATH:$DIR/include
export PYTHONHOME=$DIR
export PYTHONPATH=$PYTHONPATH:/usr/lib/python2.7:$DIR/Python