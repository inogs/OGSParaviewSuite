#!/bin/bash
#
# Load OGS ParaView Suite into your environment
#
DIR=$PWD

export PATH=$DIR/bin:$PATH
export DYLD_FALLBACK_LIBRARY_PATH=$DIR/Libraries:$DIR/Libraries/paraview-5.6:$DYLD_FALLBACK_LIBRARY_PATH
#export C_INCLUDE_PATH=$DIR/include:$C_INCLUDE_PATH
export PYTHONPATH=$DIR/Python:$PYTHONPATH