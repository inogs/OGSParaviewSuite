#!/bin/bash
#
# Load OGS ParaView Suite into your environment
#
DIR=$PWD

export PATH=$PATH:$DIR/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DIR/lib:$DIR/lib/paraview-5.6
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$DIR/include:$DIR/include/python2.7:$DIR/include/paraview-5.6
export PYTHONPATH=$PYTHONPATH:$DIR/lib/python2.7/site-packages