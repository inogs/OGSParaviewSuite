#!/bin/bash
#
# Set environment for ParaView

OGSSUITEDIR=$PWD

module purge
module load profile/advanced
module load gnu/6.1.0 
module load openmpi/3.0.0--gnu--6.1.0

export PATH=$PATH:$OGSSUITEDIR/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OGSSUITEDIR/lib:$OGSSUITEDIR/lib/paraview-5.6
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$DIR/include:$DIR/include/python2.7:$DIR/include/paraview-5.6
export PYTHONPATH=$PYTHONPATH:$OGSSUITEDIR/lib/python2.7/site-packages