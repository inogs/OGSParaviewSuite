#!/bin/bash
#
# Set environment for ParaView

OGSSUITEDIR=$PWD

module purge
module load profile/advanced
module load gnu/6.1.0
module load openmpi/3.1.1--gnu--6.1.0

export PATH=$OGSSUITEDIR/bin:$PATH
export LD_LIBRARY_PATH=$OGSSUITEDIR/lib:$OGSSUITEDIR/lib/paraview-5.6:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$DIR/include:$DIR/include/python2.7:$DIR/include/paraview-5.6:$C_INCLUDE_PATH
export PYTHONPATH=$OGSSUITEDIR/lib/python2.7/site-packages:$PYTHONPATH