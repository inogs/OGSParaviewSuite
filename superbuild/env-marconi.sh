#!/bin/bash
#
# Set environment for ParaView

OGSSUITEDIR=$PWD

module purge
module load gnu/7.3.0 
module load openmpi/3.0.0--gnu--7.3.0

export PATH=$PATH:$OGSSUITEDIR/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OGSSUITEDIR/lib:$OGSSUITEDIR/lib/paraview-5.6
export PYTHONPATH=$PYTHONPATH:$OGSSUITEDIR/lib/python2.7/site-packages