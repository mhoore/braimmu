#!/bin/bash

module load daint-gpu
module switch PrgEnv-cray PrgEnv-gnu
module load craype-accel-nvidia60
module switch cudatoolkit cudatoolkit/9.2.148_3.19-6.0.7.1_2.1__g3d9acc8
module load PyExtensions/3.6.5.1-CrayGNU-18.08

export CPATH=$CRAY_MPICH_DIR/include
