#!/bin/bash

module load daint-gpu
module load cdt/19.08
module switch PrgEnv-cray PrgEnv-pgi
#module switch PrgEnv-cray PrgEnv-gnu
module load craype-accel-nvidia60
#module load PyExtensions/3.6.5.1-CrayGNU-18.08
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
