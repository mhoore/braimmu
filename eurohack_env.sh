#!/bin/bash

module load daint-gpu
module switch PrgEnv-cray PrgEnv-gnu
module load craype-accel-nvidia60
module load PyExtensions/3.6.5.1-CrayGNU-18.08

