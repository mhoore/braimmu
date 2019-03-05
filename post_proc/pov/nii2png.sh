#!/bin/bash
make clean
make -j 4

mkdir images

minval=0
maxval=0
section=54

input="../../dumps/dump.1000000.nii"

./n2p.exe ${input} images/image.mic.pov 0 ${minval} ${maxval} ${section} &
./n2p.exe ${input} images/image.neu.pov 1 1.85e7 3.7e7 ${section} &
./n2p.exe ${input} images/image.sAb.pov 2 ${minval} ${maxval} ${section} &
./n2p.exe ${input} images/image.fAb.pov 3 ${minval} ${maxval} ${section} &
./n2p.exe ${input} images/image.act.pov 4 ${minval} ${maxval} ${section}

povray -H2000 -W2000 +I"images/image.mic.pov" +A Display=0 &
povray -H2000 -W2000 +I"images/image.neu.pov" +A Display=0 &
povray -H2000 -W2000 +I"images/image.sAb.pov" +A Display=0 &
povray -H2000 -W2000 +I"images/image.fAb.pov" +A Display=0 &
povray -H2000 -W2000 +I"images/image.act.pov" +A Display=0

rm images/image.*.pov

rm n2p.exe
make clean
