#!/bin/bash

time=1
freq=100000
dt=0.001

sec0=(90)
sec1=(90)
sec2=(80)

for ((i=0; i<1; i+=1))
do

  tdum=10000
  for ((t=0; t<=${time}; t+=${freq}))
  do
    echo "processing time ${t} sec ${sec0[$i]}_${sec1[$i]}_${sec2[$i]}"
    realtime=$(echo ${t} ${dt} | awk '{printf "%4.3f\n",$1*$2}')
    input="../../example/dump.${t}.nii"
    input0="../../example/dump.0.nii"

    python3 multiplaner.py ${input} ${sec0} ${sec1} ${sec2} ${tdum} ${realtime} ${input0}
    tdum=$((${tdum}+1))
done
  
  # movies
  mkdir gifs
  
  #ffmpeg -framerate 2 -start_number 10000 -i multi/m_Z_neu_s${sec0[$i]}_${sec1[$i]}_${sec2[$i]}_t%d.png -qscale 24 -r 2 -pix_fmt yuv420p gifs/m_Z_neu_s${sec0[$i]}_${sec1[$i]}_${sec2[$i]}.gif
  #ffmpeg -framerate 2 -start_number 10000 -i multi/m_Z_mic_s${sec0[$i]}_${sec1[$i]}_${sec2[$i]}_t%d.png -qscale 24 -r 2 -pix_fmt yuv420p gifs/m_Z_mic_s${sec0[$i]}_${sec1[$i]}_${sec2[$i]}.gif
  #ffmpeg -framerate 2 -start_number 10000 -i multi/m_Z_ast_s${sec0[$i]}_${sec1[$i]}_${sec2[$i]}_t%d.png -qscale 24 -r 2 -pix_fmt yuv420p gifs/m_Z_ast_s${sec0[$i]}_${sec1[$i]}_${sec2[$i]}.gif
  #ffmpeg -framerate 2 -start_number 10000 -i multi/m_Z_fAb_s${sec0[$i]}_${sec1[$i]}_${sec2[$i]}_t%d.png -qscale 24 -r 2 -pix_fmt yuv420p gifs/m_Z_fAb_s${sec0[$i]}_${sec1[$i]}_${sec2[$i]}.gif
  #ffmpeg -framerate 2 -start_number 10000 -i multi/m_Z_sAb_s${sec0[$i]}_${sec1[$i]}_${sec2[$i]}_t%d.png -qscale 24 -r 2 -pix_fmt yuv420p gifs/m_Z_sAb_s${sec0[$i]}_${sec1[$i]}_${sec2[$i]}.gif
  #ffmpeg -framerate 2 -start_number 10000 -i multi/m_atrophy_s${sec0[$i]}_${sec1[$i]}_${sec2[$i]}_t%d.png -qscale 24 -r 2 -pix_fmt yuv420p gifs/m_atrophy_s${sec0[$i]}_${sec1[$i]}_${sec2[$i]}.gif

done
