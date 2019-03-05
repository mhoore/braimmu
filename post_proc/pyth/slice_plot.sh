#!/bin/bash

time=6000000
freq=100000
dt=0.001

dim=1
sec=(70 90 100)

for ((i=0; i<4; i+=1))
do

  tdum=10000
  for ((t=0; t<=${time}; t+=${freq}))
  do
    echo "processing time ${t} sec ${sec[$i]} dim ${dim}"
    realtime=$(echo ${t} ${dt} | awk '{printf "%4.3f\n",$1*$2}')
    input="../../dumps/dump.${t}.nii"
    
    python3 slice_plot.py ${input} ${dim} ${sec[$i]} ${tdum} ${realtime}
    tdum=$((${tdum}+1))
  done
  
  # movies
  mkdir gifs
  
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_mic_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_mic_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_neu_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_neu_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_sAb_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_sAb_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_fAb_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_fAb_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_ast_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_ast_d${dim}_s${sec[$i]}.gif
  
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_Z_mic_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_Z_mic_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_Z_neu_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_Z_neu_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_Z_sAb_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_Z_sAb_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_Z_fAb_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_Z_fAb_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_Z_ast_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_Z_ast_d${dim}_s${sec[$i]}.gif

  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_corr_mic_neu_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_corr_mic_neu_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_corr_neu_sAb_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_corr_neu_sAb_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_corr_neu_fAb_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_corr_neu_fAb_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_corr_neu_ast_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_corr_neu_ast_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_contours/s_cont_corr_mic_fAb_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_cont_corr_mic_fAb_d${dim}_s${sec[$i]}.gif

  ffmpeg -framerate 5 -start_number 10000 -i slice_hists/s_hist_mic_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_hist_mic_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_hists/s_hist_neu_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_hist_neu_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_hists/s_hist_sAb_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_hist_sAb_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_hists/s_hist_fAb_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_hist_fAb_d${dim}_s${sec[$i]}.gif
  ffmpeg -framerate 5 -start_number 10000 -i slice_hists/s_hist_ast_d${dim}_s${sec[$i]}_t%d.png -qscale 24 -r 5 -pix_fmt yuv420p gifs/s_hist_ast_d${dim}_s${sec[$i]}.gif

done
