#!/bin/bash

# path to ffmpeg
FFMPEG=/usr/bin/ffmpeg
#FFMPEG=c:/usr/ffmpeg/ffmpeg.exe

# create video from png files
${FFMPEG} -r 30 -i "output_L_%03d.png" -an -vcodec libx264 -pix_fmt yuv420p left.mp4  
${FFMPEG} -r 30 -i "output_R_%03d.png" -an -vcodec libx264 -pix_fmt yuv420p right.mp4  
# concatenate video in forward and reverse directions
${FFMPEG} -i left.mp4 -filter_complex "[0:v]reverse,fifo[r];[0:v][r] concat=n=2:v=1 [v]" -map "[v]" left_reverse.mp4
${FFMPEG} -i right.mp4 -filter_complex "[0:v]reverse,fifo[r];[0:v][r] concat=n=2:v=1 [v]" -map "[v]" right_reverse.mp4
# looped
rm -f list.txt
for i in {1..8}; do printf "file '%s'\n" left_reverse.mp4 >> list.txt; done
${FFMPEG} -f concat -i list.txt -c copy left_looped.mp4
rm -f list.txt
for i in {1..8}; do printf "file '%s'\n" right_reverse.mp4 >> list.txt; done
${FFMPEG} -f concat -i list.txt -c copy right_looped.mp4
# side-by-side
${FFMPEG} -i left_looped.mp4 -vf "movie=right_looped.mp4 [in1]; [in]pad=iw*2:ih:iw:0[in0]; [in0][in1] overlay=0:0 [out]" sbs.mp4
