#!/bin/sh
# Unique calling argument: frame number
#
set -x

if [ "$#" -ne 1 ]; then
  frames=25
else
  frames=$1
fi

prof=0
if [ "$prof" == "2" ]; then
  bitdepth="--bit-depth=10"
fi

if [ "$prof" == "0" ]; then
  bitdepth=
fi

# In a release direcotry
cd ../aom
#git pull
commithash=`git log --pretty=%h -1`
cd ../release

#../aom/configure --disable-docs --disable-unit-tests
#make clean;make

rm -fr *
cmake ../aom -DCONFIG_UNIT_TESTS=0 -DENABLE_DOCS=0
make

date_str=`date +%b_%d_%Y`
bitstream=night_720p30_av1_$commithash.$date_str.f$frames.webm

fps=30
bitrate=2500
laginframes=19

taskset -c 0 ./aomenc -o /run/shm/$bitstream ~/Dev/samples/videos/speed-set/night_720p30.y4m --codec=av1 --fps=$fps/1 --skip=0 -p 2 --good --cpu-used=0 --target-bitrate=$bitrate --lag-in-frames=$laginframes --profile=$prof $bitdepth --limit=$frames --enable-cdef=0 --min-q=0 --max-q=63 --auto-alt-ref=1 --kf-max-dist=150 --kf-min-dist=0 --drop-frame=0 --static-thresh=0 --bias-pct=50 --minsection-pct=0 --maxsection-pct=2000 --arnr-maxframes=7 --arnr-strength=5 --sharpness=0 --undershoot-pct=100 --overshoot-pct=100 --frame-parallel=0 -t 1 --psnr --test-decode=warn -D
