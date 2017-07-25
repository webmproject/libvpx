#!/bin/sh
# File:
#  list_exp_speed.sh
# Decription:
#  Configure, build, and run encoder/decoder for each experimental tool.
#  Display the encoder/decode run time
# Preassumption:
#  1) Assume all script files are in ~/Dev/sandbox/libvpx/scripts
# Note:
#  See encoder config output if set,
#  verbose=-v
#set -x

# LBD or HBD
# Note:
#  Standard bit depth:
#   1) profile=0
#   2) bitdepth=
#  High bit depth:
#   1) profile=2
#   2) bitdepth="--bit-depth=10/12"

if [ "$#" -ne 2 ]; then
  root_dir=~/Dev/vp9d
  profile=0
else
  root_dir=$1
  profile=$2
fi

if [ "$profile" == "2" ]; then
  bitdepth="--bit-depth=10"
fi

if [ "$profile" == "0" ]; then
  bitdepth=
fi

code_dir=$root_dir/libvpx
build_dir=$root_dir/release
test_dir=~/Dev/nightly
script_dir=~/Dev/sandbox/libvpx/scripts

. $script_dir/video_sequence.sh

# General options
codec="--codec=vp9"
verbose=
core_id=1

cd $test_dir

bstream=vp9_profile_$profile.webm
elog=vp9enc_log_p_$profile.txt
dlog=vp9dec_log_p_$profile.txt

taskset -c $core_id ./vpxenc $verbose -o /dev/shm/"$bstream" $video $codec --limit=$frames --profile=$profile $bitdepth --fps=$fps --target-bitrate=$bitrate --skip=0 -p 2 --good --cpu-used=0 --lag-in-frames=25 --min-q=0 --max-q=63 --auto-alt-ref=1 --kf-max-dist=150 --kf-min-dist=0 --drop-frame=0 --static-thresh=0 --bias-pct=50 --minsection-pct=0 --maxsection-pct=2000 --arnr-maxframes=7 --arnr-strength=5 --sharpness=0 --undershoot-pct=100 --overshoot-pct=100 --frame-parallel=0 --test-decode=warn --psnr &>> $elog

# Note: $2 is the time unit, ms or us
#etime=`cat $elog | grep 'Pass 2/2' | grep 'fps)' | sed -e 's/^.*b\/s//' | awk '{print $1" "$2}'`
efps=`cat $elog | grep 'Pass 2/2' | grep 'fps)' | sed -e 's/^.*b\/s//' | awk '{print $3}'`
efps=`echo $efps | sed 's/(//'`

psnr=`cat $elog | grep 'PSNR' | awk '{print $5, $6, $7, $8, $9}'`
tmp=`cat $elog | grep mismatch`
if [ "$?" -ne 0 ]; then
  eflag=e_ok
else
  eflag=mismatch
fi

echo "VP9: $(basename $video), bitrate=$bitrate profile=$profile frames=$frames"

taskset -c $core_id ./vpxdec /dev/shm/"$bstream" $codec --i420 --noblit --summary 2>&1 &>> $dlog
if [ "$?" -ne 0 ]; then
  dflag=fault
else
  dflag=d_ok
fi

# Note: $8 is the time unit ms or us
#dtime=`awk '{print $7" "$8}' < $dlog`
dfps=`awk '{print $9}' < $dlog`
dfps=`echo $dfps | sed 's/(//'`

echo -e '\t'"Enc fps   Dec fps    PSNR"'\t\t\t\t\t\t\t'"Enc status   Dec status"
echo -e '\t'$efps"        "$dfps"     "$psnr'\t'$eflag"            "$dflag
printf "\n"

