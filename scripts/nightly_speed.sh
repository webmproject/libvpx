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

if [ "$#" -ne 1 ]; then
  root_dir=~/Dev/av1k
else
  root_dir=$1
fi

code_dir=$root_dir/aom
build_dir=$root_dir/release
test_dir=~/Dev/nightly
script_dir=~/Dev/sandbox/libvpx/scripts

. video_sequence.sh

# General options
codec="--codec=av1"
verbose=

# LBD or HBD
# Note:
#  Standard bit depth:
#   1) profile=0
#   2) remove $bitdepth in encoder command line
#   3) Change runconfig.sh, bitdepth=
#  High bit depth:
#   1) profile=2
#   2) Add $bitdepth in encoder command line, e.g. bitdepth="--bit-depth=10"
#   3) Change runconfig.sh, bitdepth=--enable-highbitdepth

profile=0

core_id=1

for exp_tool in experimental

do
  cd $code_dir
  git checkout -q master
  git pull -q
  git log -1 --oneline
  
  cd $build_dir
  make clean > /dev/null
  $script_dir/nightly_config.sh $exp_tool
  make -j > /dev/null
  if [ $? -ne 0 ]; then
    echo Build failed on experiment: $exp_tool
  fi

  cp -f ./aomenc $test_dir/.
  cp -f ./aomdec $test_dir/.

  cd $test_dir
  rm *.txt
  
  elog=e_$exp_tool.txt
  dlog=d_$exp_tool.txt
  bstream="$exp_tool"_nightly_av1.webm
  
  if [ $exp_tool == intrabc ] || [ $exp_tool == palette ] || [ $exp_tool == palette_delta_encoding ] || [ $exp_tool == palette_throughput ]; then
    tune_content="--tune-content=screen"
  else
    tune_content=
  fi

  if [ $exp_tool == ext_tile ]; then
    col_num=1
  else
    col_num=0
  fi
  
  taskset -c $core_id ./aomenc $verbose -o /dev/shm/"$bstream" $video $codec --limit=$frames --profile=$profile --fps=$fps $tune_content --target-bitrate=$bitrate --skip=0 -p 2 --good --cpu-used=0 --lag-in-frames=25 --min-q=0 --max-q=63 --auto-alt-ref=1 --kf-max-dist=150 --kf-min-dist=0 --drop-frame=0 --static-thresh=0 --bias-pct=50 --minsection-pct=0 --maxsection-pct=2000 --arnr-maxframes=7 --arnr-strength=5 --sharpness=0 --undershoot-pct=100 --overshoot-pct=100 --frame-parallel=0 --tile-columns=$col_num --test-decode=warn --psnr &>> $elog

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

  echo "AV1 bitstream: " "$bstream"
  
  taskset -c $core_id ./aomdec /dev/shm/"$bstream" $codec --i420 --noblit --summary 2>&1 &>> $dlog
  if [ "$?" -ne 0 ]; then
    dflag=fault
  else
    dflag=d_ok
  fi

  # Note: $8 is the time unit ms or us
  dfps=`awk '{print $9}' < $dlog`
  dfps=`echo $dfps | sed 's/(//'`

  echo -e '\t'"Enc fps   Dec fps    PSNR"'\t\t\t\t\t\t\t'"Enc status   Dec status"
  echo -e '\t'$efps"        "$dfps"     "$psnr'\t'$eflag"            "$dflag
  printf "\n"
done

