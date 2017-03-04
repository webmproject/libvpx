#!/bin/sh
# File:
#  list_exp_speed.sh
# Decription:
#  Configure, build, and run encoder/decoder for each experimental tool.
#  Display the encoder/decode run time
# Note:
#  See encoder config output if set,
#  verbose=-v

root_dir=~/Dev/aomedia
code_dir=$root_dir/aom
build_dir=$root_dir/release
test_dir=~/Dev/field
script_dir=~/Dev/sandbox/libvpx/scripts

# video=~/Dev/samples/videos/yaowu/soccer_cif.y4m
# wi=352
# he=288
# frames=5
# bitrate=500
# fps="30/1"

# video=~/Dev/samples/videos/speed-set/touchdown_pass_480p.y4m
# wi=854
# he=480
# frames=150
# bitrate=2400
# fps="30000/1001"

video=~/Dev/samples/videos/speed-set/BasketballDrive_1920x1080_50.y4m
wi=1920
he=1080
frames=300
bitrate=5000
fps="50/1"

# General options
bs=bs
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
#   2) Add $bitdepth in encoder command line
#   3) Change runconfig.sh, bitdepth=--enable-aom-highbitdepth

profile=0
bitdepth="--bit-depth=10"

# Not independt: bidir-pred rect-tx
#  experimental var-tx ans entropy ext-intra filter-intra supertx ext-interp motion-var new-quant dual-filter ext-partition-types ext-partition ext-inter ext-refs ext-tx rect-tx global-motion loop-restoration alt-intra cb4x4

rm *.txt

for exp_tool in experimental var-tx ans entropy ext-intra filter-intra supertx ext-interp motion-var new-quant dual-filter ext-partition-types ext-partition ext-inter ext-refs ext-tx rect-tx global-motion loop-restoration alt-intra cb4x4
do
  cd $build_dir
  make clean > /dev/null
  $script_dir/runconfig.sh $exp_tool
  make -j > /dev/null
  if [ $? -ne 0 ]; then
    echo Build failed on experiment: $exp_tool
  fi
  cp ./aomenc $test_dir/.
  cp ./aomdec $test_dir/.

  cd $test_dir
  #rm $elog $dlog
  elog=e_$exp_tool.txt
  dlog=d_$exp_tool.txt
  bs=$bs_$exp_tool.webm
  
  ./aomenc $verbose -o $bs $video $codec --limit=$frames --profile=$profile --fps=$fps --skip=0 -p 2 --good --cpu-used=0 --target-bitrate=$bitrate --lag-in-frames=25 --min-q=0 --max-q=63 --auto-alt-ref=1 --kf-max-dist=150 --kf-min-dist=0 --drop-frame=0 --static-thresh=0 --bias-pct=50 --minsection-pct=0 --maxsection-pct=2000 --arnr-maxframes=7 --arnr-strength=5 --sharpness=0 --undershoot-pct=100 --overshoot-pct=100 --frame-parallel=0 --tile-columns=0 --test-decode=warn --psnr &>> $elog
  
  etime=`cat $elog | grep 'Pass 2/2' | grep 'fps' | sed -e 's/^.*b\/s//' | awk '{print $1" "$2}'`
  psnr=`cat $elog | grep 'PSNR' | awk '{print $5, $6, $7, $8, $9}'`

  ./aomdec $bs $codec --i420 -o frm_%wx%h_$exp_tool.yuv %w $wi %h $he --summary 2>&1 &>> $dlog

  dtime=`awk '{print $7" "$8}' < $dlog`
  
  echo -e $exp_tool '\t'$etime'\t'$dtime'\t'$psnr
done

