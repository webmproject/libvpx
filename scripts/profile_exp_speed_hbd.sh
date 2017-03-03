#!/bin/sh
# File:
#  profile_exp_speed_hbd.sh
# Decription:
#  Configure, build, and run encoder/decoder perf for each experimental tool.
#  Save the decoder perf data
# Note:
#  - use perf 4.4.0
#  - see encoder config output if set, verbose=-v

root_dir=~/Dev/aomedia
code_dir=$root_dir/aom
build_dir=$root_dir/build
test_dir=~/Dev/field
script_dir=~/Dev/refine/libvpx/scripts

# video=~/Dev/samples/videos/yaowu/soccer_cif.y4m
# wi=352
# he=288
# frames=25
# bitrate=500
# fps="30/1"

# video=~/Dev/samples/videos/speed-set/touchdown_pass_480p.y4m
# wi=854
# he=480
# frames=100
# bitrate=1200
# fps="30000/1001"

video=~/Dev/samples/videos/hbd/crowd_run_1080p_10.y4m
wi=1920
he=1080
frames=300
bitrate=6000
fps="50/1"

codec="--codec=av1"
verbose=

# Note:
#  Standard bit depth:
#   1) profile=0
#   2) remove $bitdepth in encoder command line
#   3) Change runconfig.sh, bitdepth=
#  High bit depth:
#   1) profile=2
#   2) Add $bitdepth in encoder command line
#   3) Change runconfig.sh, bitdepth=--enable-aom-highbitdepth

profile=2
bitdepth="--bit-depth=10"

rm *.txt

#experimental var-tx ans entropy
#ext-intra filter-intra
#supertx ext-interp motion-var
#new-quant dual-filter ext-partition-types
#ext-partition ext-inter ext-refs
#ext-tx rect-tx global-motion
#loop-restoration alt-intra cb4x4

for exp_tool in supertx ext-interp motion-var new-quant dual-filter ext-partition-types ext-partition ext-inter ext-refs ext-tx rect-tx global-motion loop-restoration alt-intra cb4x4

do
  cd $build_dir
  make clean > /dev/null
  $script_dir/runconfig_hbd.sh $exp_tool
  make -j > /dev/null
  if [ $? -ne 0 ]; then
    echo Build failed on experiment: $exp_tool
  fi
  cp ./aomenc $test_dir/.
  cp ./aomdec $test_dir/.

  cd $test_dir
  elog=e_$exp_tool.txt
  dlog=d_$exp_tool.txt

  echo $exp_tool >> $elog
  bitstream=bs_$exp_tool.webm
  echo $exp_tool
  
  ./aomenc -o $bitstream $video $codec --limit=$frames --fps=$fps $verbose --profile=$profile $bitdepth --target-bitrate=$bitrate --skip=0 -p 2 --good --cpu-used=0 --lag-in-frames=25 --min-q=0 --max-q=63 --auto-alt-ref=1 --kf-max-dist=150 --kf-min-dist=0 --drop-frame=0 --static-thresh=0 --bias-pct=50 --minsection-pct=0 --maxsection-pct=2000 --arnr-maxframes=7 --arnr-strength=5 --sharpness=0 --undershoot-pct=100 --overshoot-pct=100 --frame-parallel=0 --tile-columns=0 --psnr
  
  #perf_file=dec_perf_$exp_tool.data
  
  #/usr/bin/perf record -e cycles -o $perf_file ./aomdec $bitstream -o too.yuv --summary

done



