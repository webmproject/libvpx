#!/bin/sh
# File:
#  encode_exp_bs.sh
# Decription:
#  This script fixes a video sequence and loops around an experimental list
#  defined by "exp_list"
#  see "encode_video_seq_bs.sh
set -x

if [ "$#" -ne 1 ]; then
  video_sequence=soccer_cif.sh
else
  video_sequence=$1
fi

root_dir=~/Dev/av1w
code_dir=$root_dir/aom
build_dir=$root_dir/release
test_dir=~/Dev/field
script_dir=~/Dev/sandbox/libvpx/scripts
bitstream_dir=~/Dev/samples/bitstreams/routine

. $script_dir/$video_sequence

videoname=$(basename $video_sequence .sh)

# General options
Codec=av1
codec="--codec=$Codec"
verbose=

cd $code_dir
#git checkout -q master
#git pull -q
commit_hash=`git log -1 --oneline | awk '{print $1}'`

cd $test_dir

d1="chroma_sub8x8 filter_7bit reference_buffer"
d2="delta_q rect_tx global_motion ext_tx"
d3="cdef ext_intra mv_compress ext_refs"
d4="dual_filter motion_var warped_motion"
d5="ext_delta_q loopfiltering_across_tiles ec_smallmul"
d6="var_tx ext_inter wedge compound_segment"
d7="interintra one_sided_compound smooth_hv"
d8="parallel_deblocking rect_intra_pred convolve_round"
d9="palette_throughput tempmv_signaling ext-comp-refs"

exp_list="$d1 $d2 $d3 $d4 $d5 $d6 $d7 $d8 $d9"
#exp_list=experimental

constant_cmdline_options="--skip=0 -p 2 --good --cpu-used=0 --lag-in-frames=25 --min-q=0 --max-q=63 --auto-alt-ref=1 --kf-max-dist=150 --kf-min-dist=0 --drop-frame=0 --static-thresh=0 --bias-pct=50 --minsection-pct=0 --maxsection-pct=2000 --arnr-maxframes=7 --arnr-strength=5 --sharpness=0 --undershoot-pct=100 --overshoot-pct=100 --frame-parallel=0 --test-decode=warn --psnr"

wi=1920
he=1080
frames=150
bitrate=4000

for exp_tool in $exp_list 

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

  profile=0
  bs="$Codec.$exp_tool.$videoname.$commit_hash.$profile.webm"
  
  ./aomenc $verbose -o $bitstream_dir/$bs $video $codec --limit=$frames --profile=$profile --fps=$fps $tune_content --target-bitrate=$bitrate --tile-columns=$col_num $constant_cmdline_options

  profile=2
  bs="$Codec.$exp_tool.$videoname.$commit_hash.$profile.webm"

  ./aomenc $verbose -o $bitstream_dir/$bs $video $codec --limit=$frames --profile=$profile --bit-depth=10 --fps=$fps $tune_content --target-bitrate=$bitrate --tile-columns=$col_num $constant_cmdline_options
 
  rm ./aomenc ./aomdec
done
