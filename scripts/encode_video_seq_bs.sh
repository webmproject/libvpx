#!/bin/sh
# File:
#  encode_exp_bs.sh
# Decription:
#  Configure, build, and run encoder/decoder for each experimental tool.
#  Display the encoder/decode run time
# Preassumption:
#  1) Assume all script files are in ~/Dev/sandbox/libvpx/scripts
# Note:
#  See encoder config output if set,
#  verbose=-v
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
set -x

# if [ "$#" -ne 1 ]; then
#   video_sequence=soccer_cif.sh
# else
#   video_sequence=$1
# fi

root_dir=~/Dev/av1w
code_dir=$root_dir/aom
build_dir=$root_dir/release
test_dir=~/Dev/field
script_dir=~/Dev/sandbox/libvpx/scripts
bitstream_dir=~/Dev/samples/bitstreams/routine

#. $script_dir/$video_sequence

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

# Note: experimental list
#exp_list="$d1 $d2 $d3 $d4 $d5 $d6 $d7 $d8 $d9"
#exp_list=experimental
exp_tool=experimental

# Note: video sequence list
#video_sequence_list="BQTerrace_1080p60.sh BasketballDrive_1080p50.sh ParkScene_1080p24.sh"
video_sequence_list="blue_sky_1080p25.sh rush_hour_1080p25.sh tennis_1080p24.sh"

constant_cmdline_options="--skip=0 -p 2 --good --cpu-used=0 --lag-in-frames=25 --min-q=0 --max-q=63 --auto-alt-ref=1 --kf-max-dist=150 --kf-min-dist=0 --drop-frame=0 --static-thresh=0 --bias-pct=50 --minsection-pct=0 --maxsection-pct=2000 --arnr-maxframes=7 --arnr-strength=5 --sharpness=0 --undershoot-pct=100 --overshoot-pct=100 --frame-parallel=0 --test-decode=warn --psnr"
col_num=1

wi=1920
he=1080
frames=150
bitrate=4000

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

for video_sequence in $video_sequence_list
do

  . $script_dir/$video_sequence
  videoname=$(basename $video_sequence .sh)
  
  profile=0
  bs="$Codec.$exp_tool.$videoname.$commit_hash.$profile.webm"
  
  ./aomenc $verbose -o $bitstream_dir/$bs $video $codec --limit=$frames --profile=$profile --fps=$fps $tune_content --target-bitrate=$bitrate --tile-columns=$col_num $constant_cmdline_options

  # profile=2
  # bs="$Codec.$exp_tool.$videoname.$commit_hash.$profile.webm"

  # ./aomenc $verbose -o $bitstream_dir/$bs $video $codec --limit=$frames --profile=$profile --bit-depth=10 --fps=$fps $tune_content --target-bitrate=$bitrate --tile-columns=$col_num $constant_cmdline_options
 
done
