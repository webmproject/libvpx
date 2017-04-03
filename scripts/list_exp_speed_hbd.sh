#!/bin/sh
# File:
#  list_exp_speed.sh
# Decription:
#  Configure, build, and run encoder/decoder for each experimental tool.
#  Display the encoder/decode run time
# Note:
#  See encoder config output if set,
#  verbose=-v

root_dir=~/Dev/experiment
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

# video=~/Videos/aom/park_joy/park_joy_360p_10.y4m
# wi=640
# he=360
# frames=5
# bitrate=800
# fps="50/1"

# video=~/Dev/samples/videos/speed-set/touchdown_pass_480p.y4m
# wi=854
# he=480
# frames=100
# bitrate=1200
# fps="30000/1001"

video=~/Dev/samples/videos/hbd/crowd_run_1080p_10.y4m
wi=1920
he=1080
frames=100
bitrate=8000
fps="50/1"

bs=bs
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
#  experimental ans var-tx entropy ext-intra filter-intra supertx ext-interp motion-var new-quant dual-filter ext-partition-types ext-partition ext-inter ext-tx ext-refs rect-tx global-motion loop-restoration

#for exp_tool in experimental var_tx tpl_mv dual_filter convolve_round ext_tx tx64x64 sub8x8_mc ext_intra intra_interp filter_intra ext_inter interintra wedge compound_segment ext_refs global_motion new_quant supertx ans new_tokenset loop_restoration ext_partition ext_partition_types unpoison_partition_ctx ext_tile motion_var ncobmc warped_motion entropy q_adapt_probs subframe_prob_update bitstream_debug rawbits pvq xiphrc cb4x4 chroma_2x2 frame_size parallel_deblocking parallel_deblocking_15tap loopfiltering_across_tiles ec_adapt tempmv_signaling coef_interleave entropy_stats masked_tx dependent_horztiles daala_dist tripred palette_throughput ref_adapt lv_map mv_compress frame_superres new_multisymbol

# loop-restoration new_tokenset
for exp_tool in experimental global-motion cb4x4 tpl-mv sub8x8_mc interintra wedge compound-segment ncobmc warped_motion q_adapt_probs convolve_round var-tx ans entropy ext-intra filter-intra supertx motion-var new-quant dual-filter ext-partition-types ext-partition ext-inter ext-refs ext-tx

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
  bs=$bs_$exp_tool.webm
  
  ./aomenc -o $bs $video $codec --limit=$frames --skip=0 -p 2 --good --cpu-used=0 --target-bitrate=$bitrate --lag-in-frames=25 --min-q=0 --max-q=63 --auto-alt-ref=1 --kf-max-dist=150 --kf-min-dist=0 --drop-frame=0 --static-thresh=0 --bias-pct=50 --minsection-pct=0 --maxsection-pct=2000 --arnr-maxframes=7 --arnr-strength=5 --sharpness=0 --undershoot-pct=100 --overshoot-pct=100 --frame-parallel=0 --tile-columns=0 --profile=$profile $bitdepth --test-decode=warn $verbose --psnr &>> $elog
  
  etime=`cat $elog | grep 'Pass 2/2' | grep 'fps' | sed -e 's/^.*b\/s//' | awk '{print $1" " $2}'`
  psnr=`cat $elog | grep 'PSNR' | awk '{print $5, $6, $7, $8, $9}'`

  ./aomdec $bs $codec --i420 -o frm_%wx%h_$exp_tool.yuv %w $wi %h $he --summary 2>&1 &>> $dlog

  dtime=`awk '{print $7" "$8}' < $dlog`
  echo -e $exp_tool'\t'$etime '\t'$dtime'\t'$psnr

done



