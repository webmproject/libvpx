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
test_dir=~/Dev/field
script_dir=~/Dev/sandbox/libvpx/scripts

# video=~/Dev/samples/videos/yaowu/soccer_cif.y4m
# wi=352
# he=288
# frames=10
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
frames=150
bitrate=4000
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
#   2) Add $bitdepth in encoder command line, e.g. bitdepth="--bit-depth=10"
#   3) Change runconfig.sh, bitdepth=--enable-highbitdepth

profile=0
rm *.txt

# mismatch: sub8x8_mc, ncobmc, ext-partition
# ext_tile?
# enc segfault: coef_interleave, pvq
# dec error: pvq
# entropy_stats: not an experiment

#ans intrabc adapt_scan convolve_round lv_map txk_sel compound_round ext_tx filter_intra intra_interp q_adapt_probs subframe_prob_update palette_delta_encoding var_tx tpl_mv chroma_2x2 global_motion new_quant compound_segment motion_var tripred ref_adapt compound_singleref ext_inter supertx loop_restoration warped_motion tx64x64 masked_tx interintra wedge
# delta_q filter_7bit frame_size reference_buffer tempmv_signaling dependent_horztiles palette_throughput ext_delta_q frame_superres cdef ext_partition_types cfl daala_ec rawbits fp_mb_stats xiphrc parallel_deblocking loopfiltering_across_tiles parallel_deblocking_15tap ec_adapt tile_groups new_multisymbol ec_smallmul daala_dist mv_compress

# experimental ans intrabc adapt_scan convolve_round lv_map txk_sel compound_round ext_tx filter_intra intra_interp q_adapt_probs subframe_prob_update palette_delta_encoding var_tx tpl_mv chroma_2x2 global_motion new_quant compound_segment

#default:
for exp_tool in experimental ans intrabc adapt_scan convolve_round lv_map txk_sel compound_round ext_tx filter_intra intra_interp q_adapt_probs subframe_prob_update palette_delta_encoding var_tx tpl_mv chroma_2x2 global_motion new_quant compound_segment ext_partition_types tripred ref_adapt compound_singleref ext_inter supertx loop_restoration warped_motion tx64x64 masked_tx interintra wedge frame_size tempmv_signaling dependent_horztiles ext_delta_q frame_superres cfl daala_ec rawbits xiphrc parallel_deblocking loopfiltering_across_tiles parallel_deblocking_15tap new_multisymbol ec_smallmul daala_dist 

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

  elog=e_$exp_tool.txt
  dlog=d_$exp_tool.txt
  bs=$bs_$exp_tool.webm

  if [ $exp_tool == intrabc ] || [ $exp_tool == palette ] || [ $exp_tool == palette_delta_encoding ] || [ $exp_tool == palette_throughput ]; then
    tune_content="--tune-content=screen"
  else
    tune_content=
  fi
  
  taskset -c 3 ./aomenc $verbose -o /dev/shm/$bs $video $codec --limit=$frames --profile=$profile --fps=$fps $tune_content --target-bitrate=$bitrate --skip=0 -p 2 --good --cpu-used=0 --lag-in-frames=25 --min-q=0 --max-q=63 --auto-alt-ref=1 --kf-max-dist=150 --kf-min-dist=0 --drop-frame=0 --static-thresh=0 --bias-pct=50 --minsection-pct=0 --maxsection-pct=2000 --arnr-maxframes=7 --arnr-strength=5 --sharpness=0 --undershoot-pct=100 --overshoot-pct=100 --frame-parallel=0 --tile-columns=0 --test-decode=warn --psnr &>> $elog
  
  etime=`cat $elog | grep 'Pass 2/2' | grep 'fps)' | sed -e 's/^.*b\/s//' | awk '{print $1" "$2}'`
  psnr=`cat $elog | grep 'PSNR' | awk '{print $5, $6, $7, $8, $9}'`
  tmp=`cat $elog | grep mismatch`
  if [ "$?" -ne 0 ]; then
    eflag=e_ok
  else
    eflag=mismatch
  fi
  
  taskset -c 3 ./aomdec /dev/shm/$bs $codec --i420 --noblit --summary 2>&1 &>> $dlog
  if [ "$?" -ne 0 ]; then
    dflag=fault
  else
    dflag=d_ok
  fi

  dtime=`awk '{print $7" "$8}' < $dlog`
  
  echo -e $exp_tool '\t'$etime'\t'$dtime'\t'$psnr'\t'$eflag'\t'$dflag

  rm ./aomenc ./aomdec
done

