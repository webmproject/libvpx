#!/bin/sh

tool=$1
platform=x86_64-linux-gcc
codec=--enable-av1
libsrc=aom
test_dir=~/Dev/nightly

common="--disable-unit-tests --disable-docs --enable-experimental"

disabled="--disable-chroma_sub8x8 --disable-filter_7bit --disable-reference_buffer --disable-delta_q --disable-rect_tx --disable-global_motion --disable-ext_tx --disable-cdef --disable-ext_intra --disable-mv_compress --disable-ext_refs --disable-dual_filter --disable-motion_var --disable-warped_motion --disable-ext_delta_q --disable-loopfiltering_across_tiles --disable-ec_smallmul --disable-var_tx --disable-ext_inter --disable-wedge --disable-compound_segment --disable-interintra --disable-one_sided_compound --disable-smooth_hv --disable-parallel_deblocking --disable-rect_intra_pred --disable-convolve_round --disable-palette_throughput --disable-tempmv_signaling"

../$libsrc/configure $common $disabled --enable-$tool > /dev/null
if [ $? -ne 0 ]; then
  echo "Error: configure fails!" > $test_dir/error_config.txt
  exit 1
fi
