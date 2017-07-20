#!/bin/sh

platform=x86_64-linux-gcc
codec=--enable-av1
libsrc=aom

common="--disable-unit-tests --disable-docs --enable-experimental"

disabled="--disable-chroma_sub8x8 --disable-filter_7bit --disable-reference_buffer
 --disable-delta_q --disable-tile_groups --disable-rect_tx --disable-global_motion --disable-ext_tx --disable-cdef --disable-ext_intra --disable-mv_compress --disable-ext_refs --disable-dual_filter --disable-motion_var --disable-warped_motion --disable-ext_delta_q --disable-loopfiltering_across_tiles --disable-ec_smallmul --disable-var_tx --disable-ext_inter --disable-wedge --disable-compound_segment --disable-interintra --disable-one_sided_compound --disable-smooth_hv --disable-parallel_deblocking --disable-palette --disable-alt_intra --disable-palette_throughput --disable-tempmv_signaling"

../$libsrc/configure $common $disabled $tool > /dev/null
