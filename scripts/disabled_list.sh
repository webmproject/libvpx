#!/bin/sh

d1="--disable-chroma_sub8x8 --disable-filter_7bit --disable-reference_buffer"
d2="--disable-delta_q --disable-rect_tx --disable-global_motion --disable-ext_tx"
d3="--disable-cdef --disable-ext_intra --disable-mv_compress --disable-ext_refs"
d4="--disable-dual_filter --disable-motion_var --disable-warped_motion"
d5="--disable-ext_delta_q --disable-loopfiltering_across_tiles"
d6="--disable-var_tx --disable-ext_inter --disable-wedge --disable-compound_segment"
d7="--disable-interintra --disable-one_sided_compound  --disable-ext-comp-refs"
d8="--disable-smooth_hv --disable-parallel_deblocking"
d9="--disable-convolve_round --disable-aom-qm --disable-dist_8x8"
d10=" --disable-palette_throughput --disable-tempmv_signaling"

disabled="$d1 $d2 $d3 $d4 $d5 $d6 $d7 $d8 $d9 $d10"
