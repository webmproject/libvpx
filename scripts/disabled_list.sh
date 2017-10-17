#!/bin/sh

d1="--disable-chroma_sub8x8"
d2="--disable-rect_tx --disable-global_motion --disable-ext_tx"
d3="--disable-cdef --disable-ext_intra --disable-intra-edge --disable-mv_compress"
d4="--disable-dual_filter --disable-motion_var --disable-warped_motion"
d5="--disable-var_tx --disable-wedge --disable-compound_segment"
d6="--disable-interintra --disable-one_sided_compound"
d7="--disable-convolve_round --disable-aom-qm --disable-dist_8x8"
d8="--disable-loop_restoration"
d9="--disable-ext-partition --disable-ext-partition-types"
d10="--disable-txmg"
d11="--disable-reference_buffer --disable-loopfiltering_across_tiles"
d12="--disable-palette_throughput --disable-smooth_hv --disable-tempmv_signaling"
d13="--disable-ext-comp-refs --disable-ext_delta_q --disable-parallel_deblocking --disable-simple_bwd_adapt"

# --disable-cdef_singlepass
disabled="$d1 $d2 $d3 $d4 $d5 $d6 $d7 $d8 $d9 $d10 $d11 $d12 $d13"
