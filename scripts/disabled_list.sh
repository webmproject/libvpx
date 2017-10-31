#!/bin/sh

d2="--disable-ext_tx"
d3="--disable-cdef --disable-ext_intra --disable-intra-edge --disable-mv_compress"
d4="--disable-dual_filter"
#--disable-motion_var --disable-warped_motion"
#d5="--disable-wedge --disable-compound_segment"
d6="--disable-interintra --disable-one_sided_compound"
d7="--disable-convolve_round --disable-aom-qm --disable-dist_8x8"
d8="--disable-loop_restoration"
d9="--disable-ext-partition --disable-ext-partition-types"
d91="--disable-new-multisymbol"
d10="--disable-reference_buffer --disable-loopfiltering_across_tiles"
d11="--disable-palette_throughput --disable-smooth_hv --disable-tempmv_signaling"
d12="--disable-ext-comp-refs --disable-ext_delta_q --disable-parallel_deblocking"
d13="--disable-simple_bwd_adapt --disable-loopfilter-level"

#d10="--disable-txmg"
# --disable-cdef_singlepass

disabled="$d2 $d3 $d4 $d6 $d7 $d8 $d9 $d10 $d11 $d12 $d13"
