#!/bin/sh

platform=x86_64-linux-gcc
codec=--enable-av1
libsrc=aom

experiment_dependency() {
  tool=--enable-$1

  if [ "$1" == ec_adapt ]; then
    tool="--enable-$1 --enable-ext-inter"
  fi

  if [ "$1" == ec_multisymbol ]; then
    tool="--enable-$1 --enable-daala-ec"
  fi

  if [ "$1" == rect_tx ]; then
    tool="--enable-$1 --enable-ext-tx"
  fi

  if [ "$1" == simp_mv_pred ]; then
    tool="--enable-$1 --enable-ref-mv"
  fi
  
  if [ "$1" == mv_compress ]; then
    tool="--enable-$1 --enable-ref-mv"
  fi

  if [ "$1" == var_tx ]; then
    tool="--enable-$1 --enable-cb4x4 --enable-rect-tx"
  fi
  
  # if [ "$1" == ext_tx ]; then
  #   tool="--enable-$1 --enable-cb4x4"
  # fi
  
  # if [ "$1" == chroma_2x2 ]; then
  #   tool="--enable-$1 --enable-cb4x4"
  # fi

  if [ "$1" == ncobmc ]; then
    tool="--enable-$1 --enable-motion-var"
  fi
  
  if [ "$1" == wedge ] || [ "$1" == compound_segment ] ; then
    tool="--enable-$1 --enable-ext-inter"
  fi

  if [ "$1" == interintra ] ; then
    tool="--enable-$1 --enable-ext-inter --enable-wedge"
  fi

  if [ "$1" == txk_sel ] ; then
    tool="--enable-$1 --enable-lv_map"
  fi

  if [ "$1" == intra_interp ] ; then
    tool="--enable-$1 --enable-ext-intra"
  fi

  if [ "$1" == ext_tile ] ; then
    tool="--enable-$1 --disable-cdef"
  fi

  if [ "$1" == lv_map ] ; then
    tool="--enable-$1 --enable-cb4x4 --enable-chroma_sub8x8"
  fi

  # if [ "$1" == chroma_sub8x8 ] ; then
  #   tool="--enable-$1 --enable-cb4x4"
  # fi

  if [ "$1" == one_sided_compound ] ; then
    tool="--enable-$1 --enable-ext-refs"
  fi
  
}

common="--disable-unit-tests --disable-docs --enable-experimental"

d1="--disable-chroma_sub8x8 --disable-filter_7bit --disable-reference_buffer"
d2="--disable-delta_q --disable-rect_tx --disable-global_motion --disable-ext_tx"
d3="--disable-cdef --disable-ext_intra --disable-mv_compress --disable-ext_refs"
d4="--disable-dual_filter --disable-motion_var --disable-warped_motion"
d5="--disable-ext_delta_q --disable-loopfiltering_across_tiles --disable-ec_smallmul"
d6="--disable-var_tx --disable-ext_inter --disable-wedge --disable-compound_segment"
d7="--disable-interintra --disable-one_sided_compound --disable-smooth_hv"
d8="--disable-parallel_deblocking --disable-rect_intra_pred --disable-convolve_round"
d9="--disable-palette_throughput --disable-tempmv_signaling --disable-ext-comp-refs"

disabled="$d1 $d2 $d3 $d4 $d5 $d6 $d7 $d8 $d9"

experiment_dependency $1

../$libsrc/configure $common $disabled $tool > /dev/null
