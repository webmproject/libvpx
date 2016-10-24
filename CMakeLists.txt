##
## Copyright (c) 2016, Alliance for Open Media. All rights reserved
##
## This source code is subject to the terms of the BSD 2 Clause License and
## the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
## was not distributed with this source code in the LICENSE file, you can
## obtain it at www.aomedia.org/license/software. If the Alliance for Open
## Media Patent License 1.0 was not distributed with this source code in the
## PATENTS file, you can obtain it at www.aomedia.org/license/patent.
##
cmake_minimum_required(VERSION 3.2)
project(AOM C CXX)

set(AOM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")

set(AOM_SRCS
    "${AOM_ROOT}/aom/aom.h"
    "${AOM_ROOT}/aom/aom_codec.h"
    "${AOM_ROOT}/aom/aom_decoder.h"
    "${AOM_ROOT}/aom/aom_encoder.h"
    "${AOM_ROOT}/aom/aom_frame_buffer.h"
    "${AOM_ROOT}/aom/aom_image.h"
    "${AOM_ROOT}/aom/aom_integer.h"
    "${AOM_ROOT}/aom/aomcx.h"
    "${AOM_ROOT}/aom/aomdx.h"
    "${AOM_ROOT}/aom/internal/aom_codec_internal.h"
    "${AOM_ROOT}/aom/src/aom_codec.c"
    "${AOM_ROOT}/aom/src/aom_decoder.c"
    "${AOM_ROOT}/aom/src/aom_encoder.c"
    "${AOM_ROOT}/aom/src/aom_image.c")

set(AOM_DSP_SRCS
    "${AOM_ROOT}/aom_dsp/aom_convolve.c"
    "${AOM_ROOT}/aom_dsp/aom_convolve.h"
    "${AOM_ROOT}/aom_dsp/aom_dsp_common.h"
    "${AOM_ROOT}/aom_dsp/aom_dsp_rtcd.c"
    "${AOM_ROOT}/aom_dsp/aom_filter.h"
    "${AOM_ROOT}/aom_dsp/aom_simd.c"
    "${AOM_ROOT}/aom_dsp/aom_simd.h"
    "${AOM_ROOT}/aom_dsp/aom_simd_inline.h"
    "${AOM_ROOT}/aom_dsp/avg.c"
    "${AOM_ROOT}/aom_dsp/bitreader.h"
    "${AOM_ROOT}/aom_dsp/bitreader_buffer.c"
    "${AOM_ROOT}/aom_dsp/bitreader_buffer.h"
    "${AOM_ROOT}/aom_dsp/bitwriter.h"
    "${AOM_ROOT}/aom_dsp/bitwriter_buffer.c"
    "${AOM_ROOT}/aom_dsp/bitwriter_buffer.h"
    "${AOM_ROOT}/aom_dsp/blend.h"
    "${AOM_ROOT}/aom_dsp/blend_a64_hmask.c"
    "${AOM_ROOT}/aom_dsp/blend_a64_mask.c"
    "${AOM_ROOT}/aom_dsp/blend_a64_vmask.c"
    "${AOM_ROOT}/aom_dsp/dkboolreader.c"
    "${AOM_ROOT}/aom_dsp/dkboolreader.h"
    "${AOM_ROOT}/aom_dsp/dkboolwriter.c"
    "${AOM_ROOT}/aom_dsp/dkboolwriter.h"
    "${AOM_ROOT}/aom_dsp/fwd_txfm.c"
    "${AOM_ROOT}/aom_dsp/fwd_txfm.h"
    "${AOM_ROOT}/aom_dsp/intrapred.c"
    "${AOM_ROOT}/aom_dsp/inv_txfm.c"
    "${AOM_ROOT}/aom_dsp/inv_txfm.h"
    "${AOM_ROOT}/aom_dsp/loopfilter.c"
    "${AOM_ROOT}/aom_dsp/prob.c"
    "${AOM_ROOT}/aom_dsp/prob.h"
    "${AOM_ROOT}/aom_dsp/psnr.c"
    "${AOM_ROOT}/aom_dsp/psnr.h"
    "${AOM_ROOT}/aom_dsp/quantize.c"
    "${AOM_ROOT}/aom_dsp/quantize.h"
    "${AOM_ROOT}/aom_dsp/sad.c"
    "${AOM_ROOT}/aom_dsp/simd/v128_intrinsics.h"
    "${AOM_ROOT}/aom_dsp/simd/v128_intrinsics_c.h"
    "${AOM_ROOT}/aom_dsp/simd/v256_intrinsics.h"
    "${AOM_ROOT}/aom_dsp/simd/v256_intrinsics_c.h"
    "${AOM_ROOT}/aom_dsp/simd/v64_intrinsics.h"
    "${AOM_ROOT}/aom_dsp/simd/v64_intrinsics_c.h"
    "${AOM_ROOT}/aom_dsp/subtract.c"
    "${AOM_ROOT}/aom_dsp/txfm_common.h"
    "${AOM_ROOT}/aom_dsp/variance.c"
    "${AOM_ROOT}/aom_dsp/variance.h")

set(AOM_MEM_SRCS
    "${AOM_ROOT}/aom_mem/aom_mem.c"
    "${AOM_ROOT}/aom_mem/aom_mem.h"
    "${AOM_ROOT}/aom_mem/include/aom_mem_intrnl.h")

set(AOM_SCALE_SRCS
    "${AOM_ROOT}/aom_scale/aom_scale.h"
    "${AOM_ROOT}/aom_scale/aom_scale_rtcd.c"
    "${AOM_ROOT}/aom_scale/generic/aom_scale.c"
    "${AOM_ROOT}/aom_scale/generic/gen_scalers.c"
    "${AOM_ROOT}/aom_scale/generic/yv12config.c"
    "${AOM_ROOT}/aom_scale/generic/yv12extend.c"
    "${AOM_ROOT}/aom_scale/yv12config.h")

# TODO(tomfinegan): Extract aom_ports from aom_util if possible.
set(AOM_UTIL_SRCS
    "${AOM_ROOT}/aom_ports/aom_once.h"
    "${AOM_ROOT}/aom_ports/aom_timer.h"
    "${AOM_ROOT}/aom_ports/bitops.h"
    "${AOM_ROOT}/aom_ports/emmintrin_compat.h"
    "${AOM_ROOT}/aom_ports/mem.h"
    "${AOM_ROOT}/aom_ports/mem_ops.h"
    "${AOM_ROOT}/aom_ports/mem_ops_aligned.h"
    "${AOM_ROOT}/aom_ports/msvc.h"
    "${AOM_ROOT}/aom_ports/system_state.h"
    "${AOM_ROOT}/aom_util/aom_thread.c"
    "${AOM_ROOT}/aom_util/aom_thread.h"
    "${AOM_ROOT}/aom_util/endian_inl.h")

set(AOM_AV1_COMMON_SRCS
    "${AOM_ROOT}/av1/av1_iface_common.h"
    "${AOM_ROOT}/av1/common/alloccommon.c"
    "${AOM_ROOT}/av1/common/alloccommon.h"
    "${AOM_ROOT}/av1/common/av1_fwd_txfm.c"
    "${AOM_ROOT}/av1/common/av1_fwd_txfm.h"
    "${AOM_ROOT}/av1/common/av1_inv_txfm.c"
    "${AOM_ROOT}/av1/common/av1_inv_txfm.h"
    "${AOM_ROOT}/av1/common/av1_rtcd.c"
    "${AOM_ROOT}/av1/common/blockd.c"
    "${AOM_ROOT}/av1/common/blockd.h"
    "${AOM_ROOT}/av1/common/common.h"
    "${AOM_ROOT}/av1/common/common_data.h"
    "${AOM_ROOT}/av1/common/convolve.c"
    "${AOM_ROOT}/av1/common/convolve.h"
    "${AOM_ROOT}/av1/common/debugmodes.c"
    "${AOM_ROOT}/av1/common/entropy.c"
    "${AOM_ROOT}/av1/common/entropy.h"
    "${AOM_ROOT}/av1/common/entropymode.c"
    "${AOM_ROOT}/av1/common/entropymode.h"
    "${AOM_ROOT}/av1/common/entropymv.c"
    "${AOM_ROOT}/av1/common/entropymv.h"
    "${AOM_ROOT}/av1/common/enums.h"
    "${AOM_ROOT}/av1/common/filter.c"
    "${AOM_ROOT}/av1/common/filter.h"
    "${AOM_ROOT}/av1/common/frame_buffers.c"
    "${AOM_ROOT}/av1/common/frame_buffers.h"
    "${AOM_ROOT}/av1/common/idct.c"
    "${AOM_ROOT}/av1/common/idct.h"
    "${AOM_ROOT}/av1/common/loopfilter.c"
    "${AOM_ROOT}/av1/common/loopfilter.h"
    "${AOM_ROOT}/av1/common/mv.h"
    "${AOM_ROOT}/av1/common/mvref_common.c"
    "${AOM_ROOT}/av1/common/mvref_common.h"
    "${AOM_ROOT}/av1/common/odintrin.c"
    "${AOM_ROOT}/av1/common/odintrin.h"
    "${AOM_ROOT}/av1/common/onyxc_int.h"
    "${AOM_ROOT}/av1/common/pred_common.c"
    "${AOM_ROOT}/av1/common/pred_common.h"
    "${AOM_ROOT}/av1/common/quant_common.c"
    "${AOM_ROOT}/av1/common/quant_common.h"
    "${AOM_ROOT}/av1/common/reconinter.c"
    "${AOM_ROOT}/av1/common/reconinter.h"
    "${AOM_ROOT}/av1/common/reconintra.c"
    "${AOM_ROOT}/av1/common/reconintra.h"
    "${AOM_ROOT}/av1/common/scale.c"
    "${AOM_ROOT}/av1/common/scale.h"
    "${AOM_ROOT}/av1/common/scan.c"
    "${AOM_ROOT}/av1/common/scan.h"
    "${AOM_ROOT}/av1/common/seg_common.c"
    "${AOM_ROOT}/av1/common/seg_common.h"
    "${AOM_ROOT}/av1/common/thread_common.c"
    "${AOM_ROOT}/av1/common/thread_common.h"
    "${AOM_ROOT}/av1/common/tile_common.c"
    "${AOM_ROOT}/av1/common/tile_common.h")

set(AOM_AV1_DECODER_SRCS
    "${AOM_ROOT}/av1/av1_dx_iface.c"
    "${AOM_ROOT}/av1/decoder/decodeframe.c"
    "${AOM_ROOT}/av1/decoder/decodeframe.h"
    "${AOM_ROOT}/av1/decoder/decodemv.c"
    "${AOM_ROOT}/av1/decoder/decodemv.h"
    "${AOM_ROOT}/av1/decoder/decoder.c"
    "${AOM_ROOT}/av1/decoder/decoder.h"
    "${AOM_ROOT}/av1/decoder/detokenize.c"
    "${AOM_ROOT}/av1/decoder/detokenize.h"
    "${AOM_ROOT}/av1/decoder/dsubexp.c"
    "${AOM_ROOT}/av1/decoder/dsubexp.h"
    "${AOM_ROOT}/av1/decoder/dthread.c"
    "${AOM_ROOT}/av1/decoder/dthread.h")

set(AOM_AV1_ENCODER_SRCS
    "${AOM_ROOT}/av1/av1_cx_iface.c"
    "${AOM_ROOT}/av1/encoder/aq_complexity.c"
    "${AOM_ROOT}/av1/encoder/aq_complexity.h"
    "${AOM_ROOT}/av1/encoder/aq_cyclicrefresh.c"
    "${AOM_ROOT}/av1/encoder/aq_cyclicrefresh.h"
    "${AOM_ROOT}/av1/encoder/aq_variance.c"
    "${AOM_ROOT}/av1/encoder/aq_variance.h"
    "${AOM_ROOT}/av1/encoder/bitstream.c"
    "${AOM_ROOT}/av1/encoder/bitstream.h"
    "${AOM_ROOT}/av1/encoder/block.h"
    "${AOM_ROOT}/av1/encoder/context_tree.c"
    "${AOM_ROOT}/av1/encoder/context_tree.h"
    "${AOM_ROOT}/av1/encoder/cost.c"
    "${AOM_ROOT}/av1/encoder/cost.h"
    "${AOM_ROOT}/av1/encoder/dct.c"
    "${AOM_ROOT}/av1/encoder/encodeframe.c"
    "${AOM_ROOT}/av1/encoder/encodeframe.h"
    "${AOM_ROOT}/av1/encoder/encodemb.c"
    "${AOM_ROOT}/av1/encoder/encodemb.h"
    "${AOM_ROOT}/av1/encoder/encodemv.c"
    "${AOM_ROOT}/av1/encoder/encodemv.h"
    "${AOM_ROOT}/av1/encoder/encoder.c"
    "${AOM_ROOT}/av1/encoder/encoder.h"
    "${AOM_ROOT}/av1/encoder/ethread.c"
    "${AOM_ROOT}/av1/encoder/ethread.h"
    "${AOM_ROOT}/av1/encoder/extend.c"
    "${AOM_ROOT}/av1/encoder/extend.h"
    "${AOM_ROOT}/av1/encoder/firstpass.c"
    "${AOM_ROOT}/av1/encoder/firstpass.h"
    "${AOM_ROOT}/av1/encoder/hybrid_fwd_txfm.c"
    "${AOM_ROOT}/av1/encoder/hybrid_fwd_txfm.h"
    "${AOM_ROOT}/av1/encoder/lookahead.c"
    "${AOM_ROOT}/av1/encoder/lookahead.h"
    "${AOM_ROOT}/av1/encoder/mbgraph.c"
    "${AOM_ROOT}/av1/encoder/mbgraph.h"
    "${AOM_ROOT}/av1/encoder/mcomp.c"
    "${AOM_ROOT}/av1/encoder/mcomp.h"
    "${AOM_ROOT}/av1/encoder/picklpf.c"
    "${AOM_ROOT}/av1/encoder/picklpf.h"
    "${AOM_ROOT}/av1/encoder/quantize.c"
    "${AOM_ROOT}/av1/encoder/quantize.h"
    "${AOM_ROOT}/av1/encoder/ratectrl.c"
    "${AOM_ROOT}/av1/encoder/ratectrl.h"
    "${AOM_ROOT}/av1/encoder/rd.c"
    "${AOM_ROOT}/av1/encoder/rd.h"
    "${AOM_ROOT}/av1/encoder/rdopt.c"
    "${AOM_ROOT}/av1/encoder/rdopt.h"
    "${AOM_ROOT}/av1/encoder/resize.c"
    "${AOM_ROOT}/av1/encoder/resize.h"
    "${AOM_ROOT}/av1/encoder/segmentation.c"
    "${AOM_ROOT}/av1/encoder/segmentation.h"
    "${AOM_ROOT}/av1/encoder/speed_features.c"
    "${AOM_ROOT}/av1/encoder/speed_features.h"
    "${AOM_ROOT}/av1/encoder/subexp.c"
    "${AOM_ROOT}/av1/encoder/subexp.h"
    "${AOM_ROOT}/av1/encoder/temporal_filter.c"
    "${AOM_ROOT}/av1/encoder/temporal_filter.h"
    "${AOM_ROOT}/av1/encoder/tokenize.c"
    "${AOM_ROOT}/av1/encoder/tokenize.h"
    "${AOM_ROOT}/av1/encoder/treewriter.c"
    "${AOM_ROOT}/av1/encoder/treewriter.h")

# Targets
add_library(aom_dsp ${AOM_DSP_SRCS})
add_library(aom_mem ${AOM_MEM_SRCS})
add_library(aom_scale ${AOM_SCALE_SRCS})
add_library(aom_util ${AOM_UTIL_SRCS})
add_library(aom_av1_decoder ${AOM_AV1_DECODER_SRCS})
add_library(aom_av1_encoder ${AOM_AV1_ENCODER_SRCS})
add_library(aom ${AOM_SRCS})
target_link_libraries(aom LINK_PUBLIC
                      aom_dsp
                      aom_mem
                      aom_scale
                      aom_util
                      aom_av1_decoder
                      aom_av1_encoder)
add_executable(simple_decoder examples/simple_decoder.c)
include_directories(${AOM_ROOT})
target_link_libraries(simple_decoder LINK_PUBLIC aom)
add_executable(simple_encoder examples/simple_encoder.c)
include_directories(${AOM_ROOT})
target_link_libraries(simple_encoder LINK_PUBLIC aom)

