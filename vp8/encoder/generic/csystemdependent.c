/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include "vp8/encoder/variance.h"
#include "vp8/encoder/onyx_int.h"


void vp8_arch_x86_encoder_init(VP8_COMP *cpi);
void vp8_arch_arm_encoder_init(VP8_COMP *cpi);

void (*vp8_yv12_copy_partial_frame_ptr)(YV12_BUFFER_CONFIG *src_ybc, YV12_BUFFER_CONFIG *dst_ybc, int Fraction);
extern void vp8_yv12_copy_partial_frame(YV12_BUFFER_CONFIG *src_ybc, YV12_BUFFER_CONFIG *dst_ybc, int Fraction);

void vp8_cmachine_specific_config(VP8_COMP *cpi) {
#if CONFIG_RUNTIME_CPU_DETECT
  cpi->rtcd.common                    = &cpi->common.rtcd;

  cpi->rtcd.fdct.short8x8                  = vp8_short_fdct8x8_c;
  cpi->rtcd.fdct.short16x16                = vp8_short_fdct16x16_c;
  cpi->rtcd.fdct.haar_short2x2             = vp8_short_fhaar2x2_c;
  cpi->rtcd.fdct.short4x4                  = vp8_short_fdct4x4_c;
  cpi->rtcd.fdct.short8x4                  = vp8_short_fdct8x4_c;
  cpi->rtcd.fdct.fast4x4                   = vp8_short_fdct4x4_c;
  cpi->rtcd.fdct.fast8x4                   = vp8_short_fdct8x4_c;
  cpi->rtcd.fdct.walsh_short4x4            = vp8_short_walsh4x4_c;

  cpi->rtcd.encodemb.berr                  = vp8_block_error_c;
  cpi->rtcd.encodemb.mberr                 = vp8_mbblock_error_c;
  cpi->rtcd.encodemb.mbuverr               = vp8_mbuverror_c;
  cpi->rtcd.encodemb.subb                  = vp8_subtract_b_c;
  cpi->rtcd.encodemb.submby                = vp8_subtract_mby_c;
  cpi->rtcd.encodemb.submbuv               = vp8_subtract_mbuv_c;

  cpi->rtcd.search.full_search             = vp8_full_search_sad;
  cpi->rtcd.search.refining_search         = vp8_refining_search_sad;
  cpi->rtcd.search.diamond_search          = vp8_diamond_search_sad;
  cpi->rtcd.temporal.apply                 = vp8_temporal_filter_apply_c;
  cpi->rtcd.fdct.short4x4                  = vp8_short_fdct4x4_c;
  cpi->rtcd.fdct.short8x4                  = vp8_short_fdct8x4_c;
  cpi->rtcd.fdct.fast4x4                   = vp8_short_fdct4x4_c;
  cpi->rtcd.fdct.fast8x4                   = vp8_short_fdct8x4_c;
  cpi->rtcd.fdct.walsh_short4x4            = vp8_short_walsh4x4_c;
#endif

  vp8_yv12_copy_partial_frame_ptr = vp8_yv12_copy_partial_frame;

#if ARCH_X86 || ARCH_X86_64
  vp8_arch_x86_encoder_init(cpi);
#endif

#if ARCH_ARM
  vp8_arch_arm_encoder_init(cpi);
#endif


}
