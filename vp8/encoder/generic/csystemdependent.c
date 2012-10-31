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


void vp9_arch_x86_encoder_init(VP9_COMP *cpi);
void vp9_arch_arm_encoder_init(VP9_COMP *cpi);

void (*vp9_yv12_copy_partial_frame_ptr)(YV12_BUFFER_CONFIG *src_ybc,
                                        YV12_BUFFER_CONFIG *dst_ybc,
                                        int fraction);
extern void vp9_yv12_copy_partial_frame(YV12_BUFFER_CONFIG *src_ybc,
                                        YV12_BUFFER_CONFIG *dst_ybc,
                                        int fraction);

void vp9_cmachine_specific_config(VP9_COMP *cpi) {
#if CONFIG_RUNTIME_CPU_DETECT
  cpi->rtcd.common                    = &cpi->common.rtcd;

  cpi->rtcd.search.full_search             = vp9_full_search_sad;
  cpi->rtcd.search.refining_search         = vp9_refining_search_sad;
  cpi->rtcd.search.diamond_search          = vp9_diamond_search_sad;
  cpi->rtcd.temporal.apply                 = vp9_temporal_filter_apply_c;
#endif

  vp9_yv12_copy_partial_frame_ptr = vp9_yv12_copy_partial_frame;

#if ARCH_X86 || ARCH_X86_64
  vp9_arch_x86_encoder_init(cpi);
#endif

#if ARCH_ARM
  vp9_arch_arm_encoder_init(cpi);
#endif


}
