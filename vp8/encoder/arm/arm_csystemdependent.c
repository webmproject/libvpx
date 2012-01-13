/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_config.h"
#include "vpx_ports/arm.h"
#include "vp8/encoder/variance.h"
#include "vp8/encoder/onyx_int.h"

extern void (*vp8_yv12_copy_partial_frame_ptr)(YV12_BUFFER_CONFIG *src_ybc, YV12_BUFFER_CONFIG *dst_ybc);
extern void vp8_yv12_copy_partial_frame(YV12_BUFFER_CONFIG *src_ybc, YV12_BUFFER_CONFIG *dst_ybc);
extern void vp8_yv12_copy_partial_frame_neon(YV12_BUFFER_CONFIG *src_ybc, YV12_BUFFER_CONFIG *dst_ybc);

void vp8_arch_arm_encoder_init(VP8_COMP *cpi)
{
#if CONFIG_RUNTIME_CPU_DETECT
    int flags = cpi->common.cpu_caps;

#if HAVE_EDSP
    if (flags & HAS_EDSP)
    {
    }
#endif

#endif /* CONFIG_RUNTIME_CPU_DETECT */

#if HAVE_NEON
#if CONFIG_RUNTIME_CPU_DETECT
    if (flags & HAS_NEON)
#endif
    {
        vp8_yv12_copy_partial_frame_ptr = vp8_yv12_copy_partial_frame_neon;
    }
#endif
}
