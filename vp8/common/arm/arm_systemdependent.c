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
#include "vp8/common/pragmas.h"
#include "vp8/common/subpixel.h"
#include "vp8/common/loopfilter.h"
#include "vp8/common/recon.h"
#include "vp8/common/onyxc_int.h"

void vp8_arch_arm_common_init(VP8_COMMON *ctx)
{
#if CONFIG_RUNTIME_CPU_DETECT
    VP8_COMMON_RTCD *rtcd = &ctx->rtcd;
    int flags = arm_cpu_caps();
    rtcd->flags = flags;

    /* Override default functions with fastest ones for this CPU. */
#if HAVE_EDSP
    if (flags & HAS_EDSP)
    {
    }
#endif

#if HAVE_MEDIA
    if (flags & HAS_MEDIA)
    {
        rtcd->subpix.sixtap16x16   = vp8_sixtap_predict16x16_armv6;
        rtcd->subpix.sixtap8x8     = vp8_sixtap_predict8x8_armv6;
        rtcd->subpix.sixtap8x4     = vp8_sixtap_predict8x4_armv6;
        rtcd->subpix.sixtap4x4     = vp8_sixtap_predict_armv6;
        rtcd->subpix.bilinear16x16 = vp8_bilinear_predict16x16_armv6;
        rtcd->subpix.bilinear8x8   = vp8_bilinear_predict8x8_armv6;
        rtcd->subpix.bilinear8x4   = vp8_bilinear_predict8x4_armv6;
        rtcd->subpix.bilinear4x4   = vp8_bilinear_predict4x4_armv6;

        rtcd->recon.copy16x16   = vp8_copy_mem16x16_v6;
        rtcd->recon.copy8x8     = vp8_copy_mem8x8_v6;
        rtcd->recon.copy8x4     = vp8_copy_mem8x4_v6;
        rtcd->recon.intra4x4_predict = vp8_intra4x4_predict_armv6;
    }
#endif

#if HAVE_NEON
    if (flags & HAS_NEON)
    {
        rtcd->subpix.sixtap16x16   = vp8_sixtap_predict16x16_neon;
        rtcd->subpix.sixtap8x8     = vp8_sixtap_predict8x8_neon;
        rtcd->subpix.sixtap8x4     = vp8_sixtap_predict8x4_neon;
        rtcd->subpix.sixtap4x4     = vp8_sixtap_predict_neon;
        rtcd->subpix.bilinear16x16 = vp8_bilinear_predict16x16_neon;
        rtcd->subpix.bilinear8x8   = vp8_bilinear_predict8x8_neon;
        rtcd->subpix.bilinear8x4   = vp8_bilinear_predict8x4_neon;
        rtcd->subpix.bilinear4x4   = vp8_bilinear_predict4x4_neon;

        rtcd->recon.copy16x16   = vp8_copy_mem16x16_neon;
        rtcd->recon.copy8x8     = vp8_copy_mem8x8_neon;
        rtcd->recon.copy8x4     = vp8_copy_mem8x4_neon;
        rtcd->recon.build_intra_predictors_mby =
            vp8_build_intra_predictors_mby_neon;
        rtcd->recon.build_intra_predictors_mby_s =
            vp8_build_intra_predictors_mby_s_neon;
    }
#endif

#endif
}
