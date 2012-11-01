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
#include "vpx_ports/arm.h"
#include "vp9/common/pragmas.h"
#include "vp9/common/subpixel.h"
#include "vp9/common/loopfilter.h"
#include "vp9/common/recon.h"
#include "vp9/common/idct.h"
#include "vp9/common/onyxc_int.h"

void vp9_arch_arm_common_init(VP9_COMMON *ctx) {
#if CONFIG_RUNTIME_CPU_DETECT
  VP9_COMMON_RTCD *rtcd = &ctx->rtcd;
  int flags = arm_cpu_caps();
  rtcd->flags = flags;

  /* Override default functions with fastest ones for this CPU. */
#if HAVE_ARMV5TE
  if (flags & HAS_EDSP) {
  }
#endif

// The commented functions need to be re-written for vpx.
#if HAVE_ARMV6
  if (flags & HAS_MEDIA) {
    rtcd->subpix.sixtap16x16   = vp9_sixtap_predict16x16_armv6;
    rtcd->subpix.sixtap8x8     = vp9_sixtap_predict8x8_armv6;
    rtcd->subpix.sixtap8x4     = vp9_sixtap_predict8x4_armv6;
    rtcd->subpix.sixtap4x4     = vp9_sixtap_predict_armv6;

    rtcd->subpix.bilinear16x16 = vp9_bilinear_predict16x16_armv6;
    rtcd->subpix.bilinear8x8   = vp9_bilinear_predict8x8_armv6;
    rtcd->subpix.bilinear8x4   = vp9_bilinear_predict8x4_armv6;
    rtcd->subpix.bilinear4x4   = vp9_bilinear_predict4x4_armv6;

    // rtcd->idct.idct1        = vp9_short_idct4x4llm_1_v6;
    // rtcd->idct.idct16       = vp9_short_idct4x4llm_v6_dual;
    // rtcd->idct.iwalsh1      = vp9_short_inv_walsh4x4_1_v6;
    // rtcd->idct.iwalsh16     = vp9_short_inv_walsh4x4_v6;

    rtcd->recon.copy16x16   = vp9_copy_mem16x16_v6;
    rtcd->recon.copy8x8     = vp9_copy_mem8x8_v6;
    rtcd->recon.copy8x4     = vp9_copy_mem8x4_v6;
    rtcd->recon.recon       = vp9_recon_b_armv6;
    rtcd->recon.recon2      = vp9_recon2b_armv6;
    rtcd->recon.recon4      = vp9_recon4b_armv6;
  }
#endif

#if HAVE_ARMV7
  if (flags & HAS_NEON) {
    rtcd->subpix.sixtap16x16   = vp9_sixtap_predict16x16_neon;
    rtcd->subpix.sixtap8x8     = vp9_sixtap_predict8x8_neon;
    rtcd->subpix.sixtap8x4     = vp9_sixtap_predict8x4_neon;
    rtcd->subpix.sixtap4x4     = vp9_sixtap_predict_neon;

    rtcd->subpix.bilinear16x16 = vp9_bilinear_predict16x16_neon;
    rtcd->subpix.bilinear8x8   = vp9_bilinear_predict8x8_neon;
    rtcd->subpix.bilinear8x4   = vp9_bilinear_predict8x4_neon;
    rtcd->subpix.bilinear4x4   = vp9_bilinear_predict4x4_neon;

    // rtcd->idct.idct1        = vp9_short_idct4x4llm_1_neon;
    // rtcd->idct.idct16       = vp9_short_idct4x4llm_neon;
    // rtcd->idct.iwalsh1      = vp9_short_inv_walsh4x4_1_neon;
    // rtcd->idct.iwalsh16     = vp9_short_inv_walsh4x4_neon;

    rtcd->recon.copy16x16   = vp9_copy_mem16x16_neon;
    rtcd->recon.copy8x8     = vp9_copy_mem8x8_neon;
    rtcd->recon.copy8x4     = vp9_copy_mem8x4_neon;
    rtcd->recon.recon       = vp9_recon_b_neon;
    rtcd->recon.recon2      = vp9_recon2b_neon;
    rtcd->recon.recon4      = vp9_recon4b_neon;
    rtcd->recon.recon_mb    = vp9_recon_mb_neon;
    rtcd->recon.build_intra_predictors_mby =
      vp9_build_intra_predictors_mby_neon;
    rtcd->recon.build_intra_predictors_mby_s =
      vp9_build_intra_predictors_mby_s_neon;
  }
#endif

#endif
}
