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
#include "vpx_rtcd.h"
#include "vp9/common/subpixel.h"
#include "vp9/common/loopfilter.h"
#include "vp9/common/idct.h"
#include "vp9/common/onyxc_int.h"

extern void vp9_arch_x86_common_init(VP9_COMMON *ctx);
extern void vp9_arch_arm_common_init(VP9_COMMON *ctx);

void vp9_machine_specific_config(VP9_COMMON *ctx) {
#if CONFIG_RUNTIME_CPU_DETECT
  VP9_COMMON_RTCD *rtcd = &ctx->rtcd;

  rtcd->idct.idct1        = vp9_short_idct4x4llm_1_c;
  rtcd->idct.idct16       = vp9_short_idct4x4llm_c;
  rtcd->idct.idct1_scalar_add = vp9_dc_only_idct_add_c;
  rtcd->idct.iwalsh1      = vp9_short_inv_walsh4x4_1_c;
  rtcd->idct.iwalsh16     = vp9_short_inv_walsh4x4_c;
  rtcd->idct.idct8        = vp9_short_idct8x8_c;
  rtcd->idct.idct1_scalar_add_8x8 = vp9_dc_only_idct_add_8x8_c;
  rtcd->idct.ihaar2       = vp9_short_ihaar2x2_c;
  rtcd->idct.idct16x16    = vp9_short_idct16x16_c;

  rtcd->subpix.eighttap16x16       = vp9_eighttap_predict16x16_c;
  rtcd->subpix.eighttap8x8         = vp9_eighttap_predict8x8_c;
  rtcd->subpix.eighttap_avg16x16   = vp9_eighttap_predict_avg16x16_c;
  rtcd->subpix.eighttap_avg8x8     = vp9_eighttap_predict_avg8x8_c;
  rtcd->subpix.eighttap_avg4x4     = vp9_eighttap_predict_avg4x4_c;
  rtcd->subpix.eighttap8x4         = vp9_eighttap_predict8x4_c;
  rtcd->subpix.eighttap4x4         = vp9_eighttap_predict_c;
  rtcd->subpix.eighttap16x16_sharp     = vp9_eighttap_predict16x16_sharp_c;
  rtcd->subpix.eighttap8x8_sharp       = vp9_eighttap_predict8x8_sharp_c;
  rtcd->subpix.eighttap_avg16x16_sharp = vp9_eighttap_predict_avg16x16_sharp_c;
  rtcd->subpix.eighttap_avg8x8_sharp   = vp9_eighttap_predict_avg8x8_sharp_c;
  rtcd->subpix.eighttap_avg4x4_sharp   = vp9_eighttap_predict_avg4x4_sharp_c;
  rtcd->subpix.eighttap8x4_sharp       = vp9_eighttap_predict8x4_sharp_c;
  rtcd->subpix.eighttap4x4_sharp       = vp9_eighttap_predict_sharp_c;

  rtcd->subpix.sixtap16x16       = vp9_sixtap_predict16x16_c;
  rtcd->subpix.sixtap8x8         = vp9_sixtap_predict8x8_c;
  rtcd->subpix.sixtap_avg16x16   = vp9_sixtap_predict_avg16x16_c;
  rtcd->subpix.sixtap_avg8x8     = vp9_sixtap_predict_avg8x8_c;
  rtcd->subpix.sixtap8x4         = vp9_sixtap_predict8x4_c;
  rtcd->subpix.sixtap4x4         = vp9_sixtap_predict_c;
  rtcd->subpix.sixtap_avg4x4     = vp9_sixtap_predict_avg_c;
  rtcd->subpix.bilinear16x16     = vp9_bilinear_predict16x16_c;
  rtcd->subpix.bilinear8x8       = vp9_bilinear_predict8x8_c;
  rtcd->subpix.bilinear_avg16x16 = vp9_bilinear_predict_avg16x16_c;
  rtcd->subpix.bilinear_avg8x8   = vp9_bilinear_predict_avg8x8_c;
  rtcd->subpix.bilinear8x4       = vp9_bilinear_predict8x4_c;
  rtcd->subpix.bilinear4x4       = vp9_bilinear_predict4x4_c;
  rtcd->subpix.bilinear_avg4x4   = vp9_bilinear_predict_avg4x4_c;

#if CONFIG_POSTPROC || (CONFIG_VP9_ENCODER && CONFIG_INTERNAL_STATS)
  rtcd->postproc.down             = vp9_mbpost_proc_down_c;
  rtcd->postproc.across           = vp9_mbpost_proc_across_ip_c;
  rtcd->postproc.downacross       = vp9_post_proc_down_and_across_c;
  rtcd->postproc.addnoise         = vp9_plane_add_noise_c;
  rtcd->postproc.blend_mb_inner   = vp9_blend_mb_inner_c;
  rtcd->postproc.blend_mb_outer   = vp9_blend_mb_outer_c;
  rtcd->postproc.blend_b          = vp9_blend_b_c;
#endif

#endif

#if ARCH_X86 || ARCH_X86_64
  vp9_arch_x86_common_init(ctx);
#endif

#if ARCH_ARM
  vp9_arch_arm_common_init(ctx);
#endif

  vpx_rtcd();
}
