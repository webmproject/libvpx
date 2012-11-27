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
#include "vpx_ports/x86.h"
#include "vp9/common/vp9_loopfilter.h"
#include "vp9/common/vp9_pragmas.h"
#include "vp9/common/vp9_onyxc_int.h"

void vp9_arch_x86_common_init(VP9_COMMON *ctx) {
#if CONFIG_RUNTIME_CPU_DETECT
  VP9_COMMON_RTCD *rtcd = &ctx->rtcd;
  int flags = x86_simd_caps();

  /* Note:
   *
   * This platform can be built without runtime CPU detection as well. If
   * you modify any of the function mappings present in this file, be sure
   * to also update them in static mapings (<arch>/filename_<arch>.h)
   */

  /* Override default functions with fastest ones for this CPU. */
#if HAVE_MMX
// The commented functions need to be re-written for vpx.
  if (flags & HAS_MMX) {

#if CONFIG_POSTPROC
    rtcd->postproc.down        = vp9_mbpost_proc_down_mmx;
    /*rtcd->postproc.across      = vp9_mbpost_proc_across_ip_c;*/
    rtcd->postproc.downacross  = vp9_post_proc_down_and_across_mmx;
    rtcd->postproc.addnoise    = vp9_plane_add_noise_mmx;
#endif
  }

#endif
#if HAVE_SSE2

  if (flags & HAS_SSE2) {


    // rtcd->idct.iwalsh16     = vp9_short_inv_walsh4x4_sse2;

#if CONFIG_POSTPROC
    rtcd->postproc.down        = vp9_mbpost_proc_down_xmm;
    rtcd->postproc.across      = vp9_mbpost_proc_across_ip_xmm;
    rtcd->postproc.downacross  = vp9_post_proc_down_and_across_xmm;
    rtcd->postproc.addnoise    = vp9_plane_add_noise_wmt;
#endif
  }

#endif

#if HAVE_SSSE3

  if (flags & HAS_SSSE3) {

    /* these are disable because of unsupported diagonal pred modes
    rtcd->recon.build_intra_predictors_mbuv =
      vp9_build_intra_predictors_mbuv_ssse3;
    rtcd->recon.build_intra_predictors_mbuv_s =
      vp9_build_intra_predictors_mbuv_s_ssse3;
      */
  }
#endif

#endif
}
