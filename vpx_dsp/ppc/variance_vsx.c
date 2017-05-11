/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/ppc/types_vsx.h"

void vpx_comp_avg_pred_vsx(uint8_t *comp_pred, const uint8_t *pred, int width,
                           int height, const uint8_t *ref, int ref_stride) {
  int i, j;
  /* comp_pred and pred must be 16 byte aligned. */
  assert(((intptr_t)comp_pred & 0xf) == 0);
  assert(((intptr_t)pred & 0xf) == 0);
  if (width >= 16) {
    for (i = 0; i < height; ++i) {
      for (j = 0; j < width; j += 16) {
        const uint8x16_t v = vec_avg(vec_vsx_ld(j, pred), vec_vsx_ld(j, ref));
        vec_vsx_st(v, j, comp_pred);
      }
      comp_pred += width;
      pred += width;
      ref += ref_stride;
    }
  } else if (width == 8) {
    // Process 2 lines at time
    for (i = 0; i < height / 2; ++i) {
      const uint8x16_t r0 = vec_vsx_ld(0, ref);
      const uint8x16_t r1 = vec_vsx_ld(0, ref + ref_stride);
      const uint8x16_t r = xxpermdi(r0, r1, 0);
      const uint8x16_t v = vec_avg(vec_vsx_ld(0, pred), r);
      vec_vsx_st(v, 0, comp_pred);
      comp_pred += 16;  // width * 2;
      pred += 16;       // width * 2;
      ref += ref_stride * 2;
    }
  } else {
    assert(width == 4);
    // process 4 lines at time
    for (i = 0; i < height / 4; ++i) {
      const uint32x4_t r0 = (uint32x4_t)vec_vsx_ld(0, ref);
      const uint32x4_t r1 = (uint32x4_t)vec_vsx_ld(0, ref + ref_stride);
      const uint32x4_t r2 = (uint32x4_t)vec_vsx_ld(0, ref + ref_stride * 2);
      const uint32x4_t r3 = (uint32x4_t)vec_vsx_ld(0, ref + ref_stride * 3);
      const uint8x16_t r =
          (uint8x16_t)xxpermdi(vec_mergeh(r0, r1), vec_mergeh(r2, r3), 0);
      const uint8x16_t v = vec_avg(vec_vsx_ld(0, pred), r);
      vec_vsx_st(v, 0, comp_pred);
      comp_pred += 16;  // width * 4;
      pred += 16;       // width * 4;
      ref += ref_stride * 4;
    }
  }
}
