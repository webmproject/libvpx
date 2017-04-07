/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/ppc/types_vsx.h"

void vpx_v_predictor_16x16_vsx(uint8_t *dst, ptrdiff_t stride,
                               const uint8_t *above, const uint8_t *left) {
  const uint8x16_t d = vec_vsx_ld(0, above);
  int i;
  (void)left;

  for (i = 0; i < 16; i++, dst += stride) {
    vec_vsx_st(d, 0, dst);
  }
}

void vpx_v_predictor_32x32_vsx(uint8_t *dst, ptrdiff_t stride,
                               const uint8_t *above, const uint8_t *left) {
  const uint8x16_t d0 = vec_vsx_ld(0, above);
  const uint8x16_t d1 = vec_vsx_ld(16, above);
  int i;
  (void)left;

  for (i = 0; i < 32; i++, dst += stride) {
    vec_vsx_st(d0, 0, dst);
    vec_vsx_st(d1, 16, dst);
  }
}
