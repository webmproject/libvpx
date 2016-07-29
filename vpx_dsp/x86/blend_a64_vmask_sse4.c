/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <smmintrin.h>  // SSE4.1

#include <assert.h>

#include "vpx/vpx_integer.h"
#include "vpx_ports/mem.h"
#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_dsp/blend.h"

#include "vpx_dsp/x86/synonyms.h"
#include "vpx_dsp/x86/blend_sse4.h"

#include "./vpx_dsp_rtcd.h"

//////////////////////////////////////////////////////////////////////////////
// Implementation - No sub-sampling
//////////////////////////////////////////////////////////////////////////////

static void blend_a64_vmask_w4_sse4_1(
    uint8_t *dst, uint32_t dst_stride,
    const uint8_t *src0, uint32_t src0_stride,
    const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w) {
  const __m128i v_maxval_w = _mm_set1_epi16(VPX_BLEND_A64_MAX_ALPHA);

  (void)w;

  do {
    const __m128i v_m0_w = _mm_set1_epi16(*mask);
    const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

    const __m128i v_res_w = blend_4(src0, src1, v_m0_w, v_m1_w);

    const __m128i v_res_b = _mm_packus_epi16(v_res_w, v_res_w);

    xx_storel_32(dst, v_res_b);

    dst += dst_stride;
    src0 += src0_stride;
    src1 += src1_stride;
    mask += 1;
  } while (--h);
}

static void blend_a64_vmask_w8_sse4_1(
    uint8_t *dst, uint32_t dst_stride,
    const uint8_t *src0, uint32_t src0_stride,
    const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w) {
  const __m128i v_maxval_w = _mm_set1_epi16(VPX_BLEND_A64_MAX_ALPHA);

  (void)w;

  do {
    const __m128i v_m0_w = _mm_set1_epi16(*mask);
    const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

    const __m128i v_res_w = blend_8(src0, src1, v_m0_w, v_m1_w);

    const __m128i v_res_b = _mm_packus_epi16(v_res_w, v_res_w);

    xx_storel_64(dst, v_res_b);

    dst += dst_stride;
    src0 += src0_stride;
    src1 += src1_stride;
    mask += 1;
  } while (--h);
}

static void blend_a64_vmask_w16n_sse4_1(
    uint8_t *dst, uint32_t dst_stride,
    const uint8_t *src0, uint32_t src0_stride,
    const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w) {
  const __m128i v_maxval_w = _mm_set1_epi16(VPX_BLEND_A64_MAX_ALPHA);

  do {
    int c;
    const __m128i v_m0_w = _mm_set1_epi16(*mask);
    const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);
    for (c = 0; c < w; c += 16) {
      const __m128i v_resl_w = blend_8(src0 + c, src1 + c,
                                       v_m0_w, v_m1_w);
      const __m128i v_resh_w = blend_8(src0 + c + 8, src1 + c + 8,
                                       v_m0_w, v_m1_w);

      const __m128i v_res_b = _mm_packus_epi16(v_resl_w, v_resh_w);

      xx_storeu_128(dst + c, v_res_b);
    }
    dst += dst_stride;
    src0 += src0_stride;
    src1 += src1_stride;
    mask += 1;
  } while (--h);
}

//////////////////////////////////////////////////////////////////////////////
// Dispatch
//////////////////////////////////////////////////////////////////////////////

void vpx_blend_a64_vmask_sse4_1(
    uint8_t *dst, uint32_t dst_stride,
    const uint8_t *src0, uint32_t src0_stride,
    const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w) {
  typedef  void (*blend_fn)(uint8_t *dst, uint32_t dst_stride,
                            const uint8_t *src0, uint32_t src0_stride,
                            const uint8_t *src1, uint32_t src1_stride,
                            const uint8_t *mask, int h, int w);

  // Dimension: width_index
  static const blend_fn blend[9] = {
    blend_a64_vmask_w16n_sse4_1,  // w % 16 == 0
    vpx_blend_a64_vmask_c,        // w == 1
    vpx_blend_a64_vmask_c,        // w == 2
    NULL,                         // INVALID
    blend_a64_vmask_w4_sse4_1,    // w == 4
    NULL,                         // INVALID
    NULL,                         // INVALID
    NULL,                         // INVALID
    blend_a64_vmask_w8_sse4_1,    // w == 8
  };

  assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
  assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

  assert(h >= 1);
  assert(w >= 1);
  assert(IS_POWER_OF_TWO(h));
  assert(IS_POWER_OF_TWO(w));

  blend[w & 0xf](dst, dst_stride,
                 src0, src0_stride,
                 src1, src1_stride,
                 mask, h, w);
}

#if CONFIG_VPX_HIGHBITDEPTH
//////////////////////////////////////////////////////////////////////////////
// Implementation - No sub-sampling
//////////////////////////////////////////////////////////////////////////////

static INLINE void blend_a64_vmask_bn_w4_sse4_1(
    uint16_t *dst, uint32_t dst_stride,
    const uint16_t *src0, uint32_t src0_stride,
    const uint16_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, blend_unit_fn blend) {
  const __m128i v_maxval_w = _mm_set1_epi16(VPX_BLEND_A64_MAX_ALPHA);

  do {
    const __m128i v_m0_w = _mm_set1_epi16(*mask);
    const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

    const __m128i v_res_w = blend(src0, src1, v_m0_w, v_m1_w);

    xx_storel_64(dst, v_res_w);

    dst += dst_stride;
    src0 += src0_stride;
    src1 += src1_stride;
    mask += 1;
  } while (--h);
}

static void blend_a64_vmask_b10_w4_sse4_1(
    uint16_t *dst, uint32_t dst_stride,
    const uint16_t *src0, uint32_t src0_stride,
    const uint16_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w) {
  (void)w;
  blend_a64_vmask_bn_w4_sse4_1(dst, dst_stride, src0, src0_stride, src1,
                               src1_stride, mask, h,
                               blend_4_b10);
}

static void blend_a64_vmask_b12_w4_sse4_1(
    uint16_t *dst, uint32_t dst_stride,
    const uint16_t *src0, uint32_t src0_stride,
    const uint16_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w) {
  (void)w;
  blend_a64_vmask_bn_w4_sse4_1(dst, dst_stride, src0, src0_stride, src1,
                               src1_stride, mask, h,
                               blend_4_b12);
}

static INLINE void blend_a64_vmask_bn_w8n_sse4_1(
    uint16_t *dst, uint32_t dst_stride,
    const uint16_t *src0, uint32_t src0_stride,
    const uint16_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w, blend_unit_fn blend) {
  const __m128i v_maxval_w = _mm_set1_epi16(VPX_BLEND_A64_MAX_ALPHA);

  do {
    int c;
    const __m128i v_m0_w = _mm_set1_epi16(*mask);
    const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);
    for (c = 0; c < w; c += 8) {
      const __m128i v_res_w = blend(src0 + c, src1 + c, v_m0_w, v_m1_w);

      xx_storeu_128(dst + c, v_res_w);
    }
    dst += dst_stride;
    src0 += src0_stride;
    src1 += src1_stride;
    mask += 1;
  } while (--h);
}

static void blend_a64_vmask_b10_w8n_sse4_1(
    uint16_t *dst, uint32_t dst_stride,
    const uint16_t *src0, uint32_t src0_stride,
    const uint16_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w) {
  blend_a64_vmask_bn_w8n_sse4_1(dst, dst_stride, src0, src0_stride, src1,
                                src1_stride, mask, h, w,
                                blend_8_b10);
}

static void blend_a64_vmask_b12_w8n_sse4_1(
    uint16_t *dst, uint32_t dst_stride,
    const uint16_t *src0, uint32_t src0_stride,
    const uint16_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w) {
  blend_a64_vmask_bn_w8n_sse4_1(dst, dst_stride, src0, src0_stride, src1,
                                src1_stride, mask, h, w,
                                blend_8_b12);
}

//////////////////////////////////////////////////////////////////////////////
// Dispatch
//////////////////////////////////////////////////////////////////////////////

void vpx_highbd_blend_a64_vmask_sse4_1(
    uint8_t *dst_8, uint32_t dst_stride,
    const uint8_t *src0_8, uint32_t src0_stride,
    const uint8_t *src1_8, uint32_t src1_stride,
    const uint8_t *mask, int h, int w, int bd) {
  typedef  void (*blend_fn)(uint16_t *dst, uint32_t dst_stride,
                            const uint16_t *src0, uint32_t src0_stride,
                            const uint16_t *src1, uint32_t src1_stride,
                            const uint8_t *mask, int h, int w);

  // Dimensions are: bd_index X width_index
  static const blend_fn blend[2][2] = {
    {     // bd == 8 or 10
      blend_a64_vmask_b10_w8n_sse4_1,  // w % 8 == 0
      blend_a64_vmask_b10_w4_sse4_1,   // w == 4
    }, {  // bd == 12
      blend_a64_vmask_b12_w8n_sse4_1,  // w % 8 == 0
      blend_a64_vmask_b12_w4_sse4_1,   // w == 4
    }
  };

  assert(IMPLIES(src0_8 == dst_8, src0_stride == dst_stride));
  assert(IMPLIES(src1_8 == dst_8, src1_stride == dst_stride));

  assert(h >= 1);
  assert(w >= 1);
  assert(IS_POWER_OF_TWO(h));
  assert(IS_POWER_OF_TWO(w));

  assert(bd == 8 || bd == 10 || bd == 12);

  if (UNLIKELY((h | w) & 3)) {  // if (w <= 2 || h <= 2)
    vpx_highbd_blend_a64_vmask_c(dst_8, dst_stride,
                                 src0_8, src0_stride,
                                 src1_8, src1_stride,
                                 mask, h, w, bd);
  } else {
    uint16_t *const dst = CONVERT_TO_SHORTPTR(dst_8);
    const uint16_t *const src0 = CONVERT_TO_SHORTPTR(src0_8);
    const uint16_t *const src1 = CONVERT_TO_SHORTPTR(src1_8);

    blend[bd == 12][(w >> 2) & 1](dst, dst_stride,
                                  src0, src0_stride,
                                  src1, src1_stride,
                                  mask, h, w);
  }
}
#endif  // CONFIG_VPX_HIGHBITDEPTH
