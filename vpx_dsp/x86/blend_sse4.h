/*
*  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
*
*  Use of this source code is governed by a BSD-style license
*  that can be found in the LICENSE file in the root of the source
*  tree. An additional intellectual property rights grant can be found
*  in the file PATENTS.  All contributing project authors may
*  be found in the AUTHORS file in the root of the source tree.
*/

#ifndef VPX_DSP_X86_BLEND_SSE4_H_
#define VPX_DSP_X86_BLEND_SSE4_H_

#include "vpx_dsp/blend.h"
#include "vpx_dsp/x86/synonyms.h"

//////////////////////////////////////////////////////////////////////////////
// Common kernels
//////////////////////////////////////////////////////////////////////////////

static INLINE __m128i blend_4(const uint8_t *src0, const uint8_t *src1,
                              const __m128i v_m0_w, const __m128i v_m1_w) {
  const __m128i v_s0_b = xx_loadl_32(src0);
  const __m128i v_s1_b = xx_loadl_32(src1);
  const __m128i v_s0_w = _mm_cvtepu8_epi16(v_s0_b);
  const __m128i v_s1_w = _mm_cvtepu8_epi16(v_s1_b);

  const __m128i v_p0_w = _mm_mullo_epi16(v_s0_w, v_m0_w);
  const __m128i v_p1_w = _mm_mullo_epi16(v_s1_w, v_m1_w);

  const __m128i v_sum_w = _mm_add_epi16(v_p0_w, v_p1_w);

  const __m128i v_res_w = xx_roundn_epu16(v_sum_w, VPX_BLEND_A64_ROUND_BITS);

  return v_res_w;
}

static INLINE __m128i blend_8(const uint8_t *src0, const uint8_t *src1,
                              const __m128i v_m0_w, const __m128i v_m1_w) {
  const __m128i v_s0_b = xx_loadl_64(src0);
  const __m128i v_s1_b = xx_loadl_64(src1);
  const __m128i v_s0_w = _mm_cvtepu8_epi16(v_s0_b);
  const __m128i v_s1_w = _mm_cvtepu8_epi16(v_s1_b);

  const __m128i v_p0_w = _mm_mullo_epi16(v_s0_w, v_m0_w);
  const __m128i v_p1_w = _mm_mullo_epi16(v_s1_w, v_m1_w);

  const __m128i v_sum_w = _mm_add_epi16(v_p0_w, v_p1_w);

  const __m128i v_res_w = xx_roundn_epu16(v_sum_w, VPX_BLEND_A64_ROUND_BITS);

  return v_res_w;
}

#if CONFIG_VP9_HIGHBITDEPTH
typedef __m128i (*blend_unit_fn)(const uint16_t *src0, const uint16_t *src1,
                                 const __m128i v_m0_w, const __m128i v_m1_w);

static INLINE __m128i blend_4_b10(const uint16_t *src0, const uint16_t *src1,
                                  const __m128i v_m0_w, const __m128i v_m1_w) {
  const __m128i v_s0_w = xx_loadl_64(src0);
  const __m128i v_s1_w = xx_loadl_64(src1);

  const __m128i v_p0_w = _mm_mullo_epi16(v_s0_w, v_m0_w);
  const __m128i v_p1_w = _mm_mullo_epi16(v_s1_w, v_m1_w);

  const __m128i v_sum_w = _mm_add_epi16(v_p0_w, v_p1_w);

  const __m128i v_res_w = xx_roundn_epu16(v_sum_w, VPX_BLEND_A64_ROUND_BITS);

  return v_res_w;
}

static INLINE __m128i blend_8_b10(const uint16_t *src0, const uint16_t *src1,
                                  const __m128i v_m0_w, const __m128i v_m1_w) {
  const __m128i v_s0_w = xx_loadu_128(src0);
  const __m128i v_s1_w = xx_loadu_128(src1);

  const __m128i v_p0_w = _mm_mullo_epi16(v_s0_w, v_m0_w);
  const __m128i v_p1_w = _mm_mullo_epi16(v_s1_w, v_m1_w);

  const __m128i v_sum_w = _mm_add_epi16(v_p0_w, v_p1_w);

  const __m128i v_res_w = xx_roundn_epu16(v_sum_w, VPX_BLEND_A64_ROUND_BITS);

  return v_res_w;
}

static INLINE __m128i blend_4_b12(const uint16_t *src0, const uint16_t *src1,
                                  const __m128i v_m0_w, const __m128i v_m1_w) {
  const __m128i v_s0_w = xx_loadl_64(src0);
  const __m128i v_s1_w = xx_loadl_64(src1);

  // Interleave
  const __m128i v_m01_w = _mm_unpacklo_epi16(v_m0_w, v_m1_w);
  const __m128i v_s01_w = _mm_unpacklo_epi16(v_s0_w, v_s1_w);

  // Multiply-Add
  const __m128i v_sum_d = _mm_madd_epi16(v_s01_w, v_m01_w);

  // Scale
  const __m128i v_ssum_d = _mm_srli_epi32(v_sum_d,
                                          VPX_BLEND_A64_ROUND_BITS - 1);

  // Pack
  const __m128i v_pssum_d = _mm_packs_epi32(v_ssum_d, v_ssum_d);

  // Round
  const __m128i v_res_w = xx_round_epu16(v_pssum_d);

  return v_res_w;
}

static INLINE __m128i blend_8_b12(const uint16_t *src0, const uint16_t *src1,
                                  const __m128i v_m0_w, const __m128i v_m1_w) {
  const __m128i v_s0_w = xx_loadu_128(src0);
  const __m128i v_s1_w = xx_loadu_128(src1);

  // Interleave
  const __m128i v_m01l_w = _mm_unpacklo_epi16(v_m0_w, v_m1_w);
  const __m128i v_m01h_w = _mm_unpackhi_epi16(v_m0_w, v_m1_w);
  const __m128i v_s01l_w = _mm_unpacklo_epi16(v_s0_w, v_s1_w);
  const __m128i v_s01h_w = _mm_unpackhi_epi16(v_s0_w, v_s1_w);

  // Multiply-Add
  const __m128i v_suml_d = _mm_madd_epi16(v_s01l_w, v_m01l_w);
  const __m128i v_sumh_d = _mm_madd_epi16(v_s01h_w, v_m01h_w);

  // Scale
  const __m128i v_ssuml_d = _mm_srli_epi32(v_suml_d,
                                           VPX_BLEND_A64_ROUND_BITS - 1);
  const __m128i v_ssumh_d = _mm_srli_epi32(v_sumh_d,
                                           VPX_BLEND_A64_ROUND_BITS - 1);

  // Pack
  const __m128i v_pssum_d = _mm_packs_epi32(v_ssuml_d, v_ssumh_d);

  // Round
  const __m128i v_res_w = xx_round_epu16(v_pssum_d);

  return v_res_w;
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

#endif  // VPX_DSP_X86_BLEND_SSE4_H_
