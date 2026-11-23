/*
 *  Copyright (c) 2026 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <immintrin.h>

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/x86/variance_common_x86.h"
#include "vpx_ports/mem.h"

// Dynamic range note for all highbd variance / subpixel variance kernels:
// - Pixel inputs are unsigned 8/10/12-bit values in [0, 4095].
// - Bilinear interpolation computes:
//     interp = a + (((b - a) * offset + 4) >> 3)
//   for offset in {1..7}, which is exact to ((a * (16 - 2*off) + b * 2*off + 8)
//   >> 4). Since |b - a| <= 4095 and offset <= 7, |(b - a) * offset + 4| <=
//   28669 < 32767, so all intermediates fit in signed int16_t without overflow.
// - Per-pixel difference diff = src - ref is in [-4095, 4095].
// - Up to 8 differences are summed in int16_t (max 8 * 4095 = 32760 <= 32767)
//   before widening to int32_t via _mm256_madd_epi16(sum16, ones).
// - Squared difference sums via _mm256_madd_epi16(diff, diff) produce at most
//   2 * 4095^2 = 33,538,050 per int32_t lane per 16-pixel vector. Flushing to
//   uint64_t every <= 16 rows (max 64 pairs/lane = 2,146,435,200 <= INT32_MAX)
//   prevents 32-bit overflow even for 12-bit 64x64 blocks.

static INLINE void highbd_var_accum_avx2(const __m256i s, const __m256i r,
                                         __m256i *const sum16,
                                         __m256i *const sse32) {
  const __m256i diff = _mm256_sub_epi16(s, r);
  *sum16 = _mm256_add_epi16(*sum16, diff);
  *sse32 = _mm256_add_epi32(*sse32, _mm256_madd_epi16(diff, diff));
}

static INLINE void flush_sum16_to_sum32_avx2(__m256i *const sum16,
                                             __m256i *const sum32) {
  const __m256i ones = _mm256_set1_epi16(1);
  *sum32 = _mm256_add_epi32(*sum32, _mm256_madd_epi16(*sum16, ones));
  *sum16 = _mm256_setzero_si256();
}

static INLINE void flush_sse32_to_sse64_avx2(__m256i *const sse32,
                                             __m256i *const sse64) {
  const __m256i lo = _mm256_cvtepu32_epi64(_mm256_castsi256_si128(*sse32));
  const __m256i hi = _mm256_cvtepu32_epi64(_mm256_extracti128_si256(*sse32, 1));
  *sse64 = _mm256_add_epi64(*sse64, _mm256_add_epi64(lo, hi));
  *sse32 = _mm256_setzero_si256();
}

static INLINE void highbd_var_final_avx2(const __m256i sum32,
                                         const __m256i sse64,
                                         uint64_t *const sse,
                                         int64_t *const sum) {
  const __m128i sum128 = _mm_add_epi32(_mm256_castsi256_si128(sum32),
                                       _mm256_extracti128_si256(sum32, 1));
  const __m128i sum64 = _mm_add_epi32(sum128, _mm_srli_si128(sum128, 8));
  const __m128i sum_final = _mm_add_epi32(sum64, _mm_srli_si128(sum64, 4));
  *sum = _mm_cvtsi128_si32(sum_final);

  const __m128i sse128 = _mm_add_epi64(_mm256_castsi256_si128(sse64),
                                       _mm256_extracti128_si256(sse64, 1));
  const __m128i sse_final = _mm_add_epi64(sse128, _mm_srli_si128(sse128, 8));
  _mm_storel_epi64((__m128i *)sse, sse_final);
}

static INLINE void highbd_variance_w4_avx2(const uint16_t *src, int src_stride,
                                           const uint16_t *ref, int ref_stride,
                                           int h, uint64_t *const sse,
                                           int64_t *const sum) {
  __m256i sum16 = _mm256_setzero_si256();
  __m256i sum32 = _mm256_setzero_si256();
  __m256i sse32 = _mm256_setzero_si256();
  __m256i sse64 = _mm256_setzero_si256();
  int i;
  assert(h % 4 == 0);
  for (i = 0; i < h; i += 4) {
    const __m128i s0 = _mm_loadl_epi64((const __m128i *)(src + 0 * src_stride));
    const __m128i s1 = _mm_loadl_epi64((const __m128i *)(src + 1 * src_stride));
    const __m128i s2 = _mm_loadl_epi64((const __m128i *)(src + 2 * src_stride));
    const __m128i s3 = _mm_loadl_epi64((const __m128i *)(src + 3 * src_stride));
    const __m128i r0 = _mm_loadl_epi64((const __m128i *)(ref + 0 * ref_stride));
    const __m128i r1 = _mm_loadl_epi64((const __m128i *)(ref + 1 * ref_stride));
    const __m128i r2 = _mm_loadl_epi64((const __m128i *)(ref + 2 * ref_stride));
    const __m128i r3 = _mm_loadl_epi64((const __m128i *)(ref + 3 * ref_stride));
    const __m256i s = _mm256_inserti128_si256(
        _mm256_castsi128_si256(_mm_unpacklo_epi64(s0, s1)),
        _mm_unpacklo_epi64(s2, s3), 1);
    const __m256i r = _mm256_inserti128_si256(
        _mm256_castsi128_si256(_mm_unpacklo_epi64(r0, r1)),
        _mm_unpacklo_epi64(r2, r3), 1);
    highbd_var_accum_avx2(s, r, &sum16, &sse32);
    src += 4 * src_stride;
    ref += 4 * ref_stride;
  }
  flush_sum16_to_sum32_avx2(&sum16, &sum32);
  flush_sse32_to_sse64_avx2(&sse32, &sse64);
  highbd_var_final_avx2(sum32, sse64, sse, sum);
}

static INLINE void highbd_variance_w8_avx2(const uint16_t *src, int src_stride,
                                           const uint16_t *ref, int ref_stride,
                                           int h, uint64_t *const sse,
                                           int64_t *const sum) {
  __m256i sum16 = _mm256_setzero_si256();
  __m256i sum32 = _mm256_setzero_si256();
  __m256i sse32 = _mm256_setzero_si256();
  __m256i sse64 = _mm256_setzero_si256();
  int i;
  assert(h % 2 == 0 && h <= 16);
  for (i = 0; i < h; i += 2) {
    const __m128i s0 = _mm_loadu_si128((const __m128i *)src);
    const __m128i s1 = _mm_loadu_si128((const __m128i *)(src + src_stride));
    const __m128i r0 = _mm_loadu_si128((const __m128i *)ref);
    const __m128i r1 = _mm_loadu_si128((const __m128i *)(ref + ref_stride));
    const __m256i s =
        _mm256_inserti128_si256(_mm256_castsi128_si256(s0), s1, 1);
    const __m256i r =
        _mm256_inserti128_si256(_mm256_castsi128_si256(r0), r1, 1);
    highbd_var_accum_avx2(s, r, &sum16, &sse32);
    src += 2 * src_stride;
    ref += 2 * ref_stride;
  }
  flush_sum16_to_sum32_avx2(&sum16, &sum32);
  flush_sse32_to_sse64_avx2(&sse32, &sse64);
  highbd_var_final_avx2(sum32, sse64, sse, sum);
}

static INLINE void highbd_variance_w16p_avx2(
    const uint16_t *src, int src_stride, const uint16_t *ref, int ref_stride,
    int w, int h, uint64_t *const sse, int64_t *const sum) {
  __m256i sum16 = _mm256_setzero_si256();
  __m256i sum32 = _mm256_setzero_si256();
  __m256i sse32 = _mm256_setzero_si256();
  __m256i sse64 = _mm256_setzero_si256();
  int i, j;
  assert(w == 16 || w == 32 || w == 64);
  assert(h >= 8 && h <= 64 && h % 8 == 0);
  for (i = 0; i < h; ++i) {
    for (j = 0; j < w; j += 16) {
      const __m256i s = _mm256_loadu_si256((const __m256i *)(src + j));
      const __m256i r = _mm256_loadu_si256((const __m256i *)(ref + j));
      highbd_var_accum_avx2(s, r, &sum16, &sse32);
    }
    // Flush sum16 every row when w == 64 (4 additions/lane), else every 4 rows.
    if (w == 64 || (i & 3) == 3) {
      flush_sum16_to_sum32_avx2(&sum16, &sum32);
    }
    if ((i & 15) == 15) {
      flush_sse32_to_sse64_avx2(&sse32, &sse64);
    }
    src += src_stride;
    ref += ref_stride;
  }
  if (h < 16) flush_sse32_to_sse64_avx2(&sse32, &sse64);
  highbd_var_final_avx2(sum32, sse64, sse, sum);
}

#define HIGHBD_GET_VAR_AVX2(S, fn_w)                                       \
  void vpx_highbd_8_get##S##x##S##var_avx2(                                \
      const uint8_t *src8, int src_stride, const uint8_t *ref8,            \
      int ref_stride, uint32_t *sse, int *sum) {                           \
    uint64_t sse_long;                                                     \
    int64_t sum_long;                                                      \
    fn_w(CONVERT_TO_SHORTPTR(src8), src_stride, CONVERT_TO_SHORTPTR(ref8), \
         ref_stride, S, &sse_long, &sum_long);                             \
    *sse = (uint32_t)sse_long;                                             \
    *sum = (int)sum_long;                                                  \
  }                                                                        \
                                                                           \
  void vpx_highbd_10_get##S##x##S##var_avx2(                               \
      const uint8_t *src8, int src_stride, const uint8_t *ref8,            \
      int ref_stride, uint32_t *sse, int *sum) {                           \
    uint64_t sse_long;                                                     \
    int64_t sum_long;                                                      \
    fn_w(CONVERT_TO_SHORTPTR(src8), src_stride, CONVERT_TO_SHORTPTR(ref8), \
         ref_stride, S, &sse_long, &sum_long);                             \
    *sse = (uint32_t)ROUND_POWER_OF_TWO(sse_long, 4);                      \
    *sum = (int)ROUND_POWER_OF_TWO(sum_long, 2);                           \
  }                                                                        \
                                                                           \
  void vpx_highbd_12_get##S##x##S##var_avx2(                               \
      const uint8_t *src8, int src_stride, const uint8_t *ref8,            \
      int ref_stride, uint32_t *sse, int *sum) {                           \
    uint64_t sse_long;                                                     \
    int64_t sum_long;                                                      \
    fn_w(CONVERT_TO_SHORTPTR(src8), src_stride, CONVERT_TO_SHORTPTR(ref8), \
         ref_stride, S, &sse_long, &sum_long);                             \
    *sse = (uint32_t)ROUND_POWER_OF_TWO(sse_long, 8);                      \
    *sum = (int)ROUND_POWER_OF_TWO(sum_long, 4);                           \
  }

static INLINE void highbd_variance_w16_avx2(const uint16_t *src, int src_stride,
                                            const uint16_t *ref, int ref_stride,
                                            int h, uint64_t *const sse,
                                            int64_t *const sum) {
  highbd_variance_w16p_avx2(src, src_stride, ref, ref_stride, 16, h, sse, sum);
}

HIGHBD_GET_VAR_AVX2(8, highbd_variance_w8_avx2)
HIGHBD_GET_VAR_AVX2(16, highbd_variance_w16_avx2)

#undef HIGHBD_GET_VAR_AVX2

#define HIGHBD_MSE_AVX2(w, h, fn_call)                          \
  uint32_t vpx_highbd_8_mse##w##x##h##_avx2(                    \
      const uint8_t *src8, int src_stride, const uint8_t *ref8, \
      int ref_stride, uint32_t *sse) {                          \
    uint64_t sse_long;                                          \
    int64_t sum_long;                                           \
    fn_call;                                                    \
    *sse = (uint32_t)sse_long;                                  \
    return *sse;                                                \
  }                                                             \
                                                                \
  uint32_t vpx_highbd_10_mse##w##x##h##_avx2(                   \
      const uint8_t *src8, int src_stride, const uint8_t *ref8, \
      int ref_stride, uint32_t *sse) {                          \
    uint64_t sse_long;                                          \
    int64_t sum_long;                                           \
    fn_call;                                                    \
    *sse = (uint32_t)ROUND_POWER_OF_TWO(sse_long, 4);           \
    return *sse;                                                \
  }                                                             \
                                                                \
  uint32_t vpx_highbd_12_mse##w##x##h##_avx2(                   \
      const uint8_t *src8, int src_stride, const uint8_t *ref8, \
      int ref_stride, uint32_t *sse) {                          \
    uint64_t sse_long;                                          \
    int64_t sum_long;                                           \
    fn_call;                                                    \
    *sse = (uint32_t)ROUND_POWER_OF_TWO(sse_long, 8);           \
    return *sse;                                                \
  }

HIGHBD_MSE_AVX2(8, 8,
                highbd_variance_w8_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                        CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                        8, &sse_long, &sum_long))
HIGHBD_MSE_AVX2(8, 16,
                highbd_variance_w8_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                        CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                        16, &sse_long, &sum_long))
HIGHBD_MSE_AVX2(16, 8,
                highbd_variance_w16p_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                          CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                          16, 8, &sse_long, &sum_long))
HIGHBD_MSE_AVX2(16, 16,
                highbd_variance_w16p_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                          CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                          16, 16, &sse_long, &sum_long))

#undef HIGHBD_MSE_AVX2

#define HIGHBD_VAR_AVX2(w, h, shift, fn_call)                      \
  uint32_t vpx_highbd_8_variance##w##x##h##_avx2(                  \
      const uint8_t *src8, int src_stride, const uint8_t *ref8,    \
      int ref_stride, uint32_t *sse) {                             \
    uint64_t sse_long;                                             \
    int64_t sum_long;                                              \
    fn_call;                                                       \
    return finalize_highbd_var_8(sse_long, sum_long, shift, sse);  \
  }                                                                \
                                                                   \
  uint32_t vpx_highbd_10_variance##w##x##h##_avx2(                 \
      const uint8_t *src8, int src_stride, const uint8_t *ref8,    \
      int ref_stride, uint32_t *sse) {                             \
    uint64_t sse_long;                                             \
    int64_t sum_long;                                              \
    fn_call;                                                       \
    return finalize_highbd_var_10(sse_long, sum_long, shift, sse); \
  }                                                                \
                                                                   \
  uint32_t vpx_highbd_12_variance##w##x##h##_avx2(                 \
      const uint8_t *src8, int src_stride, const uint8_t *ref8,    \
      int ref_stride, uint32_t *sse) {                             \
    uint64_t sse_long;                                             \
    int64_t sum_long;                                              \
    fn_call;                                                       \
    return finalize_highbd_var_12(sse_long, sum_long, shift, sse); \
  }

HIGHBD_VAR_AVX2(4, 4, 4,
                highbd_variance_w4_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                        CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                        4, &sse_long, &sum_long))
HIGHBD_VAR_AVX2(4, 8, 5,
                highbd_variance_w4_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                        CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                        8, &sse_long, &sum_long))
HIGHBD_VAR_AVX2(8, 4, 5,
                highbd_variance_w8_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                        CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                        4, &sse_long, &sum_long))
HIGHBD_VAR_AVX2(8, 8, 6,
                highbd_variance_w8_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                        CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                        8, &sse_long, &sum_long))
HIGHBD_VAR_AVX2(8, 16, 7,
                highbd_variance_w8_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                        CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                        16, &sse_long, &sum_long))
HIGHBD_VAR_AVX2(16, 8, 7,
                highbd_variance_w16p_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                          CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                          16, 8, &sse_long, &sum_long))
HIGHBD_VAR_AVX2(16, 16, 8,
                highbd_variance_w16p_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                          CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                          16, 16, &sse_long, &sum_long))
HIGHBD_VAR_AVX2(16, 32, 9,
                highbd_variance_w16p_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                          CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                          16, 32, &sse_long, &sum_long))
HIGHBD_VAR_AVX2(32, 16, 9,
                highbd_variance_w16p_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                          CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                          32, 16, &sse_long, &sum_long))
HIGHBD_VAR_AVX2(32, 32, 10,
                highbd_variance_w16p_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                          CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                          32, 32, &sse_long, &sum_long))
HIGHBD_VAR_AVX2(32, 64, 11,
                highbd_variance_w16p_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                          CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                          32, 64, &sse_long, &sum_long))
HIGHBD_VAR_AVX2(64, 32, 11,
                highbd_variance_w16p_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                          CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                          64, 32, &sse_long, &sum_long))
HIGHBD_VAR_AVX2(64, 64, 12,
                highbd_variance_w16p_avx2(CONVERT_TO_SHORTPTR(src8), src_stride,
                                          CONVERT_TO_SHORTPTR(ref8), ref_stride,
                                          64, 64, &sse_long, &sum_long))

#undef HIGHBD_VAR_AVX2

// Bilinear interpolation helper for uint16_t pixels:
// Computes ROUND_POWER_OF_TWO(a * (8 - off) + b * off, 3) =
// a + (((b - a) * off + 4) >> 3).
static INLINE __m256i highbd_bilin_interp_avx2(const __m256i a, const __m256i b,
                                               const __m256i *const v_off,
                                               const __m256i *const pw4) {
  const __m256i diff = _mm256_sub_epi16(b, a);
  const __m256i prod = _mm256_mullo_epi16(diff, *v_off);
  const __m256i rounded = _mm256_add_epi16(prod, *pw4);
  return _mm256_add_epi16(a, _mm256_srai_epi16(rounded, 3));
}

// 4-wide subpixel variance (4x4, 4x8): processes 4 rows per iteration in YMM.
static void highbd_spv_w4_avx2(const uint16_t *src, int src_stride,
                               int x_offset, int y_offset, const uint16_t *ref,
                               int ref_stride, const uint16_t *sec, int do_sec,
                               int h, uint64_t *const sse, int64_t *const sum) {
  __m256i sum16 = _mm256_setzero_si256();
  __m256i sum32 = _mm256_setzero_si256();
  __m256i sse32 = _mm256_setzero_si256();
  __m256i sse64 = _mm256_setzero_si256();
  const __m256i pw4 = _mm256_set1_epi16(4);
  const __m256i vx_off = _mm256_set1_epi16((int16_t)x_offset);
  const __m256i vy_off = _mm256_set1_epi16((int16_t)y_offset);
  int i;
  assert(h % 4 == 0 && h <= 8);
  assert((unsigned int)x_offset < 8 && (unsigned int)y_offset < 8);

  for (i = 0; i < h; i += 4) {
    const __m128i r0 = _mm_loadl_epi64((const __m128i *)(ref + 0 * ref_stride));
    const __m128i r1 = _mm_loadl_epi64((const __m128i *)(ref + 1 * ref_stride));
    const __m128i r2 = _mm_loadl_epi64((const __m128i *)(ref + 2 * ref_stride));
    const __m128i r3 = _mm_loadl_epi64((const __m128i *)(ref + 3 * ref_stride));
    const __m256i r = _mm256_inserti128_si256(
        _mm256_castsi128_si256(_mm_unpacklo_epi64(r0, r1)),
        _mm_unpacklo_epi64(r2, r3), 1);

    const __m128i s0a =
        _mm_loadl_epi64((const __m128i *)(src + 0 * src_stride));
    const __m128i s1a =
        _mm_loadl_epi64((const __m128i *)(src + 1 * src_stride));
    const __m128i s2a =
        _mm_loadl_epi64((const __m128i *)(src + 2 * src_stride));
    const __m128i s3a =
        _mm_loadl_epi64((const __m128i *)(src + 3 * src_stride));
    __m256i top = _mm256_inserti128_si256(
        _mm256_castsi128_si256(_mm_unpacklo_epi64(s0a, s1a)),
        _mm_unpacklo_epi64(s2a, s3a), 1);

    if (x_offset != 0) {
      const __m128i s0b =
          _mm_loadl_epi64((const __m128i *)(src + 0 * src_stride + 1));
      const __m128i s1b =
          _mm_loadl_epi64((const __m128i *)(src + 1 * src_stride + 1));
      const __m128i s2b =
          _mm_loadl_epi64((const __m128i *)(src + 2 * src_stride + 1));
      const __m128i s3b =
          _mm_loadl_epi64((const __m128i *)(src + 3 * src_stride + 1));
      const __m256i top_r = _mm256_inserti128_si256(
          _mm256_castsi128_si256(_mm_unpacklo_epi64(s0b, s1b)),
          _mm_unpacklo_epi64(s2b, s3b), 1);
      if (x_offset == 4) {
        top = _mm256_avg_epu16(top, top_r);
      } else {
        top = highbd_bilin_interp_avx2(top, top_r, &vx_off, &pw4);
      }
    }

    if (y_offset != 0) {
      // Shift top down by 1 row (64 bits) and insert row 4 into highest 64
      // bits.
      const __m128i s4a =
          _mm_loadl_epi64((const __m128i *)(src + 4 * src_stride));
      __m128i bot4 = s4a;
      if (x_offset != 0) {
        const __m128i s4b =
            _mm_loadl_epi64((const __m128i *)(src + 4 * src_stride + 1));
        if (x_offset == 4) {
          bot4 = _mm_avg_epu16(s4a, s4b);
        } else {
          const __m128i diff = _mm_sub_epi16(s4b, s4a);
          const __m128i prod =
              _mm_mullo_epi16(diff, _mm256_castsi256_si128(vx_off));
          const __m128i rnd = _mm_add_epi16(prod, _mm256_castsi256_si128(pw4));
          bot4 = _mm_add_epi16(s4a, _mm_srai_epi16(rnd, 3));
        }
      }
      const __m256i shifted =
          _mm256_permute4x64_epi64(top, 0x39);  // [q1,q2,q3,q0]
      const __m256i b4_256 =
          _mm256_permute4x64_epi64(_mm256_castsi128_si256(bot4), 0x00);
      const __m256i bot = _mm256_blend_epi32(shifted, b4_256, 0xC0);
      if (y_offset == 4) {
        top = _mm256_avg_epu16(top, bot);
      } else {
        top = highbd_bilin_interp_avx2(top, bot, &vy_off, &pw4);
      }
    }

    if (do_sec) {
      const __m256i s2 = _mm256_loadu_si256((const __m256i *)sec);
      top = _mm256_avg_epu16(top, s2);
      sec += 16;
    }

    highbd_var_accum_avx2(top, r, &sum16, &sse32);
    src += 4 * src_stride;
    ref += 4 * ref_stride;
  }
  flush_sum16_to_sum32_avx2(&sum16, &sum32);
  flush_sse32_to_sse64_avx2(&sse32, &sse64);
  highbd_var_final_avx2(sum32, sse64, sse, sum);
}

// 8-wide subpixel variance (8x4, 8x8, 8x16): processes 2 rows per iter in YMM.
static VPX_FORCE_INLINE void highbd_spv_w8_core_avx2(
    const uint16_t *src, int src_stride, int x_offset, int y_offset,
    const uint16_t *ref, int ref_stride, const uint16_t *sec, int sec_stride,
    const int do_sec, int w, int h, const int x_mode, const int y_mode,
    uint64_t *const sse, int64_t *const sum) {
  (void)w;
  (void)sec_stride;
  __m256i sum16 = _mm256_setzero_si256();
  __m256i sum32 = _mm256_setzero_si256();
  __m256i sse32 = _mm256_setzero_si256();
  __m256i sse64 = _mm256_setzero_si256();
  const __m256i pw4 = _mm256_set1_epi16(4);
  const __m256i vx_off = _mm256_set1_epi16((int16_t)x_offset);
  const __m256i vy_off = _mm256_set1_epi16((int16_t)y_offset);
  int i;
  assert(h % 2 == 0 && h <= 16);
  assert((unsigned int)x_offset < 8 && (unsigned int)y_offset < 8);

  if (y_mode == 0) {
    for (i = 0; i < h; i += 2) {
      const __m128i s0 = _mm_loadu_si128((const __m128i *)src);
      const __m128i s1 = _mm_loadu_si128((const __m128i *)(src + src_stride));
      __m256i pred = _mm256_inserti128_si256(_mm256_castsi128_si256(s0), s1, 1);
      if (x_mode != 0) {
        const __m128i s0r = _mm_loadu_si128((const __m128i *)(src + 1));
        const __m128i s1r =
            _mm_loadu_si128((const __m128i *)(src + src_stride + 1));
        const __m256i pred_r =
            _mm256_inserti128_si256(_mm256_castsi128_si256(s0r), s1r, 1);
        pred = (x_mode == 1)
                   ? _mm256_avg_epu16(pred, pred_r)
                   : highbd_bilin_interp_avx2(pred, pred_r, &vx_off, &pw4);
      }
      if (do_sec) {
        const __m256i s2 = _mm256_loadu_si256((const __m256i *)sec);
        pred = _mm256_avg_epu16(pred, s2);
        sec += 16;
      }
      const __m128i r0 = _mm_loadu_si128((const __m128i *)ref);
      const __m128i r1 = _mm_loadu_si128((const __m128i *)(ref + ref_stride));
      const __m256i r =
          _mm256_inserti128_si256(_mm256_castsi128_si256(r0), r1, 1);
      highbd_var_accum_avx2(pred, r, &sum16, &sse32);
      src += 2 * src_stride;
      ref += 2 * ref_stride;
    }
  } else {
    __m128i prev_r0 = _mm_loadu_si128((const __m128i *)src);
    if (x_mode != 0) {
      const __m128i prev_r0_r = _mm_loadu_si128((const __m128i *)(src + 1));
      if (x_mode == 1) {
        prev_r0 = _mm_avg_epu16(prev_r0, prev_r0_r);
      } else {
        const __m128i diff = _mm_sub_epi16(prev_r0_r, prev_r0);
        const __m128i prod =
            _mm_mullo_epi16(diff, _mm256_castsi256_si128(vx_off));
        const __m128i rnd = _mm_add_epi16(prod, _mm256_castsi256_si128(pw4));
        prev_r0 = _mm_add_epi16(prev_r0, _mm_srai_epi16(rnd, 3));
      }
    }
    __m256i prev_pack =
        _mm256_inserti128_si256(_mm256_setzero_si256(), prev_r0, 1);
    src += src_stride;

    for (i = 0; i < h; i += 2) {
      const __m128i s1 = _mm_loadu_si128((const __m128i *)src);
      const __m128i s2 = _mm_loadu_si128((const __m128i *)(src + src_stride));
      __m256i curr_pack =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s1), s2, 1);
      if (x_mode != 0) {
        const __m128i s1r = _mm_loadu_si128((const __m128i *)(src + 1));
        const __m128i s2r =
            _mm_loadu_si128((const __m128i *)(src + src_stride + 1));
        const __m256i curr_r =
            _mm256_inserti128_si256(_mm256_castsi128_si256(s1r), s2r, 1);
        curr_pack = (x_mode == 1) ? _mm256_avg_epu16(curr_pack, curr_r)
                                  : highbd_bilin_interp_avx2(curr_pack, curr_r,
                                                             &vx_off, &pw4);
      }
      const __m256i top_pack =
          _mm256_permute2x128_si256(prev_pack, curr_pack, 0x21);
      __m256i pred = (y_mode == 1) ? _mm256_avg_epu16(top_pack, curr_pack)
                                   : highbd_bilin_interp_avx2(
                                         top_pack, curr_pack, &vy_off, &pw4);
      prev_pack = curr_pack;

      if (do_sec) {
        const __m256i s2_sec = _mm256_loadu_si256((const __m256i *)sec);
        pred = _mm256_avg_epu16(pred, s2_sec);
        sec += 16;
      }
      const __m128i r0 = _mm_loadu_si128((const __m128i *)ref);
      const __m128i r1 = _mm_loadu_si128((const __m128i *)(ref + ref_stride));
      const __m256i r =
          _mm256_inserti128_si256(_mm256_castsi128_si256(r0), r1, 1);
      highbd_var_accum_avx2(pred, r, &sum16, &sse32);
      src += 2 * src_stride;
      ref += 2 * ref_stride;
    }
  }
  flush_sum16_to_sum32_avx2(&sum16, &sum32);
  flush_sse32_to_sse64_avx2(&sse32, &sse64);
  highbd_var_final_avx2(sum32, sse64, sse, sum);
}

static void highbd_spv_w8_var_avx2(const uint16_t *src, int src_stride,
                                   int x_offset, int y_offset,
                                   const uint16_t *ref, int ref_stride,
                                   const uint16_t *sec, const int do_sec, int h,
                                   uint64_t *const sse, int64_t *const sum) {
  (void)sec;
  (void)do_sec;
  VPX_DISPATCH_SPV_HBD(highbd_spv_w8_core_avx2, src, src_stride, x_offset,
                       y_offset, ref, ref_stride, NULL, 8, 0, 8, h, sse, sum);
}

static void highbd_spv_w8_avg_avx2(const uint16_t *src, int src_stride,
                                   int x_offset, int y_offset,
                                   const uint16_t *ref, int ref_stride,
                                   const uint16_t *sec, const int do_sec, int h,
                                   uint64_t *const sse, int64_t *const sum) {
  (void)do_sec;
  VPX_DISPATCH_SPV_HBD(highbd_spv_w8_core_avx2, src, src_stride, x_offset,
                       y_offset, ref, ref_stride, sec, 8, 1, 8, h, sse, sum);
}

// Horizontally filter 16 uint16_t pixels in a single YMM row (x_mode: 0=zero,
// 1=half, 2=bilin).
static VPX_FORCE_INLINE __m256i
highbd_horiz_w16_avx2(const uint16_t *src, const int x_mode,
                      const __m256i *const vx_off, const __m256i *const pw4) {
  const __m256i s0 = _mm256_loadu_si256((const __m256i *)src);
  if (x_mode == 0) return s0;
  const __m256i s1 = _mm256_loadu_si256((const __m256i *)(src + 1));
  if (x_mode == 1) return _mm256_avg_epu16(s0, s1);
  return highbd_bilin_interp_avx2(s0, s1, vx_off, pw4);
}

// Parameterized row-by-row subpixel variance kernel for w in {16, 32} (w == 64
// runs as two 32-wide strips with sec_stride == 64).
// x_mode, y_mode in {0=zero, 1=half, 2=bilin}; do_sec in {0, 1}.
static VPX_FORCE_INLINE void highbd_spv_w16p_core_avx2(
    const uint16_t *src, int src_stride, int x_offset, int y_offset,
    const uint16_t *ref, int ref_stride, const uint16_t *sec, int sec_stride,
    int do_sec, int w, int h, const int x_mode, const int y_mode,
    uint64_t *const sse, int64_t *const sum) {
  __m256i sum16 = _mm256_setzero_si256();
  __m256i sum32 = _mm256_setzero_si256();
  __m256i sse32 = _mm256_setzero_si256();
  __m256i sse64 = _mm256_setzero_si256();
  const __m256i pw4 = _mm256_set1_epi16(4);
  const __m256i vx_off = _mm256_set1_epi16((int16_t)x_offset);
  const __m256i vy_off = _mm256_set1_epi16((int16_t)y_offset);
  int i;
  assert(w == 16 || w == 32);
  assert(h >= 8 && h <= 64 && h % 8 == 0);
  assert((unsigned int)x_offset < 8 && (unsigned int)y_offset < 8);

  if (y_mode == 0) {
    for (i = 0; i < h; ++i) {
      int j;
      for (j = 0; j < w; j += 16) {
        __m256i pred = highbd_horiz_w16_avx2(src + j, x_mode, &vx_off, &pw4);
        if (do_sec) {
          const __m256i s2 = _mm256_loadu_si256((const __m256i *)(sec + j));
          pred = _mm256_avg_epu16(pred, s2);
        }
        const __m256i r = _mm256_loadu_si256((const __m256i *)(ref + j));
        highbd_var_accum_avx2(pred, r, &sum16, &sse32);
      }
      if ((i & 3) == 3) flush_sum16_to_sum32_avx2(&sum16, &sum32);
      if ((i & 15) == 15) flush_sse32_to_sse64_avx2(&sse32, &sse64);
      src += src_stride;
      ref += ref_stride;
      if (do_sec) sec += sec_stride;
    }
  } else {
    __m256i prev0 = highbd_horiz_w16_avx2(src + 0, x_mode, &vx_off, &pw4);
    __m256i prev1 = (w >= 32)
                        ? highbd_horiz_w16_avx2(src + 16, x_mode, &vx_off, &pw4)
                        : _mm256_setzero_si256();
    src += src_stride;
    for (i = 0; i < h; ++i) {
      const __m256i curr0 =
          highbd_horiz_w16_avx2(src + 0, x_mode, &vx_off, &pw4);
      __m256i pred0 =
          (y_mode == 1) ? _mm256_avg_epu16(prev0, curr0)
                        : highbd_bilin_interp_avx2(prev0, curr0, &vy_off, &pw4);
      prev0 = curr0;
      if (do_sec) {
        pred0 =
            _mm256_avg_epu16(pred0, _mm256_loadu_si256((const __m256i *)sec));
      }
      highbd_var_accum_avx2(pred0, _mm256_loadu_si256((const __m256i *)ref),
                            &sum16, &sse32);

      if (w >= 32) {
        const __m256i curr1 =
            highbd_horiz_w16_avx2(src + 16, x_mode, &vx_off, &pw4);
        __m256i pred1 = (y_mode == 1) ? _mm256_avg_epu16(prev1, curr1)
                                      : highbd_bilin_interp_avx2(prev1, curr1,
                                                                 &vy_off, &pw4);
        prev1 = curr1;
        if (do_sec) {
          pred1 = _mm256_avg_epu16(
              pred1, _mm256_loadu_si256((const __m256i *)(sec + 16)));
        }
        highbd_var_accum_avx2(pred1,
                              _mm256_loadu_si256((const __m256i *)(ref + 16)),
                              &sum16, &sse32);
      }

      if ((i & 3) == 3) flush_sum16_to_sum32_avx2(&sum16, &sum32);
      if ((i & 15) == 15) flush_sse32_to_sse64_avx2(&sse32, &sse64);
      src += src_stride;
      ref += ref_stride;
      if (do_sec) sec += sec_stride;
    }
  }
  if (h < 16) flush_sse32_to_sse64_avx2(&sse32, &sse64);
  highbd_var_final_avx2(sum32, sse64, sse, sum);
}

static void highbd_spv_w16_var_avx2(const uint16_t *src, int src_stride,
                                    int x_offset, int y_offset,
                                    const uint16_t *ref, int ref_stride,
                                    const uint16_t *sec, const int do_sec,
                                    int h, uint64_t *const sse,
                                    int64_t *const sum) {
  (void)sec;
  (void)do_sec;
  VPX_DISPATCH_SPV_HBD(highbd_spv_w16p_core_avx2, src, src_stride, x_offset,
                       y_offset, ref, ref_stride, NULL, 16, 0, 16, h, sse, sum);
}

static void highbd_spv_w16_avg_avx2(const uint16_t *src, int src_stride,
                                    int x_offset, int y_offset,
                                    const uint16_t *ref, int ref_stride,
                                    const uint16_t *sec, const int do_sec,
                                    int h, uint64_t *const sse,
                                    int64_t *const sum) {
  (void)do_sec;
  VPX_DISPATCH_SPV_HBD(highbd_spv_w16p_core_avx2, src, src_stride, x_offset,
                       y_offset, ref, ref_stride, sec, 16, 1, 16, h, sse, sum);
}

static void highbd_spv_w32_var_stride_avx2(const uint16_t *src, int src_stride,
                                           int x_offset, int y_offset,
                                           const uint16_t *ref, int ref_stride,
                                           int h, uint64_t *const sse,
                                           int64_t *const sum) {
  VPX_DISPATCH_SPV_HBD(highbd_spv_w16p_core_avx2, src, src_stride, x_offset,
                       y_offset, ref, ref_stride, NULL, 32, 0, 32, h, sse, sum);
}

static void highbd_spv_w32_avg_stride_avx2(const uint16_t *src, int src_stride,
                                           int x_offset, int y_offset,
                                           const uint16_t *ref, int ref_stride,
                                           const uint16_t *sec, int sec_stride,
                                           int h, uint64_t *const sse,
                                           int64_t *const sum) {
  VPX_DISPATCH_SPV_HBD(highbd_spv_w16p_core_avx2, src, src_stride, x_offset,
                       y_offset, ref, ref_stride, sec, sec_stride, 1, 32, h,
                       sse, sum);
}

static INLINE void highbd_spv_w32_var_avx2(
    const uint16_t *src, int src_stride, int x_offset, int y_offset,
    const uint16_t *ref, int ref_stride, const uint16_t *sec, const int do_sec,
    int h, uint64_t *const sse, int64_t *const sum) {
  (void)sec;
  (void)do_sec;
  highbd_spv_w32_var_stride_avx2(src, src_stride, x_offset, y_offset, ref,
                                 ref_stride, h, sse, sum);
}

static INLINE void highbd_spv_w32_avg_avx2(
    const uint16_t *src, int src_stride, int x_offset, int y_offset,
    const uint16_t *ref, int ref_stride, const uint16_t *sec, const int do_sec,
    int h, uint64_t *const sse, int64_t *const sum) {
  (void)do_sec;
  highbd_spv_w32_avg_stride_avx2(src, src_stride, x_offset, y_offset, ref,
                                 ref_stride, sec, 32, h, sse, sum);
}

static void highbd_spv_w64_var_avx2(const uint16_t *src, int src_stride,
                                    int x_offset, int y_offset,
                                    const uint16_t *ref, int ref_stride,
                                    const uint16_t *sec, const int do_sec,
                                    int h, uint64_t *const sse,
                                    int64_t *const sum) {
  uint64_t sse1, sse2;
  int64_t sum1, sum2;
  (void)sec;
  (void)do_sec;
  highbd_spv_w32_var_stride_avx2(src, src_stride, x_offset, y_offset, ref,
                                 ref_stride, h, &sse1, &sum1);
  highbd_spv_w32_var_stride_avx2(src + 32, src_stride, x_offset, y_offset,
                                 ref + 32, ref_stride, h, &sse2, &sum2);
  *sse = sse1 + sse2;
  *sum = sum1 + sum2;
}

static void highbd_spv_w64_avg_avx2(const uint16_t *src, int src_stride,
                                    int x_offset, int y_offset,
                                    const uint16_t *ref, int ref_stride,
                                    const uint16_t *sec, const int do_sec,
                                    int h, uint64_t *const sse,
                                    int64_t *const sum) {
  uint64_t sse1, sse2;
  int64_t sum1, sum2;
  (void)do_sec;
  highbd_spv_w32_avg_stride_avx2(src, src_stride, x_offset, y_offset, ref,
                                 ref_stride, sec, 64, h, &sse1, &sum1);
  highbd_spv_w32_avg_stride_avx2(src + 32, src_stride, x_offset, y_offset,
                                 ref + 32, ref_stride, sec + 32, 64, h, &sse2,
                                 &sum2);
  *sse = sse1 + sse2;
  *sum = sum1 + sum2;
}

#define HIGHBD_SPV_AVX2(w, h, shift, fn_var, fn_avg)                         \
  uint32_t vpx_highbd_8_sub_pixel_variance##w##x##h##_avx2(                  \
      const uint8_t *src8, int src_stride, int x_offset, int y_offset,       \
      const uint8_t *ref8, int ref_stride, uint32_t *sse) {                  \
    uint64_t sse_long;                                                       \
    int64_t sum_long;                                                        \
    fn_var(CONVERT_TO_SHORTPTR(src8), src_stride, x_offset, y_offset,        \
           CONVERT_TO_SHORTPTR(ref8), ref_stride, NULL, 0, h, &sse_long,     \
           &sum_long);                                                       \
    return finalize_highbd_var_8(sse_long, sum_long, shift, sse);            \
  }                                                                          \
                                                                             \
  uint32_t vpx_highbd_10_sub_pixel_variance##w##x##h##_avx2(                 \
      const uint8_t *src8, int src_stride, int x_offset, int y_offset,       \
      const uint8_t *ref8, int ref_stride, uint32_t *sse) {                  \
    uint64_t sse_long;                                                       \
    int64_t sum_long;                                                        \
    fn_var(CONVERT_TO_SHORTPTR(src8), src_stride, x_offset, y_offset,        \
           CONVERT_TO_SHORTPTR(ref8), ref_stride, NULL, 0, h, &sse_long,     \
           &sum_long);                                                       \
    return finalize_highbd_var_10(sse_long, sum_long, shift, sse);           \
  }                                                                          \
                                                                             \
  uint32_t vpx_highbd_12_sub_pixel_variance##w##x##h##_avx2(                 \
      const uint8_t *src8, int src_stride, int x_offset, int y_offset,       \
      const uint8_t *ref8, int ref_stride, uint32_t *sse) {                  \
    uint64_t sse_long;                                                       \
    int64_t sum_long;                                                        \
    fn_var(CONVERT_TO_SHORTPTR(src8), src_stride, x_offset, y_offset,        \
           CONVERT_TO_SHORTPTR(ref8), ref_stride, NULL, 0, h, &sse_long,     \
           &sum_long);                                                       \
    return finalize_highbd_var_12(sse_long, sum_long, shift, sse);           \
  }                                                                          \
                                                                             \
  uint32_t vpx_highbd_8_sub_pixel_avg_variance##w##x##h##_avx2(              \
      const uint8_t *src8, int src_stride, int x_offset, int y_offset,       \
      const uint8_t *ref8, int ref_stride, uint32_t *sse,                    \
      const uint8_t *sec8) {                                                 \
    uint64_t sse_long;                                                       \
    int64_t sum_long;                                                        \
    fn_avg(CONVERT_TO_SHORTPTR(src8), src_stride, x_offset, y_offset,        \
           CONVERT_TO_SHORTPTR(ref8), ref_stride, CONVERT_TO_SHORTPTR(sec8), \
           1, h, &sse_long, &sum_long);                                      \
    return finalize_highbd_var_8(sse_long, sum_long, shift, sse);            \
  }                                                                          \
                                                                             \
  uint32_t vpx_highbd_10_sub_pixel_avg_variance##w##x##h##_avx2(             \
      const uint8_t *src8, int src_stride, int x_offset, int y_offset,       \
      const uint8_t *ref8, int ref_stride, uint32_t *sse,                    \
      const uint8_t *sec8) {                                                 \
    uint64_t sse_long;                                                       \
    int64_t sum_long;                                                        \
    fn_avg(CONVERT_TO_SHORTPTR(src8), src_stride, x_offset, y_offset,        \
           CONVERT_TO_SHORTPTR(ref8), ref_stride, CONVERT_TO_SHORTPTR(sec8), \
           1, h, &sse_long, &sum_long);                                      \
    return finalize_highbd_var_10(sse_long, sum_long, shift, sse);           \
  }                                                                          \
                                                                             \
  uint32_t vpx_highbd_12_sub_pixel_avg_variance##w##x##h##_avx2(             \
      const uint8_t *src8, int src_stride, int x_offset, int y_offset,       \
      const uint8_t *ref8, int ref_stride, uint32_t *sse,                    \
      const uint8_t *sec8) {                                                 \
    uint64_t sse_long;                                                       \
    int64_t sum_long;                                                        \
    fn_avg(CONVERT_TO_SHORTPTR(src8), src_stride, x_offset, y_offset,        \
           CONVERT_TO_SHORTPTR(ref8), ref_stride, CONVERT_TO_SHORTPTR(sec8), \
           1, h, &sse_long, &sum_long);                                      \
    return finalize_highbd_var_12(sse_long, sum_long, shift, sse);           \
  }

HIGHBD_SPV_AVX2(4, 4, 4, highbd_spv_w4_avx2, highbd_spv_w4_avx2)
HIGHBD_SPV_AVX2(4, 8, 5, highbd_spv_w4_avx2, highbd_spv_w4_avx2)
HIGHBD_SPV_AVX2(8, 4, 5, highbd_spv_w8_var_avx2, highbd_spv_w8_avg_avx2)
HIGHBD_SPV_AVX2(8, 8, 6, highbd_spv_w8_var_avx2, highbd_spv_w8_avg_avx2)
HIGHBD_SPV_AVX2(8, 16, 7, highbd_spv_w8_var_avx2, highbd_spv_w8_avg_avx2)
HIGHBD_SPV_AVX2(16, 8, 7, highbd_spv_w16_var_avx2, highbd_spv_w16_avg_avx2)
HIGHBD_SPV_AVX2(16, 16, 8, highbd_spv_w16_var_avx2, highbd_spv_w16_avg_avx2)
HIGHBD_SPV_AVX2(16, 32, 9, highbd_spv_w16_var_avx2, highbd_spv_w16_avg_avx2)
HIGHBD_SPV_AVX2(32, 16, 9, highbd_spv_w32_var_avx2, highbd_spv_w32_avg_avx2)
HIGHBD_SPV_AVX2(32, 32, 10, highbd_spv_w32_var_avx2, highbd_spv_w32_avg_avx2)
HIGHBD_SPV_AVX2(32, 64, 11, highbd_spv_w32_var_avx2, highbd_spv_w32_avg_avx2)
HIGHBD_SPV_AVX2(64, 32, 11, highbd_spv_w64_var_avx2, highbd_spv_w64_avg_avx2)
HIGHBD_SPV_AVX2(64, 64, 12, highbd_spv_w64_var_avx2, highbd_spv_w64_avg_avx2)

#undef HIGHBD_SPV_AVX2
