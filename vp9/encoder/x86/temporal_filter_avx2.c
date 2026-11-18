/*
 *  Copyright (c) 2024 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <immintrin.h>

#include "./vp9_rtcd.h"
#include "vp9/encoder/vp9_temporal_filter.h"
#include "vp9/encoder/vp9_temporal_filter_constants.h"

DECLARE_ALIGNED(32, static const uint8_t,
                shuffle_src_mask1_avx2[32]) = { 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5,
                                                6, 6, 7, 7, 8, 0, 1, 1, 2, 2, 3,
                                                3, 4, 4, 5, 5, 6, 6, 7, 7, 8 };

DECLARE_ALIGNED(32, static const uint8_t, shuffle_src_mask2_avx2[32]) = {
  2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10,
  2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10
};

DECLARE_ALIGNED(32, static const uint8_t, shuffle_src_mask3_avx2[32]) = {
  4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12,
  4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12
};

DECLARE_ALIGNED(32, static const uint8_t, shuffle_src_mask4_avx2[32]) = {
  6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14,
  6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14
};

static INLINE void shuffle_12tap_filter_avx2(const int16_t *filter,
                                             __m256i *f) {
  const __m256i f_low =
      _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)filter));
  const __m256i f_high = _mm256_broadcastsi128_si256(
      _mm_loadl_epi64((const __m128i *)(filter + 8)));

  f[0] = _mm256_shuffle_epi8(f_low, _mm256_set1_epi16(0x0200u));
  f[1] = _mm256_shuffle_epi8(f_low, _mm256_set1_epi16(0x0604u));
  f[2] = _mm256_shuffle_epi8(f_low, _mm256_set1_epi16(0x0a08u));
  f[3] = _mm256_shuffle_epi8(f_low, _mm256_set1_epi16(0x0e0cu));
  f[4] = _mm256_shuffle_epi8(f_high, _mm256_set1_epi16(0x0200u));
  f[5] = _mm256_shuffle_epi8(f_high, _mm256_set1_epi16(0x0604u));
}

static INLINE void shuffle_src_data_avx2(const __m256i *r1, const __m256i *r2,
                                         const __m256i *f, __m256i *s) {
  s[0] = _mm256_shuffle_epi8(*r1, f[0]);
  s[1] = _mm256_shuffle_epi8(*r1, f[1]);
  s[2] = _mm256_shuffle_epi8(*r1, f[2]);
  s[3] = _mm256_shuffle_epi8(*r1, f[3]);
  s[4] = _mm256_shuffle_epi8(*r2, f[0]);
  s[5] = _mm256_shuffle_epi8(*r2, f[1]);
}

static INLINE void reuse_src_data_avx2(const __m256i *src, __m256i *des) {
  des[0] = src[0];
  des[1] = src[1];
  des[2] = src[2];
  des[3] = src[3];
  des[4] = src[4];
}

static INLINE __m256i convolve12_16_avx2(const __m256i *s, const __m256i *f) {
  // multiply 2 adjacent elements with the filter and add the result
  const __m256i k_64 = _mm256_set1_epi16(1 << (FILTER_BITS - 1));
  const __m256i x0 = _mm256_maddubs_epi16(s[0], f[0]);
  const __m256i x1 = _mm256_maddubs_epi16(s[1], f[1]);
  const __m256i x2 = _mm256_maddubs_epi16(s[2], f[2]);
  const __m256i x3 = _mm256_maddubs_epi16(s[3], f[3]);
  const __m256i x4 = _mm256_maddubs_epi16(s[4], f[4]);
  const __m256i x5 = _mm256_maddubs_epi16(s[5], f[5]);
  __m256i sum1, sum2, sum3;

  sum1 = _mm256_add_epi16(x0, x2);
  sum2 = _mm256_add_epi16(x3, x5);
  sum3 = _mm256_add_epi16(x1, x4);
  sum3 = _mm256_add_epi16(sum3, k_64);

  const __m256i s0 = _mm256_cvtepi16_epi32(_mm256_castsi256_si128(sum1));
  const __m256i s1 = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(sum1, 1));
  const __m256i s2 = _mm256_cvtepi16_epi32(_mm256_castsi256_si128(sum2));
  const __m256i s3 = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(sum2, 1));
  const __m256i s4 = _mm256_cvtepi16_epi32(_mm256_castsi256_si128(sum3));
  const __m256i s5 = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(sum3, 1));

  sum1 = _mm256_add_epi32(s0, s2);
  sum2 = _mm256_add_epi32(s1, s3);
  sum1 = _mm256_add_epi32(sum1, s4);
  sum2 = _mm256_add_epi32(sum2, s5);

  // round and shift by 7 bit each 32 bit
  // 0 1 2 3 4 5 6 7
  sum1 = _mm256_srai_epi32(sum1, FILTER_BITS);
  // 8 9 10 11 12 13 14 15
  sum2 = _mm256_srai_epi32(sum2, FILTER_BITS);

  // 0 1 2 3 8 9 10 11 4 5 6 7 12 13 14 15
  // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
  __m256i const res =
      _mm256_permute4x64_epi64(_mm256_packus_epi32(sum1, sum2), 0xD8);
  return res;
}

void vpx_convolve12_horiz_avx2(const uint8_t *src, ptrdiff_t src_stride,
                               uint8_t *dst, ptrdiff_t dst_stride,
                               const InterpKernel12 *filter, int x0_q4,
                               int x_step_q4, int y0_q4, int y_step_q4, int w,
                               int h) {
  assert(x_step_q4 == 16);
  assert(w == 32 || w == 16 || w == 8);
  (void)y0_q4;
  (void)x_step_q4;
  (void)y_step_q4;
  const uint8_t *src_ptr = src;
  src_ptr -= MAX_FILTER_TAP / 2 - 1;
  __m256i s[6], f[6], src_mask[4];

  shuffle_12tap_filter_avx2(filter[x0_q4], f);
  src_mask[0] = _mm256_load_si256((__m256i const *)shuffle_src_mask1_avx2);
  src_mask[1] = _mm256_load_si256((__m256i const *)shuffle_src_mask2_avx2);
  src_mask[2] = _mm256_load_si256((__m256i const *)shuffle_src_mask3_avx2);
  src_mask[3] = _mm256_load_si256((__m256i const *)shuffle_src_mask4_avx2);
  if (w == 8) {
    for (int i = 0; i < h; i += 4) {
      // s00 s01 s02 s03 s04 s05 s06 s07 s08 s09 s010 s011 s012 s013 s014 s015
      const __m128i row0 =
          _mm_loadu_si128((const __m128i *)&src_ptr[i * src_stride]);
      // s08 s09 s010 s011 s012 s013 s014 s015 s016 s017 s018 s019 s020 s021
      // s022 s023
      const __m128i row0_8 =
          _mm_loadu_si128((const __m128i *)&src_ptr[i * src_stride + 8]);
      // s10 s11 s12 s13 s14 s15 s16 s17 s18 s19 s110 s111 s112 s113 s114 s115
      const __m128i row1 =
          _mm_loadu_si128((const __m128i *)&src_ptr[(i + 1) * src_stride]);
      const __m128i row1_8 =
          _mm_loadu_si128((const __m128i *)&src_ptr[(i + 1) * src_stride + 8]);
      // s20 s21 s22 s23 s24 s25 s26 s27 s28 s29 s210 s211 s212 s213 s214 s215
      const __m128i row2 =
          _mm_loadu_si128((const __m128i *)&src_ptr[(i + 2) * src_stride]);
      const __m128i row2_8 =
          _mm_loadu_si128((const __m128i *)&src_ptr[(i + 2) * src_stride + 8]);
      // s30 s31 s32 s33 s34 s35 s36 s37 s38 s39 s310 s311 s312 s313 s314 s115
      const __m128i row3 =
          _mm_loadu_si128((const __m128i *)&src_ptr[(i + 3) * src_stride]);
      const __m128i row3_8 =
          _mm_loadu_si128((const __m128i *)&src_ptr[(i + 3) * src_stride + 8]);
      // s00 s01 s02 s03 s04 s05 s06 s07 s08 s09 s010 s011 s012 s013 s014 s015 |
      // s10 s11 s12 s13 s14 s15 s16 s17 s18 s19 s110 s111 s112 s113 s114 s115
      const __m256i row01 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(row0), row1, 1);
      // s20 s21 s22 s23 s24 s25 s26 s27 s28 s29 s210 s211 s212 s213 s214 s215 |
      // s30 s31 s32 s33 s34 s35 s36 s37 s38 s39 s310 s311 s312 s313 s314 s115
      const __m256i row23 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(row2), row3, 1);
      // s08 s09 s010 s011 s012 s013 s014 s015 s016 s017 s018 s019 s020 s021
      // s022 s023 | s18 s19 s110 s111 s112 s113 s114 s115 s116 s117 s118 s119
      // s120 s121 s122 s123
      const __m256i row01_8 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(row0_8), row1_8, 1);
      const __m256i row23_8 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(row2_8), row3_8, 1);

      shuffle_src_data_avx2(&row01, &row01_8, src_mask, s);
      const __m256i res_0 = convolve12_16_avx2(s, f);

      shuffle_src_data_avx2(&row23, &row23_8, src_mask, s);
      const __m256i res_1 = convolve12_16_avx2(s, f);

      // 00 01 02 03 04 05 06 07 | 10 11 12 13 14 15 16 17 | 08 09 010 011 012
      // 013 014 015 | 18 19 110 111 112 113 114 115
      const __m256i res = _mm256_packus_epi16(res_0, res_1);
      const __m128i res_lo = _mm256_castsi256_si128(res);
      const __m128i res_hi = _mm256_extracti128_si256(res, 1);

      _mm_storel_epi64((__m128i *)&dst[i * dst_stride], res_lo);
      _mm_storel_epi64((__m128i *)&dst[(i + 1) * dst_stride], res_hi);
      _mm_storel_epi64((__m128i *)&dst[(i + 2) * dst_stride],
                       _mm_srli_si128(res_lo, 8));
      _mm_storel_epi64((__m128i *)&dst[(i + 3) * dst_stride],
                       _mm_srli_si128(res_hi, 8));
    }
  } else {
    for (int j = 0; j < w; j += 16) {
      for (int i = 0; i < h; i += 2) {
        // s00 s01 s02 s03 s04 s05 s06 s07 s08 s09 s010 s011 s012 s013 s014 s015
        const __m128i row0 =
            _mm_loadu_si128((const __m128i *)&src_ptr[i * src_stride + j]);
        // s016 s017 s018 s019 s020 s021 s022 s023 s024 s025 s026 s027 s028 s029
        // s030 s031
        const __m128i row0_16 =
            _mm_loadu_si128((const __m128i *)&src_ptr[i * src_stride + j + 16]);
        // s10 s11 s12 s13 s14 s15 s16 s17 s18 s19 s110 s111 s112 s113 s114
        // s115
        const __m128i row1 = _mm_loadu_si128(
            (const __m128i *)&src_ptr[(i + 1) * src_stride + j]);
        // s116 s117 s118 s119 s120 s121 s122 s123 s124 s125 s126 s127 s128
        // s129 s130 s131
        const __m128i row1_16 = _mm_loadu_si128(
            (const __m128i *)&src_ptr[(i + 1) * src_stride + j + 16]);

        // s00 s01 s02 s03 s04 s05 s06 s07 s08 s09 s010 s011 s012 s013 s014 s015
        // | s10 s11 s12 s13 s14 s15 s16 s17 s18 s19 s110 s111 s112 s113 s114
        // s115
        const __m256i r0 =
            _mm256_inserti128_si256(_mm256_castsi128_si256(row0), row1, 1);
        // s016 s017 s018 s019 s020 s021 s022 s023 s024 s025 s026 s027 s028 s029
        // s030 s031 | s116 s117 s118 s119 s120 s121 s122 s123 s124 s125 s126
        // s127 s128 s129 s130 s131
        const __m256i r2 = _mm256_inserti128_si256(
            _mm256_castsi128_si256(row0_16), row1_16, 1);

        // s08 s09 s010 s011 s012 s013 s014 s015 s016 s017 s018 s019 s020 s021
        // s022 s023 | s18 s19 s110 s111 s112 s113 s114 s115 s116 s117 s118 s119
        // s120 s121 s122 s123
        const __m256i r1 = _mm256_alignr_epi8(r2, r0, 8);

        shuffle_src_data_avx2(&r0, &r1, src_mask, s);
        const __m256i res_0 = convolve12_16_avx2(s, f);

        shuffle_src_data_avx2(&r1, &r2, src_mask, s);
        const __m256i res_1 = convolve12_16_avx2(s, f);

        const __m256i res = _mm256_packus_epi16(res_0, res_1);

        _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j],
                         _mm256_castsi256_si128(res));
        if (i + 1 < h) {
          _mm_storeu_si128((__m128i *)&dst[(i + 1) * dst_stride + j],
                           _mm256_extracti128_si256(res, 1));
        }
      }
    }
  }
}

void vpx_convolve12_vert_avx2(const uint8_t *src, ptrdiff_t src_stride,
                              uint8_t *dst, ptrdiff_t dst_stride,
                              const InterpKernel12 *filter, int x0_q4,
                              int x_step_q4, int y0_q4, int y_step_q4, int w,
                              int h) {
  assert(y_step_q4 == 16);
  assert(h == 32 || h == 16 || h == 8);
  assert(w == 32 || w == 16 || w == 8);
  (void)x0_q4;
  (void)x_step_q4;
  (void)y_step_q4;
  const uint8_t *src_ptr = src;
  src_ptr -= src_stride * (MAX_FILTER_TAP / 2 - 1);
  __m256i s[12], f[6];

  shuffle_12tap_filter_avx2(filter[y0_q4], f);
  if (w == 8) {
    const __m128i s0 =
        _mm_loadl_epi64((const __m128i *)(src_ptr + 0 * src_stride));
    const __m128i s1 =
        _mm_loadl_epi64((const __m128i *)(src_ptr + 1 * src_stride));
    const __m128i s2 =
        _mm_loadl_epi64((const __m128i *)(src_ptr + 2 * src_stride));
    const __m128i s3 =
        _mm_loadl_epi64((const __m128i *)(src_ptr + 3 * src_stride));
    const __m128i s4 =
        _mm_loadl_epi64((const __m128i *)(src_ptr + 4 * src_stride));
    const __m128i s5 =
        _mm_loadl_epi64((const __m128i *)(src_ptr + 5 * src_stride));
    const __m128i s6 =
        _mm_loadl_epi64((const __m128i *)(src_ptr + 6 * src_stride));
    const __m128i s7 =
        _mm_loadl_epi64((const __m128i *)(src_ptr + 7 * src_stride));
    const __m128i s8 =
        _mm_loadl_epi64((const __m128i *)(src_ptr + 8 * src_stride));
    const __m128i s9 =
        _mm_loadl_epi64((const __m128i *)(src_ptr + 9 * src_stride));
    const __m128i s10t =
        _mm_loadl_epi64((const __m128i *)(src_ptr + 10 * src_stride));

    const __m256i r01 =
        _mm256_inserti128_si256(_mm256_castsi128_si256(s0), s1, 1);
    const __m256i r12 =
        _mm256_inserti128_si256(_mm256_castsi128_si256(s1), s2, 1);
    const __m256i r23 =
        _mm256_inserti128_si256(_mm256_castsi128_si256(s2), s3, 1);
    const __m256i r34 =
        _mm256_inserti128_si256(_mm256_castsi128_si256(s3), s4, 1);
    const __m256i r45 =
        _mm256_inserti128_si256(_mm256_castsi128_si256(s4), s5, 1);
    const __m256i r56 =
        _mm256_inserti128_si256(_mm256_castsi128_si256(s5), s6, 1);
    const __m256i r67 =
        _mm256_inserti128_si256(_mm256_castsi128_si256(s6), s7, 1);
    const __m256i r78 =
        _mm256_inserti128_si256(_mm256_castsi128_si256(s7), s8, 1);
    const __m256i r89 =
        _mm256_inserti128_si256(_mm256_castsi128_si256(s8), s9, 1);
    const __m256i r910 =
        _mm256_inserti128_si256(_mm256_castsi128_si256(s9), s10t, 1);

    s[0] = _mm256_unpacklo_epi8(r01, r12);
    s[1] = _mm256_unpacklo_epi8(r23, r34);
    s[2] = _mm256_unpacklo_epi8(r45, r56);
    s[3] = _mm256_unpacklo_epi8(r67, r78);
    s[4] = _mm256_unpacklo_epi8(r89, r910);
    for (int i = 0; i < h; i += 2) {
      const __m128i s10 =
          _mm_loadl_epi64((const __m128i *)(src_ptr + (i + 10) * src_stride));
      const __m128i s11 =
          _mm_loadl_epi64((const __m128i *)(src_ptr + (i + 11) * src_stride));
      const __m128i s12 =
          _mm_loadl_epi64((const __m128i *)(src_ptr + (i + 12) * src_stride));

      const __m256i r1011 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s10), s11, 1);
      const __m256i r1112 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s11), s12, 1);
      s[5] = _mm256_unpacklo_epi8(r1011, r1112);
      const __m256i res_0 = convolve12_16_avx2(s, f);

      __m256i res = _mm256_packus_epi16(res_0, res_0);

      _mm_storel_epi64((__m128i *)&dst[i * dst_stride],
                       _mm256_castsi256_si128(res));
      _mm_storel_epi64((__m128i *)&dst[(i + 1) * dst_stride],
                       _mm256_extracti128_si256(res, 1));

      reuse_src_data_avx2(s + 1, s);
    }
  } else {
    for (int j = 0; j < w; j += 16) {
      const __m128i s0 =
          _mm_loadu_si128((const __m128i *)(src_ptr + 0 * src_stride + j));
      const __m128i s1 =
          _mm_loadu_si128((const __m128i *)(src_ptr + 1 * src_stride + j));
      const __m128i s2 =
          _mm_loadu_si128((const __m128i *)(src_ptr + 2 * src_stride + j));
      const __m128i s3 =
          _mm_loadu_si128((const __m128i *)(src_ptr + 3 * src_stride + j));
      const __m128i s4 =
          _mm_loadu_si128((const __m128i *)(src_ptr + 4 * src_stride + j));
      const __m128i s5 =
          _mm_loadu_si128((const __m128i *)(src_ptr + 5 * src_stride + j));
      const __m128i s6 =
          _mm_loadu_si128((const __m128i *)(src_ptr + 6 * src_stride + j));
      const __m128i s7 =
          _mm_loadu_si128((const __m128i *)(src_ptr + 7 * src_stride + j));
      const __m128i s8 =
          _mm_loadu_si128((const __m128i *)(src_ptr + 8 * src_stride + j));
      const __m128i s9 =
          _mm_loadu_si128((const __m128i *)(src_ptr + 9 * src_stride + j));
      const __m128i s10t =
          _mm_loadu_si128((const __m128i *)(src_ptr + 10 * src_stride + j));

      const __m256i r01 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s0), s1, 1);
      const __m256i r12 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s1), s2, 1);
      const __m256i r23 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s2), s3, 1);
      const __m256i r34 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s3), s4, 1);
      const __m256i r45 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s4), s5, 1);
      const __m256i r56 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s5), s6, 1);
      const __m256i r67 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s6), s7, 1);
      const __m256i r78 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s7), s8, 1);
      const __m256i r89 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s8), s9, 1);
      const __m256i r910 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s9), s10t, 1);

      s[0] = _mm256_unpacklo_epi8(r01, r12);
      s[1] = _mm256_unpacklo_epi8(r23, r34);
      s[2] = _mm256_unpacklo_epi8(r45, r56);
      s[3] = _mm256_unpacklo_epi8(r67, r78);
      s[4] = _mm256_unpacklo_epi8(r89, r910);

      s[6] = _mm256_unpackhi_epi8(r01, r12);
      s[7] = _mm256_unpackhi_epi8(r23, r34);
      s[8] = _mm256_unpackhi_epi8(r45, r56);
      s[9] = _mm256_unpackhi_epi8(r67, r78);
      s[10] = _mm256_unpackhi_epi8(r89, r910);
      for (int i = 0; i < h; i += 2) {
        const __m128i s10 = _mm_loadu_si128(
            (const __m128i *)(src_ptr + (i + 10) * src_stride + j));
        const __m128i s11 = _mm_loadu_si128(
            (const __m128i *)(src_ptr + (i + 11) * src_stride + j));
        const __m128i s12 = _mm_loadu_si128(
            (const __m128i *)(src_ptr + (i + 12) * src_stride + j));

        const __m256i r1011 =
            _mm256_inserti128_si256(_mm256_castsi128_si256(s10), s11, 1);
        const __m256i r1112 =
            _mm256_inserti128_si256(_mm256_castsi128_si256(s11), s12, 1);

        s[5] = _mm256_unpacklo_epi8(r1011, r1112);
        s[11] = _mm256_unpackhi_epi8(r1011, r1112);

        const __m256i res_0 = convolve12_16_avx2(s, f);
        const __m256i res_1 = convolve12_16_avx2(s + 6, f);

        __m256i res = _mm256_packus_epi16(res_0, res_1);

        _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j],
                         _mm256_castsi256_si128(res));
        _mm_storeu_si128((__m128i *)&dst[(i + 1) * dst_stride + j],
                         _mm256_extracti128_si256(res, 1));

        reuse_src_data_avx2(s + 1, s);
        reuse_src_data_avx2(s + 7, s + 6);
      }
    }
  }
}

void vpx_convolve12_avx2(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
                         ptrdiff_t dst_stride, const InterpKernel12 *filter,
                         int x0_q4, int x_step_q4, int y0_q4, int y_step_q4,
                         int w, int h) {
  assert(x_step_q4 == 16 && y_step_q4 == 16);
  assert(h == 32 || h == 16 || h == 8);
  assert(w == 32 || w == 16 || w == 8);
  DECLARE_ALIGNED(32, uint8_t, temp[BW * (BH + MAX_FILTER_TAP - 1)]);
  const int temp_stride = BW;
  const int intermediate_height =
      (((h - 1) * y_step_q4 + y0_q4) >> SUBPEL_BITS) + MAX_FILTER_TAP;
  vpx_convolve12_horiz_avx2(src - src_stride * (MAX_FILTER_TAP / 2 - 1),
                            src_stride, temp, temp_stride, filter, x0_q4,
                            x_step_q4, y0_q4, y_step_q4, w,
                            intermediate_height);
  vpx_convolve12_vert_avx2(temp + temp_stride * (MAX_FILTER_TAP / 2 - 1),
                           temp_stride, dst, dst_stride, filter, x0_q4,
                           x_step_q4, y0_q4, y_step_q4, w, h);
}

static INLINE __m256i combine_128_pair(const __m128i lo, const __m128i hi) {
  return _mm256_inserti128_si256(_mm256_castsi128_si256(lo), hi, 1);
}

static INLINE void store_dist_8_avx2(const uint8_t *a, const uint8_t *b,
                                     uint16_t *dst) {
  const __m128i a_reg = _mm_loadl_epi64((const __m128i *)a);
  const __m128i b_reg = _mm_loadl_epi64((const __m128i *)b);
  const __m128i diff =
      _mm_sub_epi16(_mm_cvtepu8_epi16(a_reg), _mm_cvtepu8_epi16(b_reg));
  _mm_storeu_si128((__m128i *)dst, _mm_mullo_epi16(diff, diff));
}

static INLINE void store_dist_16_avx2(const uint8_t *a, const uint8_t *b,
                                      uint16_t *dst) {
  const __m128i a_reg = _mm_loadu_si128((const __m128i *)a);
  const __m128i b_reg = _mm_loadu_si128((const __m128i *)b);
  const __m256i diff = _mm256_sub_epi16(_mm256_cvtepu8_epi16(a_reg),
                                        _mm256_cvtepu8_epi16(b_reg));
  _mm256_storeu_si256((__m256i *)dst, _mm256_mullo_epi16(diff, diff));
}

static INLINE void store_dist_32_avx2(const uint8_t *a, const uint8_t *b,
                                      uint16_t *dst) {
  const __m256i a_reg = _mm256_loadu_si256((const __m256i *)a);
  const __m256i b_reg = _mm256_loadu_si256((const __m256i *)b);
  const __m256i diff_lo =
      _mm256_sub_epi16(_mm256_cvtepu8_epi16(_mm256_castsi256_si128(a_reg)),
                       _mm256_cvtepu8_epi16(_mm256_castsi256_si128(b_reg)));
  const __m256i diff_hi = _mm256_sub_epi16(
      _mm256_cvtepu8_epi16(_mm256_extracti128_si256(a_reg, 1)),
      _mm256_cvtepu8_epi16(_mm256_extracti128_si256(b_reg, 1)));
  _mm256_storeu_si256((__m256i *)dst, _mm256_mullo_epi16(diff_lo, diff_lo));
  _mm256_storeu_si256((__m256i *)(dst + 16),
                      _mm256_mullo_epi16(diff_hi, diff_hi));
}

static INLINE __m256i get_sum_16_avx2(const uint16_t *dist) {
  const __m256i dist_reg = _mm256_loadu_si256((const __m256i *)dist);
  const __m256i dist_left = _mm256_loadu_si256((const __m256i *)(dist - 1));
  const __m256i dist_right = _mm256_loadu_si256((const __m256i *)(dist + 1));
  return _mm256_adds_epu16(_mm256_adds_epu16(dist_reg, dist_left), dist_right);
}

static INLINE __m128i get_sum_8_sse4(const uint16_t *dist) {
  const __m128i dist_reg = _mm_loadu_si128((const __m128i *)dist);
  const __m128i dist_left = _mm_loadu_si128((const __m128i *)(dist - 1));
  const __m128i dist_right = _mm_loadu_si128((const __m128i *)(dist + 1));
  return _mm_adds_epu16(_mm_adds_epu16(dist_reg, dist_left), dist_right);
}

static INLINE void read_chroma_dist_row_16_avx2(int ss_x,
                                                const uint16_t *u_dist,
                                                const uint16_t *v_dist,
                                                __m256i *u_out,
                                                __m256i *v_out) {
  if (!ss_x) {
    *u_out = _mm256_loadu_si256((const __m256i *)u_dist);
    *v_out = _mm256_loadu_si256((const __m256i *)v_dist);
  } else {
    const __m128i u_reg = _mm_loadu_si128((const __m128i *)u_dist);
    const __m128i v_reg = _mm_loadu_si128((const __m128i *)v_dist);
    *u_out = combine_128_pair(_mm_unpacklo_epi16(u_reg, u_reg),
                              _mm_unpackhi_epi16(u_reg, u_reg));
    *v_out = combine_128_pair(_mm_unpacklo_epi16(v_reg, v_reg),
                              _mm_unpackhi_epi16(v_reg, v_reg));
  }
}

static INLINE __m256i average_16_avx2(__m256i sum, const __m256i *mul_constants,
                                      int strength, int rounding,
                                      const __m256i *weight) {
  const __m128i strength_u128 = _mm_set_epi32(0, 0, 0, strength);
  const __m256i rounding_u16 = _mm256_set1_epi16(rounding);
  const __m256i sixteen = _mm256_set1_epi16(16);

  sum = _mm256_mulhi_epu16(sum, *mul_constants);
  sum = _mm256_adds_epu16(sum, rounding_u16);
  sum = _mm256_srl_epi16(sum, strength_u128);
  sum = _mm256_min_epu16(sum, sixteen);
  sum = _mm256_sub_epi16(sixteen, sum);
  return _mm256_mullo_epi16(sum, *weight);
}

static INLINE __m128i average_8_sse4(__m128i sum, const __m128i *mul_constants,
                                     int strength, int rounding,
                                     const __m128i *weight) {
  const __m128i strength_u128 = _mm_set_epi32(0, 0, 0, strength);
  const __m128i rounding_u16 = _mm_set1_epi16(rounding);
  const __m128i sixteen = _mm_set1_epi16(16);

  sum = _mm_mulhi_epu16(sum, *mul_constants);
  sum = _mm_adds_epu16(sum, rounding_u16);
  sum = _mm_srl_epi16(sum, strength_u128);
  sum = _mm_min_epu16(sum, sixteen);
  sum = _mm_sub_epi16(sixteen, sum);
  return _mm_mullo_epi16(sum, *weight);
}

static INLINE void accumulate_and_store_16_avx2(const __m256i sum_u16,
                                                const uint8_t *pred,
                                                uint16_t *count,
                                                uint32_t *accumulator) {
  const __m128i pred_u8 = _mm_loadu_si128((const __m128i *)pred);
  __m256i count_u16 = _mm256_loadu_si256((const __m256i *)count);
  const __m256i pred_u16 = _mm256_cvtepu8_epi16(pred_u8);
  __m256i prod_u16, prod_0_u32, prod_1_u32, accum_0_u32, accum_1_u32;

  count_u16 = _mm256_adds_epu16(count_u16, sum_u16);
  _mm256_storeu_si256((__m256i *)count, count_u16);

  prod_u16 = _mm256_mullo_epi16(sum_u16, pred_u16);
  prod_0_u32 = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(prod_u16));
  prod_1_u32 = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(prod_u16, 1));

  accum_0_u32 = _mm256_loadu_si256((const __m256i *)accumulator);
  accum_1_u32 = _mm256_loadu_si256((const __m256i *)(accumulator + 8));

  accum_0_u32 = _mm256_add_epi32(accum_0_u32, prod_0_u32);
  accum_1_u32 = _mm256_add_epi32(accum_1_u32, prod_1_u32);

  _mm256_storeu_si256((__m256i *)accumulator, accum_0_u32);
  _mm256_storeu_si256((__m256i *)(accumulator + 8), accum_1_u32);
}

static INLINE void accumulate_and_store_8_sse4(const __m128i sum_u16,
                                               const uint8_t *pred,
                                               uint16_t *count,
                                               uint32_t *accumulator) {
  const __m128i pred_u8 = _mm_loadl_epi64((const __m128i *)pred);
  __m128i count_u16 = _mm_loadu_si128((const __m128i *)count);
  const __m128i pred_u16 = _mm_cvtepu8_epi16(pred_u8);
  const __m128i prod_u16 = _mm_mullo_epi16(sum_u16, pred_u16);
  const __m256i prod_u32 = _mm256_cvtepu16_epi32(prod_u16);
  const __m256i accum_u32 = _mm256_loadu_si256((const __m256i *)accumulator);

  _mm_storeu_si128((__m128i *)count, _mm_adds_epu16(count_u16, sum_u16));
  _mm256_storeu_si256((__m256i *)accumulator,
                      _mm256_add_epi32(accum_u32, prod_u32));
}

static INLINE __m256i hadd_epu16_to_u16_avx2(const __m256i a, const __m256i b) {
  const __m256i low_mask = _mm256_set1_epi32(0x0000FFFF);
  const __m256i a_even = _mm256_and_si256(a, low_mask);
  const __m256i a_odd = _mm256_srli_epi32(a, 16);
  const __m256i b_even = _mm256_and_si256(b, low_mask);
  const __m256i b_odd = _mm256_srli_epi32(b, 16);
  const __m256i sum_a = _mm256_add_epi32(a_even, a_odd);
  const __m256i sum_b = _mm256_add_epi32(b_even, b_odd);
  return _mm256_permute4x64_epi64(_mm256_packus_epi32(sum_a, sum_b), 0xD8);
}

static INLINE void add_luma_dist_to_16_chroma_mod_avx2(const uint16_t *y_dist,
                                                       int ss_x, int ss_y,
                                                       __m256i *u_mod,
                                                       __m256i *v_mod) {
  __m256i y_reg;
  if (!ss_x) {
    y_reg = _mm256_loadu_si256((const __m256i *)y_dist);
    if (ss_y == 1) {
      y_reg = _mm256_adds_epu16(
          y_reg, _mm256_loadu_si256((const __m256i *)(y_dist + DIST_STRIDE)));
    }
  } else {
    __m256i y_first = _mm256_loadu_si256((const __m256i *)y_dist);
    __m256i y_second = _mm256_loadu_si256((const __m256i *)(y_dist + 16));
    if (ss_y == 1) {
      y_first = _mm256_adds_epu16(
          y_first, _mm256_loadu_si256((const __m256i *)(y_dist + DIST_STRIDE)));
      y_second = _mm256_adds_epu16(
          y_second,
          _mm256_loadu_si256((const __m256i *)(y_dist + 16 + DIST_STRIDE)));
    }
    y_reg = hadd_epu16_to_u16_avx2(y_first, y_second);
  }

  *u_mod = _mm256_adds_epu16(*u_mod, y_reg);
  *v_mod = _mm256_adds_epu16(*v_mod, y_reg);
}

static INLINE void add_luma_dist_to_8_chroma_mod_sse4(const uint16_t *y_dist,
                                                      int ss_x, int ss_y,
                                                      __m128i *u_mod,
                                                      __m128i *v_mod) {
  __m128i y_reg;
  if (!ss_x) {
    y_reg = _mm_loadu_si128((const __m128i *)y_dist);
    if (ss_y == 1) {
      y_reg = _mm_adds_epu16(
          y_reg, _mm_loadu_si128((const __m128i *)(y_dist + DIST_STRIDE)));
    }
  } else {
    const __m128i low_mask = _mm_set1_epi32(0x0000FFFF);
    __m128i y_first = _mm_loadu_si128((const __m128i *)y_dist);
    __m128i y_second = _mm_loadu_si128((const __m128i *)(y_dist + 8));
    if (ss_y == 1) {
      y_first = _mm_adds_epu16(
          y_first, _mm_loadu_si128((const __m128i *)(y_dist + DIST_STRIDE)));
      y_second = _mm_adds_epu16(
          y_second,
          _mm_loadu_si128((const __m128i *)(y_dist + 8 + DIST_STRIDE)));
    }
    y_first = _mm_add_epi32(_mm_and_si128(y_first, low_mask),
                            _mm_srli_epi32(y_first, 16));
    y_second = _mm_add_epi32(_mm_and_si128(y_second, low_mask),
                             _mm_srli_epi32(y_second, 16));
    y_reg = _mm_packus_epi32(y_first, y_second);
  }
  *u_mod = _mm_adds_epu16(*u_mod, y_reg);
  *v_mod = _mm_adds_epu16(*v_mod, y_reg);
}

static void vp9_apply_temporal_filter_luma_16_avx2(
    const uint8_t *y_pre, int y_pre_stride, unsigned int block_height, int ss_x,
    int ss_y, int strength, int use_whole_blk, uint32_t *y_accum,
    uint16_t *y_count, const uint16_t *y_dist, const uint16_t *u_dist,
    const uint16_t *v_dist, const int16_t *const *neighbors_first,
    const int16_t *const *neighbors_second, int top_weight, int bottom_weight,
    const int *blk_fw) {
  const int rounding = (1 << strength) >> 1;
  __m256i weight, mul;
  __m256i sum_row_1, sum_row_2, sum_row_3;
  __m256i u_dist_reg, v_dist_reg, sum_row;
  unsigned int h;

  if (blk_fw) {
    weight =
        combine_128_pair(_mm_set1_epi16(blk_fw[0]), _mm_set1_epi16(blk_fw[1]));
  } else {
    weight = _mm256_set1_epi16(top_weight);
  }

  mul = combine_128_pair(_mm_load_si128((const __m128i *)neighbors_first[0]),
                         _mm_load_si128((const __m128i *)neighbors_second[0]));

  sum_row_2 = get_sum_16_avx2(y_dist);
  sum_row_3 = get_sum_16_avx2(y_dist + DIST_STRIDE);
  sum_row = _mm256_adds_epu16(sum_row_2, sum_row_3);

  read_chroma_dist_row_16_avx2(ss_x, u_dist, v_dist, &u_dist_reg, &v_dist_reg);
  sum_row = _mm256_adds_epu16(sum_row, u_dist_reg);
  sum_row = _mm256_adds_epu16(sum_row, v_dist_reg);

  sum_row = average_16_avx2(sum_row, &mul, strength, rounding, &weight);
  accumulate_and_store_16_avx2(sum_row, y_pre, y_count, y_accum);

  y_pre += y_pre_stride;
  y_count += y_pre_stride;
  y_accum += y_pre_stride;
  y_dist += DIST_STRIDE;
  u_dist += DIST_STRIDE;
  v_dist += DIST_STRIDE;

  mul = combine_128_pair(_mm_load_si128((const __m128i *)neighbors_first[1]),
                         _mm_load_si128((const __m128i *)neighbors_second[1]));

  for (h = 1; h < block_height - 1; ++h) {
    if (!use_whole_blk && h == block_height / 2) {
      if (blk_fw) {
        weight = combine_128_pair(_mm_set1_epi16(blk_fw[2]),
                                  _mm_set1_epi16(blk_fw[3]));
      } else {
        weight = _mm256_set1_epi16(bottom_weight);
      }
    }
    sum_row_1 = sum_row_2;
    sum_row_2 = sum_row_3;
    sum_row = _mm256_adds_epu16(sum_row_1, sum_row_2);
    sum_row_3 = get_sum_16_avx2(y_dist + DIST_STRIDE);
    sum_row = _mm256_adds_epu16(sum_row, sum_row_3);

    if (ss_y == 0 || h % 2 == 0) {
      read_chroma_dist_row_16_avx2(ss_x, u_dist, v_dist, &u_dist_reg,
                                   &v_dist_reg);
      u_dist += DIST_STRIDE;
      v_dist += DIST_STRIDE;
    }

    sum_row = _mm256_adds_epu16(sum_row, u_dist_reg);
    sum_row = _mm256_adds_epu16(sum_row, v_dist_reg);

    sum_row = average_16_avx2(sum_row, &mul, strength, rounding, &weight);
    accumulate_and_store_16_avx2(sum_row, y_pre, y_count, y_accum);

    y_pre += y_pre_stride;
    y_count += y_pre_stride;
    y_accum += y_pre_stride;
    y_dist += DIST_STRIDE;
  }

  mul = combine_128_pair(_mm_load_si128((const __m128i *)neighbors_first[0]),
                         _mm_load_si128((const __m128i *)neighbors_second[0]));

  sum_row_1 = sum_row_2;
  sum_row_2 = sum_row_3;
  sum_row = _mm256_adds_epu16(sum_row_1, sum_row_2);

  if (ss_y == 0) {
    read_chroma_dist_row_16_avx2(ss_x, u_dist, v_dist, &u_dist_reg,
                                 &v_dist_reg);
  }

  sum_row = _mm256_adds_epu16(sum_row, u_dist_reg);
  sum_row = _mm256_adds_epu16(sum_row, v_dist_reg);

  sum_row = average_16_avx2(sum_row, &mul, strength, rounding, &weight);
  accumulate_and_store_16_avx2(sum_row, y_pre, y_count, y_accum);
}

static void vp9_apply_temporal_filter_luma_avx2(
    const uint8_t *y_pre, int y_pre_stride, unsigned int block_width,
    unsigned int block_height, int ss_x, int ss_y, int strength,
    const int *blk_fw, int use_whole_blk, uint32_t *y_accum, uint16_t *y_count,
    const uint16_t *y_dist, const uint16_t *u_dist, const uint16_t *v_dist) {
  unsigned int blk_col = 0, uv_blk_col = 0;
  const unsigned int blk_col_step = 16, uv_blk_col_step = 16 >> ss_x;
  const unsigned int mid_width = block_width >> 1,
                     last_width = block_width - blk_col_step;
  int top_weight = blk_fw[0],
      bottom_weight = use_whole_blk ? blk_fw[0] : blk_fw[2];
  const int16_t *const *neighbors_first;
  const int16_t *const *neighbors_second;

  if (block_width == 16) {
    neighbors_first = LUMA_LEFT_COLUMN_NEIGHBORS;
    neighbors_second = LUMA_RIGHT_COLUMN_NEIGHBORS;
    vp9_apply_temporal_filter_luma_16_avx2(
        y_pre, y_pre_stride, block_height, ss_x, ss_y, strength, use_whole_blk,
        y_accum, y_count, y_dist, u_dist, v_dist, neighbors_first,
        neighbors_second, top_weight, bottom_weight,
        use_whole_blk ? NULL : blk_fw);
    return;
  }

  neighbors_first = LUMA_LEFT_COLUMN_NEIGHBORS;
  neighbors_second = LUMA_MIDDLE_COLUMN_NEIGHBORS;
  vp9_apply_temporal_filter_luma_16_avx2(
      y_pre + blk_col, y_pre_stride, block_height, ss_x, ss_y, strength,
      use_whole_blk, y_accum + blk_col, y_count + blk_col, y_dist + blk_col,
      u_dist + uv_blk_col, v_dist + uv_blk_col, neighbors_first,
      neighbors_second, top_weight, bottom_weight, NULL);

  blk_col += blk_col_step;
  uv_blk_col += uv_blk_col_step;

  neighbors_first = LUMA_MIDDLE_COLUMN_NEIGHBORS;
  for (; blk_col < mid_width;
       blk_col += blk_col_step, uv_blk_col += uv_blk_col_step) {
    vp9_apply_temporal_filter_luma_16_avx2(
        y_pre + blk_col, y_pre_stride, block_height, ss_x, ss_y, strength,
        use_whole_blk, y_accum + blk_col, y_count + blk_col, y_dist + blk_col,
        u_dist + uv_blk_col, v_dist + uv_blk_col, neighbors_first,
        neighbors_second, top_weight, bottom_weight, NULL);
  }

  if (!use_whole_blk) {
    top_weight = blk_fw[1];
    bottom_weight = blk_fw[3];
  }

  for (; blk_col < last_width;
       blk_col += blk_col_step, uv_blk_col += uv_blk_col_step) {
    vp9_apply_temporal_filter_luma_16_avx2(
        y_pre + blk_col, y_pre_stride, block_height, ss_x, ss_y, strength,
        use_whole_blk, y_accum + blk_col, y_count + blk_col, y_dist + blk_col,
        u_dist + uv_blk_col, v_dist + uv_blk_col, neighbors_first,
        neighbors_second, top_weight, bottom_weight, NULL);
  }

  neighbors_second = LUMA_RIGHT_COLUMN_NEIGHBORS;
  vp9_apply_temporal_filter_luma_16_avx2(
      y_pre + blk_col, y_pre_stride, block_height, ss_x, ss_y, strength,
      use_whole_blk, y_accum + blk_col, y_count + blk_col, y_dist + blk_col,
      u_dist + uv_blk_col, v_dist + uv_blk_col, neighbors_first,
      neighbors_second, top_weight, bottom_weight, NULL);
}

static void vp9_apply_temporal_filter_chroma_8_sse4(
    const uint8_t *u_pre, const uint8_t *v_pre, int uv_pre_stride,
    unsigned int uv_block_height, int ss_x, int ss_y, int strength,
    uint32_t *u_accum, uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count,
    const uint16_t *y_dist, const uint16_t *u_dist, const uint16_t *v_dist,
    const int16_t *const *neighbors, int top_weight, int bottom_weight,
    const int *blk_fw) {
  const int rounding = (1 << strength) >> 1;
  __m128i weight, mul;
  __m128i u_sum_row_1, u_sum_row_2, u_sum_row_3;
  __m128i v_sum_row_1, v_sum_row_2, v_sum_row_3;
  __m128i u_sum_row, v_sum_row;
  unsigned int h;

  if (blk_fw) {
    weight = _mm_setr_epi16(blk_fw[0], blk_fw[0], blk_fw[0], blk_fw[0],
                            blk_fw[1], blk_fw[1], blk_fw[1], blk_fw[1]);
  } else {
    weight = _mm_set1_epi16(top_weight);
  }

  mul = _mm_load_si128((const __m128i *)neighbors[0]);
  u_sum_row_2 = get_sum_8_sse4(u_dist);
  u_sum_row_3 = get_sum_8_sse4(u_dist + DIST_STRIDE);
  u_sum_row = _mm_adds_epu16(u_sum_row_2, u_sum_row_3);

  v_sum_row_2 = get_sum_8_sse4(v_dist);
  v_sum_row_3 = get_sum_8_sse4(v_dist + DIST_STRIDE);
  v_sum_row = _mm_adds_epu16(v_sum_row_2, v_sum_row_3);

  add_luma_dist_to_8_chroma_mod_sse4(y_dist, ss_x, ss_y, &u_sum_row,
                                     &v_sum_row);

  u_sum_row = average_8_sse4(u_sum_row, &mul, strength, rounding, &weight);
  v_sum_row = average_8_sse4(v_sum_row, &mul, strength, rounding, &weight);
  accumulate_and_store_8_sse4(u_sum_row, u_pre, u_count, u_accum);
  accumulate_and_store_8_sse4(v_sum_row, v_pre, v_count, v_accum);

  u_pre += uv_pre_stride;
  u_dist += DIST_STRIDE;
  v_pre += uv_pre_stride;
  v_dist += DIST_STRIDE;
  u_count += uv_pre_stride;
  u_accum += uv_pre_stride;
  v_count += uv_pre_stride;
  v_accum += uv_pre_stride;
  y_dist += DIST_STRIDE * (1 + ss_y);

  mul = _mm_load_si128((const __m128i *)neighbors[1]);

  for (h = 1; h < uv_block_height - 1; ++h) {
    if (h == uv_block_height / 2) {
      if (blk_fw) {
        weight = _mm_setr_epi16(blk_fw[2], blk_fw[2], blk_fw[2], blk_fw[2],
                                blk_fw[3], blk_fw[3], blk_fw[3], blk_fw[3]);
      } else {
        weight = _mm_set1_epi16(bottom_weight);
      }
    }
    u_sum_row_1 = u_sum_row_2;
    u_sum_row_2 = u_sum_row_3;
    v_sum_row_1 = v_sum_row_2;
    v_sum_row_2 = v_sum_row_3;

    u_sum_row = _mm_adds_epu16(u_sum_row_1, u_sum_row_2);
    u_sum_row_3 = get_sum_8_sse4(u_dist + DIST_STRIDE);
    u_sum_row = _mm_adds_epu16(u_sum_row, u_sum_row_3);

    v_sum_row = _mm_adds_epu16(v_sum_row_1, v_sum_row_2);
    v_sum_row_3 = get_sum_8_sse4(v_dist + DIST_STRIDE);
    v_sum_row = _mm_adds_epu16(v_sum_row, v_sum_row_3);

    add_luma_dist_to_8_chroma_mod_sse4(y_dist, ss_x, ss_y, &u_sum_row,
                                       &v_sum_row);

    u_sum_row = average_8_sse4(u_sum_row, &mul, strength, rounding, &weight);
    v_sum_row = average_8_sse4(v_sum_row, &mul, strength, rounding, &weight);
    accumulate_and_store_8_sse4(u_sum_row, u_pre, u_count, u_accum);
    accumulate_and_store_8_sse4(v_sum_row, v_pre, v_count, v_accum);

    u_pre += uv_pre_stride;
    u_dist += DIST_STRIDE;
    v_pre += uv_pre_stride;
    v_dist += DIST_STRIDE;
    u_count += uv_pre_stride;
    u_accum += uv_pre_stride;
    v_count += uv_pre_stride;
    v_accum += uv_pre_stride;
    y_dist += DIST_STRIDE * (1 + ss_y);
  }

  mul = _mm_load_si128((const __m128i *)neighbors[0]);

  u_sum_row = _mm_adds_epu16(u_sum_row_2, u_sum_row_3);
  v_sum_row = _mm_adds_epu16(v_sum_row_2, v_sum_row_3);
  add_luma_dist_to_8_chroma_mod_sse4(y_dist, ss_x, ss_y, &u_sum_row,
                                     &v_sum_row);

  u_sum_row = average_8_sse4(u_sum_row, &mul, strength, rounding, &weight);
  v_sum_row = average_8_sse4(v_sum_row, &mul, strength, rounding, &weight);
  accumulate_and_store_8_sse4(u_sum_row, u_pre, u_count, u_accum);
  accumulate_and_store_8_sse4(v_sum_row, v_pre, v_count, v_accum);
}

static void vp9_apply_temporal_filter_chroma_16_avx2(
    const uint8_t *u_pre, const uint8_t *v_pre, int uv_pre_stride,
    unsigned int uv_block_height, int ss_x, int ss_y, int strength,
    uint32_t *u_accum, uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count,
    const uint16_t *y_dist, const uint16_t *u_dist, const uint16_t *v_dist,
    const int16_t *const *neighbors_first,
    const int16_t *const *neighbors_second, int top_weight, int bottom_weight,
    const int *blk_fw) {
  const int rounding = (1 << strength) >> 1;
  __m256i weight, mul;
  __m256i u_sum_row_1, u_sum_row_2, u_sum_row_3;
  __m256i v_sum_row_1, v_sum_row_2, v_sum_row_3;
  __m256i u_sum_row, v_sum_row;
  unsigned int h;

  if (blk_fw) {
    weight =
        combine_128_pair(_mm_set1_epi16(blk_fw[0]), _mm_set1_epi16(blk_fw[1]));
  } else {
    weight = _mm256_set1_epi16(top_weight);
  }

  mul = combine_128_pair(_mm_load_si128((const __m128i *)neighbors_first[0]),
                         _mm_load_si128((const __m128i *)neighbors_second[0]));

  u_sum_row_2 = get_sum_16_avx2(u_dist);
  u_sum_row_3 = get_sum_16_avx2(u_dist + DIST_STRIDE);
  u_sum_row = _mm256_adds_epu16(u_sum_row_2, u_sum_row_3);

  v_sum_row_2 = get_sum_16_avx2(v_dist);
  v_sum_row_3 = get_sum_16_avx2(v_dist + DIST_STRIDE);
  v_sum_row = _mm256_adds_epu16(v_sum_row_2, v_sum_row_3);

  add_luma_dist_to_16_chroma_mod_avx2(y_dist, ss_x, ss_y, &u_sum_row,
                                      &v_sum_row);

  u_sum_row = average_16_avx2(u_sum_row, &mul, strength, rounding, &weight);
  v_sum_row = average_16_avx2(v_sum_row, &mul, strength, rounding, &weight);
  accumulate_and_store_16_avx2(u_sum_row, u_pre, u_count, u_accum);
  accumulate_and_store_16_avx2(v_sum_row, v_pre, v_count, v_accum);

  u_pre += uv_pre_stride;
  u_dist += DIST_STRIDE;
  v_pre += uv_pre_stride;
  v_dist += DIST_STRIDE;
  u_count += uv_pre_stride;
  u_accum += uv_pre_stride;
  v_count += uv_pre_stride;
  v_accum += uv_pre_stride;
  y_dist += DIST_STRIDE * (1 + ss_y);

  mul = combine_128_pair(_mm_load_si128((const __m128i *)neighbors_first[1]),
                         _mm_load_si128((const __m128i *)neighbors_second[1]));

  for (h = 1; h < uv_block_height - 1; ++h) {
    if (h == uv_block_height / 2) {
      if (blk_fw) {
        weight = combine_128_pair(_mm_set1_epi16(blk_fw[2]),
                                  _mm_set1_epi16(blk_fw[3]));
      } else {
        weight = _mm256_set1_epi16(bottom_weight);
      }
    }
    u_sum_row_1 = u_sum_row_2;
    u_sum_row_2 = u_sum_row_3;
    v_sum_row_1 = v_sum_row_2;
    v_sum_row_2 = v_sum_row_3;

    u_sum_row = _mm256_adds_epu16(u_sum_row_1, u_sum_row_2);
    u_sum_row_3 = get_sum_16_avx2(u_dist + DIST_STRIDE);
    u_sum_row = _mm256_adds_epu16(u_sum_row, u_sum_row_3);

    v_sum_row = _mm256_adds_epu16(v_sum_row_1, v_sum_row_2);
    v_sum_row_3 = get_sum_16_avx2(v_dist + DIST_STRIDE);
    v_sum_row = _mm256_adds_epu16(v_sum_row, v_sum_row_3);

    add_luma_dist_to_16_chroma_mod_avx2(y_dist, ss_x, ss_y, &u_sum_row,
                                        &v_sum_row);

    u_sum_row = average_16_avx2(u_sum_row, &mul, strength, rounding, &weight);
    v_sum_row = average_16_avx2(v_sum_row, &mul, strength, rounding, &weight);
    accumulate_and_store_16_avx2(u_sum_row, u_pre, u_count, u_accum);
    accumulate_and_store_16_avx2(v_sum_row, v_pre, v_count, v_accum);

    u_pre += uv_pre_stride;
    u_dist += DIST_STRIDE;
    v_pre += uv_pre_stride;
    v_dist += DIST_STRIDE;
    u_count += uv_pre_stride;
    u_accum += uv_pre_stride;
    v_count += uv_pre_stride;
    v_accum += uv_pre_stride;
    y_dist += DIST_STRIDE * (1 + ss_y);
  }

  mul = combine_128_pair(_mm_load_si128((const __m128i *)neighbors_first[0]),
                         _mm_load_si128((const __m128i *)neighbors_second[0]));

  u_sum_row = _mm256_adds_epu16(u_sum_row_2, u_sum_row_3);
  v_sum_row = _mm256_adds_epu16(v_sum_row_2, v_sum_row_3);
  add_luma_dist_to_16_chroma_mod_avx2(y_dist, ss_x, ss_y, &u_sum_row,
                                      &v_sum_row);

  u_sum_row = average_16_avx2(u_sum_row, &mul, strength, rounding, &weight);
  v_sum_row = average_16_avx2(v_sum_row, &mul, strength, rounding, &weight);
  accumulate_and_store_16_avx2(u_sum_row, u_pre, u_count, u_accum);
  accumulate_and_store_16_avx2(v_sum_row, v_pre, v_count, v_accum);
}

static void vp9_apply_temporal_filter_chroma_avx2(
    const uint8_t *u_pre, const uint8_t *v_pre, int uv_pre_stride,
    unsigned int block_width, unsigned int block_height, int ss_x, int ss_y,
    int strength, const int *blk_fw, int use_whole_blk, uint32_t *u_accum,
    uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count,
    const uint16_t *y_dist, const uint16_t *u_dist, const uint16_t *v_dist) {
  const unsigned int uv_width = block_width >> ss_x;
  const unsigned int uv_height = block_height >> ss_y;
  int top_weight = blk_fw[0];
  int bottom_weight = use_whole_blk ? blk_fw[0] : blk_fw[2];
  const int16_t *const *left_n;
  const int16_t *const *mid_n;
  const int16_t *const *right_n;

  if (uv_width == 8) {
    const int16_t *const *neighbors =
        ss_y ? CHROMA_DOUBLE_SS_SINGLE_COLUMN_NEIGHBORS
             : CHROMA_SINGLE_SS_SINGLE_COLUMN_NEIGHBORS;
    assert(ss_x);
    vp9_apply_temporal_filter_chroma_8_sse4(
        u_pre, v_pre, uv_pre_stride, uv_height, ss_x, ss_y, strength, u_accum,
        u_count, v_accum, v_count, y_dist, u_dist, v_dist, neighbors,
        top_weight, bottom_weight, use_whole_blk ? NULL : blk_fw);
    return;
  }

  if (ss_x && ss_y) {
    left_n = CHROMA_DOUBLE_SS_LEFT_COLUMN_NEIGHBORS;
    mid_n = CHROMA_DOUBLE_SS_MIDDLE_COLUMN_NEIGHBORS;
    right_n = CHROMA_DOUBLE_SS_RIGHT_COLUMN_NEIGHBORS;
  } else if (ss_x || ss_y) {
    left_n = CHROMA_SINGLE_SS_LEFT_COLUMN_NEIGHBORS;
    mid_n = CHROMA_SINGLE_SS_MIDDLE_COLUMN_NEIGHBORS;
    right_n = CHROMA_SINGLE_SS_RIGHT_COLUMN_NEIGHBORS;
  } else {
    left_n = CHROMA_NO_SS_LEFT_COLUMN_NEIGHBORS;
    mid_n = CHROMA_NO_SS_MIDDLE_COLUMN_NEIGHBORS;
    right_n = CHROMA_NO_SS_RIGHT_COLUMN_NEIGHBORS;
  }

  if (uv_width == 16) {
    vp9_apply_temporal_filter_chroma_16_avx2(
        u_pre, v_pre, uv_pre_stride, uv_height, ss_x, ss_y, strength, u_accum,
        u_count, v_accum, v_count, y_dist, u_dist, v_dist, left_n, right_n,
        top_weight, bottom_weight, use_whole_blk ? NULL : blk_fw);
    return;
  }

  // uv_width >= 32
  {
    unsigned int blk_col = 0, uv_blk_col = 0;
    const unsigned int uv_blk_col_step = 16, blk_col_step = 16 << ss_x;
    const unsigned int uv_mid_width = uv_width >> 1,
                       uv_last_width = uv_width - uv_blk_col_step;

    if (uv_width == 32) {
      vp9_apply_temporal_filter_chroma_16_avx2(
          u_pre, v_pre, uv_pre_stride, uv_height, ss_x, ss_y, strength, u_accum,
          u_count, v_accum, v_count, y_dist, u_dist, v_dist, left_n, mid_n,
          top_weight, bottom_weight, NULL);
      if (!use_whole_blk) {
        top_weight = blk_fw[1];
        bottom_weight = blk_fw[3];
      }
      vp9_apply_temporal_filter_chroma_16_avx2(
          u_pre + 16, v_pre + 16, uv_pre_stride, uv_height, ss_x, ss_y,
          strength, u_accum + 16, u_count + 16, v_accum + 16, v_count + 16,
          y_dist + (16 << ss_x), u_dist + 16, v_dist + 16, mid_n, right_n,
          top_weight, bottom_weight, NULL);
      return;
    }

    vp9_apply_temporal_filter_chroma_16_avx2(
        u_pre, v_pre, uv_pre_stride, uv_height, ss_x, ss_y, strength, u_accum,
        u_count, v_accum, v_count, y_dist, u_dist, v_dist, left_n, mid_n,
        top_weight, bottom_weight, NULL);
    blk_col += blk_col_step;
    uv_blk_col += uv_blk_col_step;

    for (; uv_blk_col < uv_mid_width;
         blk_col += blk_col_step, uv_blk_col += uv_blk_col_step) {
      vp9_apply_temporal_filter_chroma_16_avx2(
          u_pre + uv_blk_col, v_pre + uv_blk_col, uv_pre_stride, uv_height,
          ss_x, ss_y, strength, u_accum + uv_blk_col, u_count + uv_blk_col,
          v_accum + uv_blk_col, v_count + uv_blk_col, y_dist + blk_col,
          u_dist + uv_blk_col, v_dist + uv_blk_col, mid_n, mid_n, top_weight,
          bottom_weight, NULL);
    }

    if (!use_whole_blk) {
      top_weight = blk_fw[1];
      bottom_weight = blk_fw[3];
    }

    for (; uv_blk_col < uv_last_width;
         blk_col += blk_col_step, uv_blk_col += uv_blk_col_step) {
      vp9_apply_temporal_filter_chroma_16_avx2(
          u_pre + uv_blk_col, v_pre + uv_blk_col, uv_pre_stride, uv_height,
          ss_x, ss_y, strength, u_accum + uv_blk_col, u_count + uv_blk_col,
          v_accum + uv_blk_col, v_count + uv_blk_col, y_dist + blk_col,
          u_dist + uv_blk_col, v_dist + uv_blk_col, mid_n, mid_n, top_weight,
          bottom_weight, NULL);
    }

    vp9_apply_temporal_filter_chroma_16_avx2(
        u_pre + uv_blk_col, v_pre + uv_blk_col, uv_pre_stride, uv_height, ss_x,
        ss_y, strength, u_accum + uv_blk_col, u_count + uv_blk_col,
        v_accum + uv_blk_col, v_count + uv_blk_col, y_dist + blk_col,
        u_dist + uv_blk_col, v_dist + uv_blk_col, mid_n, right_n, top_weight,
        bottom_weight, NULL);
  }
}

void vp9_apply_temporal_filter_avx2(
    const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre,
    int y_pre_stride, const uint8_t *u_src, const uint8_t *v_src,
    int uv_src_stride, const uint8_t *u_pre, const uint8_t *v_pre,
    int uv_pre_stride, unsigned int block_width, unsigned int block_height,
    int ss_x, int ss_y, int strength, const int *const blk_fw,
    int use_whole_blk, uint32_t *y_accum, uint16_t *y_count, uint32_t *u_accum,
    uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count) {
  const unsigned int chroma_height = block_height >> ss_y;
  const unsigned int chroma_width = block_width >> ss_x;

  DECLARE_ALIGNED(32, uint16_t, y_dist[BH * DIST_STRIDE]) = { 0 };
  DECLARE_ALIGNED(32, uint16_t, u_dist[BH * DIST_STRIDE]) = { 0 };
  DECLARE_ALIGNED(32, uint16_t, v_dist[BH * DIST_STRIDE]) = { 0 };

  uint16_t *y_dist_ptr = y_dist + 1;
  uint16_t *u_dist_ptr = u_dist + 1;
  uint16_t *v_dist_ptr = v_dist + 1;
  const uint8_t *y_src_ptr = y_src, *u_src_ptr = u_src, *v_src_ptr = v_src;
  const uint8_t *y_pre_ptr = y_pre, *u_pre_ptr = u_pre, *v_pre_ptr = v_pre;
  unsigned int row, blk_col;

  assert(block_width <= BW && "block width too large");
  assert(block_height <= BH && "block height too large");
  assert(block_width % 16 == 0 && "block width must be multiple of 16");
  assert(block_height % 2 == 0 && "block height must be even");
  assert((ss_x == 0 || ss_x == 1) && (ss_y == 0 || ss_y == 1) &&
         "invalid chroma subsampling");
  assert(strength >= 0 && strength <= 6 && "invalid temporal filter strength");
  assert(blk_fw[0] >= 0 && "filter weight must be positive");
  assert(
      (use_whole_blk || (blk_fw[1] >= 0 && blk_fw[2] >= 0 && blk_fw[3] >= 0)) &&
      "subblock filter weight must be positive");
  assert(blk_fw[0] <= 2 && "subblock filter weight must be less than 2");
  assert(
      (use_whole_blk || (blk_fw[1] <= 2 && blk_fw[2] <= 2 && blk_fw[3] <= 2)) &&
      "subblock filter weight must be less than 2");

  if (block_width % 32 == 0) {
    for (row = 0; row < block_height; row++) {
      for (blk_col = 0; blk_col < block_width; blk_col += 32) {
        store_dist_32_avx2(y_src_ptr + blk_col, y_pre_ptr + blk_col,
                           y_dist_ptr + blk_col);
      }
      y_src_ptr += y_src_stride;
      y_pre_ptr += y_pre_stride;
      y_dist_ptr += DIST_STRIDE;
    }
  } else {
    for (row = 0; row < block_height; row++) {
      for (blk_col = 0; blk_col < block_width; blk_col += 16) {
        store_dist_16_avx2(y_src_ptr + blk_col, y_pre_ptr + blk_col,
                           y_dist_ptr + blk_col);
      }
      y_src_ptr += y_src_stride;
      y_pre_ptr += y_pre_stride;
      y_dist_ptr += DIST_STRIDE;
    }
  }

  if (chroma_width % 16 == 0) {
    for (row = 0; row < chroma_height; row++) {
      for (blk_col = 0; blk_col < chroma_width; blk_col += 16) {
        store_dist_16_avx2(u_src_ptr + blk_col, u_pre_ptr + blk_col,
                           u_dist_ptr + blk_col);
        store_dist_16_avx2(v_src_ptr + blk_col, v_pre_ptr + blk_col,
                           v_dist_ptr + blk_col);
      }
      u_src_ptr += uv_src_stride;
      u_pre_ptr += uv_pre_stride;
      u_dist_ptr += DIST_STRIDE;
      v_src_ptr += uv_src_stride;
      v_pre_ptr += uv_pre_stride;
      v_dist_ptr += DIST_STRIDE;
    }
  } else {
    for (row = 0; row < chroma_height; row++) {
      for (blk_col = 0; blk_col < chroma_width; blk_col += 8) {
        store_dist_8_avx2(u_src_ptr + blk_col, u_pre_ptr + blk_col,
                          u_dist_ptr + blk_col);
        store_dist_8_avx2(v_src_ptr + blk_col, v_pre_ptr + blk_col,
                          v_dist_ptr + blk_col);
      }
      u_src_ptr += uv_src_stride;
      u_pre_ptr += uv_pre_stride;
      u_dist_ptr += DIST_STRIDE;
      v_src_ptr += uv_src_stride;
      v_pre_ptr += uv_pre_stride;
      v_dist_ptr += DIST_STRIDE;
    }
  }

  y_dist_ptr = y_dist + 1;
  u_dist_ptr = u_dist + 1;
  v_dist_ptr = v_dist + 1;

  vp9_apply_temporal_filter_luma_avx2(y_pre, y_pre_stride, block_width,
                                      block_height, ss_x, ss_y, strength,
                                      blk_fw, use_whole_blk, y_accum, y_count,
                                      y_dist_ptr, u_dist_ptr, v_dist_ptr);

  vp9_apply_temporal_filter_chroma_avx2(
      u_pre, v_pre, uv_pre_stride, block_width, block_height, ss_x, ss_y,
      strength, blk_fw, use_whole_blk, u_accum, u_count, v_accum, v_count,
      y_dist_ptr, u_dist_ptr, v_dist_ptr);
}
