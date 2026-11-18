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
#include "./vpx_dsp_rtcd.h"
#include "vp9/encoder/vp9_temporal_filter.h"
#include "vp9/encoder/vp9_temporal_filter_constants.h"

static INLINE void highbd_shuffle_12tap_filter_avx2(const int16_t *filter,
                                                    __m256i *f) {
  const __m256i f_low =
      _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)filter));
  const __m256i f_high = _mm256_broadcastsi128_si256(
      _mm_loadl_epi64((const __m128i *)(filter + 8)));

  f[0] = _mm256_shuffle_epi32(f_low, 0x00);
  f[1] = _mm256_shuffle_epi32(f_low, 0x55);
  f[2] = _mm256_shuffle_epi32(f_low, 0xaa);
  f[3] = _mm256_shuffle_epi32(f_low, 0xff);
  f[4] = _mm256_shuffle_epi32(f_high, 0x00);
  f[5] = _mm256_shuffle_epi32(f_high, 0x55);
}

static INLINE __m256i highbd_convolve_12tap(const __m256i *s,
                                            const __m256i *f) {
  const __m256i res_0 = _mm256_madd_epi16(s[0], f[0]);
  const __m256i res_1 = _mm256_madd_epi16(s[1], f[1]);
  const __m256i res_2 = _mm256_madd_epi16(s[2], f[2]);
  const __m256i res_3 = _mm256_madd_epi16(s[3], f[3]);
  const __m256i res_4 = _mm256_madd_epi16(s[4], f[4]);
  const __m256i res_5 = _mm256_madd_epi16(s[5], f[5]);

  const __m256i res =
      _mm256_add_epi32(_mm256_add_epi32(res_0, res_1),
                       _mm256_add_epi32(_mm256_add_epi32(res_2, res_3),
                                        _mm256_add_epi32(res_4, res_5)));
  return res;
}

static INLINE void reuse_src_data_avx2(const __m256i *src, __m256i *des) {
  des[0] = src[0];
  des[1] = src[1];
  des[2] = src[2];
  des[3] = src[3];
  des[4] = src[4];
}

void vpx_highbd_convolve12_horiz_avx2(const uint16_t *src, ptrdiff_t src_stride,
                                      uint16_t *dst, ptrdiff_t dst_stride,
                                      const InterpKernel12 *filter, int x0_q4,
                                      int x_step_q4, int y0_q4, int y_step_q4,
                                      int w, int h, int bd) {
  assert(x_step_q4 == 16);
  (void)y0_q4;
  (void)x_step_q4;
  (void)y_step_q4;
  const uint16_t *src_ptr = src;
  src_ptr -= MAX_FILTER_TAP / 2 - 1;
  __m256i s[6], f[6];
  const __m256i rounding = _mm256_set1_epi32(1 << (FILTER_BITS - 1));
  const __m256i max = _mm256_set1_epi16((1 << bd) - 1);
  highbd_shuffle_12tap_filter_avx2(filter[x0_q4], f);

  for (int j = 0; j < w; j += 8) {
    for (int i = 0; i < h; i += 2) {
      // s00 s01 s02 s03 s04 s05 s06 s07 s08 s09 s010 s011 s012 s013 s014 s015
      const __m256i row0 =
          _mm256_loadu_si256((const __m256i *)&src_ptr[i * src_stride + j]);
      // s10 s11 s12 s13 s14 s15 s16 s17 s18 s19 s110 s111 s112 s113 s114
      // s115
      const __m256i row1 = _mm256_loadu_si256(
          (const __m256i *)&src_ptr[(i + 1) * src_stride + j]);
      // s016 s017 s018 s019 s020 s021 s022 s023
      const __m128i row0_16 =
          _mm_loadu_si128((const __m128i *)&src_ptr[i * src_stride + j + 16]);
      // s116 s117 s118 s119 s120 s121 s122 s123
      const __m128i row1_16 = _mm_loadu_si128(
          (const __m128i *)&src_ptr[(i + 1) * src_stride + j + 16]);

      // s00 s01 s02 s03 s04 s05 s06 s07 | s10 s11 s12 s13 s14 s15 s16 s17
      const __m256i r0 = _mm256_permute2x128_si256(row0, row1, 0x20);
      // s08 s09 s010 s011 s012 s013 s014 s015 | s18 s19 s110 s111 s112 s113
      // s114 s115
      const __m256i r1 = _mm256_permute2x128_si256(row0, row1, 0x31);
      // s016 s017 s018 s019 s020 s021 s022 s023 | s116 s117 s118 s119 s120 s121
      // s122 s123
      const __m256i r2 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(row0_16), row1_16, 1);

      // even pixels
      s[0] = r0;
      s[1] = _mm256_alignr_epi8(r1, r0, 4);
      s[2] = _mm256_alignr_epi8(r1, r0, 8);
      s[3] = _mm256_alignr_epi8(r1, r0, 12);
      s[4] = r1;
      s[5] = _mm256_alignr_epi8(r2, r1, 4);

      // 00 02 04 06 | 10 12 14 16
      __m256i res_even = highbd_convolve_12tap(s, f);
      res_even =
          _mm256_srai_epi32(_mm256_add_epi32(res_even, rounding), FILTER_BITS);

      // odd pixels
      s[0] = _mm256_alignr_epi8(r1, r0, 2);
      s[1] = _mm256_alignr_epi8(r1, r0, 6);
      s[2] = _mm256_alignr_epi8(r1, r0, 10);
      s[3] = _mm256_alignr_epi8(r1, r0, 14);
      s[4] = _mm256_alignr_epi8(r2, r1, 2);
      s[5] = _mm256_alignr_epi8(r2, r1, 6);

      // 01 03 05 07 | 11 13 15 17
      __m256i res_odd = highbd_convolve_12tap(s, f);
      res_odd =
          _mm256_srai_epi32(_mm256_add_epi32(res_odd, rounding), FILTER_BITS);

      // 00 01 02 03 | 10 11 12 13
      const __m256i res_0 = _mm256_unpacklo_epi32(res_even, res_odd);
      // 04 05 06 07 | 14 15 16 17
      const __m256i res_1 = _mm256_unpackhi_epi32(res_even, res_odd);
      // 00 01 02 03 | 04 05 06 07 | 10 11 12 13 | 14 15 16 17
      const __m256i res_2 = _mm256_packus_epi32(res_0, res_1);
      const __m256i res = _mm256_min_epi16(res_2, max);
      _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j],
                       _mm256_castsi256_si128(res));
      if (i + 1 < h) {
        _mm_storeu_si128((__m128i *)(&dst[(i + 1) * dst_stride + j]),
                         _mm256_extractf128_si256(res, 1));
      }
    }
  }
}

void vpx_highbd_convolve12_vert_avx2(const uint16_t *src, ptrdiff_t src_stride,
                                     uint16_t *dst, ptrdiff_t dst_stride,
                                     const InterpKernel12 *filter, int x0_q4,
                                     int x_step_q4, int y0_q4, int y_step_q4,
                                     int w, int h, int bd) {
  assert(y_step_q4 == 16);
  (void)x0_q4;
  (void)x_step_q4;
  (void)y_step_q4;
  const uint16_t *src_ptr = src;
  src_ptr -= src_stride * (MAX_FILTER_TAP / 2 - 1);
  __m256i s[12], f[6];
  const __m256i rounding = _mm256_set1_epi32(((1 << FILTER_BITS) >> 1));
  const __m256i max = _mm256_set1_epi16((1 << bd) - 1);
  highbd_shuffle_12tap_filter_avx2(filter[y0_q4], f);

  for (int j = 0; j < w; j += 8) {
    __m128i s0 =
        _mm_loadu_si128((const __m128i *)(src_ptr + 0 * src_stride + j));
    __m128i s1 =
        _mm_loadu_si128((const __m128i *)(src_ptr + 1 * src_stride + j));
    __m128i s2 =
        _mm_loadu_si128((const __m128i *)(src_ptr + 2 * src_stride + j));
    __m128i s3 =
        _mm_loadu_si128((const __m128i *)(src_ptr + 3 * src_stride + j));
    __m128i s4 =
        _mm_loadu_si128((const __m128i *)(src_ptr + 4 * src_stride + j));
    __m128i s5 =
        _mm_loadu_si128((const __m128i *)(src_ptr + 5 * src_stride + j));
    __m128i s6 =
        _mm_loadu_si128((const __m128i *)(src_ptr + 6 * src_stride + j));
    __m128i s7 =
        _mm_loadu_si128((const __m128i *)(src_ptr + 7 * src_stride + j));
    __m128i s8 =
        _mm_loadu_si128((const __m128i *)(src_ptr + 8 * src_stride + j));
    __m128i s9 =
        _mm_loadu_si128((const __m128i *)(src_ptr + 9 * src_stride + j));
    __m128i s10t =
        _mm_loadu_si128((const __m128i *)(src_ptr + 10 * src_stride + j));

    __m256i r01 = _mm256_inserti128_si256(_mm256_castsi128_si256(s0), s1, 1);
    __m256i r12 = _mm256_inserti128_si256(_mm256_castsi128_si256(s1), s2, 1);
    __m256i r23 = _mm256_inserti128_si256(_mm256_castsi128_si256(s2), s3, 1);
    __m256i r34 = _mm256_inserti128_si256(_mm256_castsi128_si256(s3), s4, 1);
    __m256i r45 = _mm256_inserti128_si256(_mm256_castsi128_si256(s4), s5, 1);
    __m256i r56 = _mm256_inserti128_si256(_mm256_castsi128_si256(s5), s6, 1);
    __m256i r67 = _mm256_inserti128_si256(_mm256_castsi128_si256(s6), s7, 1);
    __m256i r78 = _mm256_inserti128_si256(_mm256_castsi128_si256(s7), s8, 1);
    __m256i r89 = _mm256_inserti128_si256(_mm256_castsi128_si256(s8), s9, 1);
    __m256i r910 = _mm256_inserti128_si256(_mm256_castsi128_si256(s9), s10t, 1);

    s[0] = _mm256_unpacklo_epi16(r01, r12);
    s[1] = _mm256_unpacklo_epi16(r23, r34);
    s[2] = _mm256_unpacklo_epi16(r45, r56);
    s[3] = _mm256_unpacklo_epi16(r67, r78);
    s[4] = _mm256_unpacklo_epi16(r89, r910);

    s[6] = _mm256_unpackhi_epi16(r01, r12);
    s[7] = _mm256_unpackhi_epi16(r23, r34);
    s[8] = _mm256_unpackhi_epi16(r45, r56);
    s[9] = _mm256_unpackhi_epi16(r67, r78);
    s[10] = _mm256_unpackhi_epi16(r89, r910);
    for (int i = 0; i < h; i += 2) {
      const __m128i s10 = _mm_loadu_si128(
          (const __m128i *)(src_ptr + (i + 10) * src_stride + j));
      const __m128i s11 = _mm_loadu_si128(
          (const __m128i *)(src_ptr + (i + 11) * src_stride + j));
      const __m128i s12 = _mm_loadu_si128(
          (const __m128i *)(src_ptr + (i + 12) * src_stride + j));
      __m256i r1011 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s10), s11, 1);
      __m256i r1112 =
          _mm256_inserti128_si256(_mm256_castsi128_si256(s11), s12, 1);

      s[5] = _mm256_unpacklo_epi16(r1011, r1112);
      s[11] = _mm256_unpackhi_epi16(r1011, r1112);

      // 00 01 02 03 | 10 11 12 13
      const __m256i res_a = highbd_convolve_12tap(s, f);
      __m256i res_a_round =
          _mm256_srai_epi32(_mm256_add_epi32(res_a, rounding), FILTER_BITS);
      // 04 05 06 07 | 14 15 16 17
      const __m256i res_b = highbd_convolve_12tap(s + 6, f);
      __m256i res_b_round =
          _mm256_srai_epi32(_mm256_add_epi32(res_b, rounding), FILTER_BITS);

      // 00 01 02 03 | 04 05 06 07 | 10 11 12 13 | 14 15 16 17
      const __m256i res_0 = _mm256_packus_epi32(res_a_round, res_b_round);
      const __m256i res = _mm256_min_epi16(res_0, max);
      _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j],
                       _mm256_castsi256_si128(res));

      _mm_storeu_si128((__m128i *)(&dst[(i + 1) * dst_stride + j]),
                       _mm256_extractf128_si256(res, 1));

      reuse_src_data_avx2(s + 1, s);
      reuse_src_data_avx2(s + 7, s + 6);
    }
  }
}

void vpx_highbd_convolve12_avx2(const uint16_t *src, ptrdiff_t src_stride,
                                uint16_t *dst, ptrdiff_t dst_stride,
                                const InterpKernel12 *filter, int x0_q4,
                                int x_step_q4, int y0_q4, int y_step_q4, int w,
                                int h, int bd) {
  assert(x_step_q4 == 16 && y_step_q4 == 16);
  assert(h == 32 || h == 16 || h == 8);
  assert(w == 32 || w == 16 || w == 8);
  DECLARE_ALIGNED(32, uint16_t, temp[BW * (BH + MAX_FILTER_TAP - 1)]);
  const int temp_stride = BW;
  const int intermediate_height =
      (((h - 1) * y_step_q4 + y0_q4) >> SUBPEL_BITS) + MAX_FILTER_TAP;

  vpx_highbd_convolve12_horiz_avx2(src - src_stride * (MAX_FILTER_TAP / 2 - 1),
                                   src_stride, temp, temp_stride, filter, x0_q4,
                                   x_step_q4, y0_q4, y_step_q4, w,
                                   intermediate_height, bd);
  vpx_highbd_convolve12_vert_avx2(temp + temp_stride * (MAX_FILTER_TAP / 2 - 1),
                                  temp_stride, dst, dst_stride, filter, x0_q4,
                                  x_step_q4, y0_q4, y_step_q4, w, h, bd);
}

static INLINE __m256i combine_128_pair_hbd(const __m128i lo, const __m128i hi) {
  return _mm256_inserti128_si256(_mm256_castsi128_si256(lo), hi, 1);
}

static INLINE void highbd_store_dist_8_avx2(const uint16_t *a,
                                            const uint16_t *b, uint32_t *dst) {
  const __m128i a_reg = _mm_loadu_si128((const __m128i *)a);
  const __m128i b_reg = _mm_loadu_si128((const __m128i *)b);
  const __m256i diff = _mm256_sub_epi32(_mm256_cvtepu16_epi32(a_reg),
                                        _mm256_cvtepu16_epi32(b_reg));
  _mm256_storeu_si256((__m256i *)dst, _mm256_mullo_epi32(diff, diff));
}

static INLINE __m256i highbd_get_sum_8_avx2(const uint32_t *dist) {
  const __m256i dist_reg = _mm256_loadu_si256((const __m256i *)dist);
  const __m256i dist_left = _mm256_loadu_si256((const __m256i *)(dist - 1));
  const __m256i dist_right = _mm256_loadu_si256((const __m256i *)(dist + 1));
  return _mm256_add_epi32(_mm256_add_epi32(dist_reg, dist_left), dist_right);
}

static INLINE __m256i highbd_average_8_avx2(__m256i sum,
                                            const __m256i *mul_constants,
                                            int strength, int rounding,
                                            const __m256i *weight) {
  const __m128i strength_u128 = _mm_set_epi32(0, 0, 0, strength);
  const __m256i rounding_u32 = _mm256_set1_epi32(rounding);
  const __m256i sixteen = _mm256_set1_epi32(16);

  const __m256i mul_even =
      _mm256_srli_epi64(_mm256_mul_epu32(sum, *mul_constants), 32);
  const __m256i mul_odd = _mm256_mul_epu32(
      _mm256_srli_epi64(sum, 32), _mm256_srli_epi64(*mul_constants, 32));
  const __m256i mul = _mm256_blend_epi32(mul_even, mul_odd, 0xAA);

  __m256i out = _mm256_add_epi32(mul, rounding_u32);
  out = _mm256_srl_epi32(out, strength_u128);
  out = _mm256_min_epu32(out, sixteen);
  out = _mm256_sub_epi32(sixteen, out);
  return _mm256_mullo_epi32(out, *weight);
}

static INLINE void highbd_accumulate_and_store_8_avx2(const __m256i sum_u32,
                                                      const uint16_t *pred,
                                                      uint16_t *count,
                                                      uint32_t *accumulator) {
  const __m128i sum_u16 = _mm256_castsi256_si128(
      _mm256_permute4x64_epi64(_mm256_packus_epi32(sum_u32, sum_u32), 0xD8));
  const __m128i pred_u16 = _mm_loadu_si128((const __m128i *)pred);
  const __m128i count_u16 = _mm_loadu_si128((const __m128i *)count);
  const __m256i pred_u32 = _mm256_cvtepu16_epi32(pred_u16);
  const __m256i accum_cur = _mm256_loadu_si256((const __m256i *)accumulator);

  _mm_storeu_si128((__m128i *)count, _mm_adds_epu16(count_u16, sum_u16));
  _mm256_storeu_si256(
      (__m256i *)accumulator,
      _mm256_add_epi32(accum_cur, _mm256_mullo_epi32(sum_u32, pred_u32)));
}

static INLINE void highbd_read_chroma_dist_row_8_avx2(int ss_x,
                                                      const uint32_t *u_dist,
                                                      const uint32_t *v_dist,
                                                      __m256i *u_out,
                                                      __m256i *v_out) {
  if (!ss_x) {
    *u_out = _mm256_loadu_si256((const __m256i *)u_dist);
    *v_out = _mm256_loadu_si256((const __m256i *)v_dist);
  } else {
    const __m128i u_reg = _mm_loadu_si128((const __m128i *)u_dist);
    const __m128i v_reg = _mm_loadu_si128((const __m128i *)v_dist);
    const __m256i dup_idx = _mm256_setr_epi32(0, 0, 1, 1, 2, 2, 3, 3);
    *u_out =
        _mm256_permutevar8x32_epi32(_mm256_castsi128_si256(u_reg), dup_idx);
    *v_out =
        _mm256_permutevar8x32_epi32(_mm256_castsi128_si256(v_reg), dup_idx);
  }
}

static INLINE void highbd_add_luma_dist_to_8_chroma_mod_avx2(
    const uint32_t *y_dist, int ss_x, int ss_y, __m256i *u_mod,
    __m256i *v_mod) {
  __m256i y_reg;
  if (!ss_x) {
    y_reg = _mm256_loadu_si256((const __m256i *)y_dist);
    if (ss_y == 1) {
      y_reg = _mm256_add_epi32(
          y_reg, _mm256_loadu_si256((const __m256i *)(y_dist + DIST_STRIDE)));
    }
  } else {
    __m256i y_fst = _mm256_loadu_si256((const __m256i *)y_dist);
    __m256i y_snd = _mm256_loadu_si256((const __m256i *)(y_dist + 8));
    if (ss_y == 1) {
      y_fst = _mm256_add_epi32(
          y_fst, _mm256_loadu_si256((const __m256i *)(y_dist + DIST_STRIDE)));
      y_snd = _mm256_add_epi32(
          y_snd,
          _mm256_loadu_si256((const __m256i *)(y_dist + 8 + DIST_STRIDE)));
    }
    y_reg = _mm256_permute4x64_epi64(_mm256_hadd_epi32(y_fst, y_snd), 0xD8);
  }

  *u_mod = _mm256_add_epi32(*u_mod, y_reg);
  *v_mod = _mm256_add_epi32(*v_mod, y_reg);
}

static void vp9_highbd_apply_temporal_filter_luma_8_avx2(
    const uint16_t *y_pre, int y_pre_stride, unsigned int block_height,
    int ss_x, int ss_y, int strength, int use_whole_blk, uint32_t *y_accum,
    uint16_t *y_count, const uint32_t *y_dist, const uint32_t *u_dist,
    const uint32_t *v_dist, const uint32_t *const *neighbors_first,
    const uint32_t *const *neighbors_second, int top_weight,
    int bottom_weight) {
  const int rounding = (1 << strength) >> 1;
  __m256i weight = _mm256_set1_epi32(top_weight);
  __m256i mul;
  __m256i sum_row_1, sum_row_2, sum_row_3;
  __m256i u_reg, v_reg, sum_row;
  unsigned int h;

  assert(strength >= 4 && strength <= 14 &&
         "invalid adjusted temporal filter strength");

  mul = combine_128_pair_hbd(
      _mm_load_si128((const __m128i *)neighbors_first[0]),
      _mm_load_si128((const __m128i *)neighbors_second[0]));

  sum_row_2 = highbd_get_sum_8_avx2(y_dist);
  sum_row_3 = highbd_get_sum_8_avx2(y_dist + DIST_STRIDE);
  sum_row = _mm256_add_epi32(sum_row_2, sum_row_3);

  highbd_read_chroma_dist_row_8_avx2(ss_x, u_dist, v_dist, &u_reg, &v_reg);
  sum_row = _mm256_add_epi32(sum_row, u_reg);
  sum_row = _mm256_add_epi32(sum_row, v_reg);

  sum_row = highbd_average_8_avx2(sum_row, &mul, strength, rounding, &weight);
  highbd_accumulate_and_store_8_avx2(sum_row, y_pre, y_count, y_accum);

  y_pre += y_pre_stride;
  y_count += y_pre_stride;
  y_accum += y_pre_stride;
  y_dist += DIST_STRIDE;
  u_dist += DIST_STRIDE;
  v_dist += DIST_STRIDE;

  mul = combine_128_pair_hbd(
      _mm_load_si128((const __m128i *)neighbors_first[1]),
      _mm_load_si128((const __m128i *)neighbors_second[1]));

  for (h = 1; h < block_height - 1; ++h) {
    if (!use_whole_blk && h == block_height / 2) {
      weight = _mm256_set1_epi32(bottom_weight);
    }
    sum_row_1 = sum_row_2;
    sum_row_2 = sum_row_3;
    sum_row = _mm256_add_epi32(sum_row_1, sum_row_2);
    sum_row_3 = highbd_get_sum_8_avx2(y_dist + DIST_STRIDE);
    sum_row = _mm256_add_epi32(sum_row, sum_row_3);

    if (ss_y == 0 || h % 2 == 0) {
      highbd_read_chroma_dist_row_8_avx2(ss_x, u_dist, v_dist, &u_reg, &v_reg);
      u_dist += DIST_STRIDE;
      v_dist += DIST_STRIDE;
    }

    sum_row = _mm256_add_epi32(sum_row, u_reg);
    sum_row = _mm256_add_epi32(sum_row, v_reg);

    sum_row = highbd_average_8_avx2(sum_row, &mul, strength, rounding, &weight);
    highbd_accumulate_and_store_8_avx2(sum_row, y_pre, y_count, y_accum);

    y_pre += y_pre_stride;
    y_count += y_pre_stride;
    y_accum += y_pre_stride;
    y_dist += DIST_STRIDE;
  }

  mul = combine_128_pair_hbd(
      _mm_load_si128((const __m128i *)neighbors_first[0]),
      _mm_load_si128((const __m128i *)neighbors_second[0]));

  sum_row_1 = sum_row_2;
  sum_row_2 = sum_row_3;
  sum_row = _mm256_add_epi32(sum_row_1, sum_row_2);

  if (ss_y == 0) {
    highbd_read_chroma_dist_row_8_avx2(ss_x, u_dist, v_dist, &u_reg, &v_reg);
  }

  sum_row = _mm256_add_epi32(sum_row, u_reg);
  sum_row = _mm256_add_epi32(sum_row, v_reg);

  sum_row = highbd_average_8_avx2(sum_row, &mul, strength, rounding, &weight);
  highbd_accumulate_and_store_8_avx2(sum_row, y_pre, y_count, y_accum);
}

static void vp9_highbd_apply_temporal_filter_luma_avx2(
    const uint16_t *y_pre, int y_pre_stride, unsigned int block_width,
    unsigned int block_height, int ss_x, int ss_y, int strength,
    const int *blk_fw, int use_whole_blk, uint32_t *y_accum, uint16_t *y_count,
    const uint32_t *y_dist, const uint32_t *u_dist, const uint32_t *v_dist) {
  unsigned int blk_col = 0, uv_blk_col = 0;
  const unsigned int blk_col_step = 8, uv_blk_col_step = 8 >> ss_x;
  const unsigned int mid_width = block_width >> 1,
                     last_width = block_width - blk_col_step;
  int top_weight = blk_fw[0],
      bottom_weight = use_whole_blk ? blk_fw[0] : blk_fw[2];
  const uint32_t *const *neighbors_first = HIGHBD_LUMA_LEFT_COLUMN_NEIGHBORS;
  const uint32_t *const *neighbors_second = HIGHBD_LUMA_MIDDLE_COLUMN_NEIGHBORS;

  vp9_highbd_apply_temporal_filter_luma_8_avx2(
      y_pre + blk_col, y_pre_stride, block_height, ss_x, ss_y, strength,
      use_whole_blk, y_accum + blk_col, y_count + blk_col, y_dist + blk_col,
      u_dist + uv_blk_col, v_dist + uv_blk_col, neighbors_first,
      neighbors_second, top_weight, bottom_weight);

  blk_col += blk_col_step;
  uv_blk_col += uv_blk_col_step;

  neighbors_first = HIGHBD_LUMA_MIDDLE_COLUMN_NEIGHBORS;
  for (; blk_col < mid_width;
       blk_col += blk_col_step, uv_blk_col += uv_blk_col_step) {
    vp9_highbd_apply_temporal_filter_luma_8_avx2(
        y_pre + blk_col, y_pre_stride, block_height, ss_x, ss_y, strength,
        use_whole_blk, y_accum + blk_col, y_count + blk_col, y_dist + blk_col,
        u_dist + uv_blk_col, v_dist + uv_blk_col, neighbors_first,
        neighbors_second, top_weight, bottom_weight);
  }

  if (!use_whole_blk) {
    top_weight = blk_fw[1];
    bottom_weight = blk_fw[3];
  }

  for (; blk_col < last_width;
       blk_col += blk_col_step, uv_blk_col += uv_blk_col_step) {
    vp9_highbd_apply_temporal_filter_luma_8_avx2(
        y_pre + blk_col, y_pre_stride, block_height, ss_x, ss_y, strength,
        use_whole_blk, y_accum + blk_col, y_count + blk_col, y_dist + blk_col,
        u_dist + uv_blk_col, v_dist + uv_blk_col, neighbors_first,
        neighbors_second, top_weight, bottom_weight);
  }

  neighbors_second = HIGHBD_LUMA_RIGHT_COLUMN_NEIGHBORS;
  vp9_highbd_apply_temporal_filter_luma_8_avx2(
      y_pre + blk_col, y_pre_stride, block_height, ss_x, ss_y, strength,
      use_whole_blk, y_accum + blk_col, y_count + blk_col, y_dist + blk_col,
      u_dist + uv_blk_col, v_dist + uv_blk_col, neighbors_first,
      neighbors_second, top_weight, bottom_weight);
}

static void vp9_highbd_apply_temporal_filter_chroma_8_avx2(
    const uint16_t *u_pre, const uint16_t *v_pre, int uv_pre_stride,
    unsigned int uv_block_height, int ss_x, int ss_y, int strength,
    uint32_t *u_accum, uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count,
    const uint32_t *y_dist, const uint32_t *u_dist, const uint32_t *v_dist,
    const uint32_t *const *neighbors_fst, const uint32_t *const *neighbors_snd,
    int top_weight, int bottom_weight, const int *blk_fw) {
  const int rounding = (1 << strength) >> 1;
  __m256i weight;
  __m256i mul;
  __m256i u_sum_row_1, u_sum_row_2, u_sum_row_3;
  __m256i v_sum_row_1, v_sum_row_2, v_sum_row_3;
  __m256i u_sum_row, v_sum_row;
  unsigned int h;

  if (blk_fw) {
    weight = combine_128_pair_hbd(_mm_set1_epi32(blk_fw[0]),
                                  _mm_set1_epi32(blk_fw[1]));
  } else {
    weight = _mm256_set1_epi32(top_weight);
  }

  mul = combine_128_pair_hbd(_mm_load_si128((const __m128i *)neighbors_fst[0]),
                             _mm_load_si128((const __m128i *)neighbors_snd[0]));

  u_sum_row_2 = highbd_get_sum_8_avx2(u_dist);
  u_sum_row_3 = highbd_get_sum_8_avx2(u_dist + DIST_STRIDE);
  u_sum_row = _mm256_add_epi32(u_sum_row_2, u_sum_row_3);

  v_sum_row_2 = highbd_get_sum_8_avx2(v_dist);
  v_sum_row_3 = highbd_get_sum_8_avx2(v_dist + DIST_STRIDE);
  v_sum_row = _mm256_add_epi32(v_sum_row_2, v_sum_row_3);

  highbd_add_luma_dist_to_8_chroma_mod_avx2(y_dist, ss_x, ss_y, &u_sum_row,
                                            &v_sum_row);

  u_sum_row =
      highbd_average_8_avx2(u_sum_row, &mul, strength, rounding, &weight);
  v_sum_row =
      highbd_average_8_avx2(v_sum_row, &mul, strength, rounding, &weight);
  highbd_accumulate_and_store_8_avx2(u_sum_row, u_pre, u_count, u_accum);
  highbd_accumulate_and_store_8_avx2(v_sum_row, v_pre, v_count, v_accum);

  u_pre += uv_pre_stride;
  u_dist += DIST_STRIDE;
  v_pre += uv_pre_stride;
  v_dist += DIST_STRIDE;
  u_count += uv_pre_stride;
  u_accum += uv_pre_stride;
  v_count += uv_pre_stride;
  v_accum += uv_pre_stride;
  y_dist += DIST_STRIDE * (1 + ss_y);

  mul = combine_128_pair_hbd(_mm_load_si128((const __m128i *)neighbors_fst[1]),
                             _mm_load_si128((const __m128i *)neighbors_snd[1]));

  for (h = 1; h < uv_block_height - 1; ++h) {
    if (h == uv_block_height / 2) {
      if (blk_fw) {
        weight = combine_128_pair_hbd(_mm_set1_epi32(blk_fw[2]),
                                      _mm_set1_epi32(blk_fw[3]));
      } else {
        weight = _mm256_set1_epi32(bottom_weight);
      }
    }
    u_sum_row_1 = u_sum_row_2;
    u_sum_row_2 = u_sum_row_3;
    v_sum_row_1 = v_sum_row_2;
    v_sum_row_2 = v_sum_row_3;

    u_sum_row = _mm256_add_epi32(u_sum_row_1, u_sum_row_2);
    u_sum_row_3 = highbd_get_sum_8_avx2(u_dist + DIST_STRIDE);
    u_sum_row = _mm256_add_epi32(u_sum_row, u_sum_row_3);

    v_sum_row = _mm256_add_epi32(v_sum_row_1, v_sum_row_2);
    v_sum_row_3 = highbd_get_sum_8_avx2(v_dist + DIST_STRIDE);
    v_sum_row = _mm256_add_epi32(v_sum_row, v_sum_row_3);

    highbd_add_luma_dist_to_8_chroma_mod_avx2(y_dist, ss_x, ss_y, &u_sum_row,
                                              &v_sum_row);

    u_sum_row =
        highbd_average_8_avx2(u_sum_row, &mul, strength, rounding, &weight);
    v_sum_row =
        highbd_average_8_avx2(v_sum_row, &mul, strength, rounding, &weight);
    highbd_accumulate_and_store_8_avx2(u_sum_row, u_pre, u_count, u_accum);
    highbd_accumulate_and_store_8_avx2(v_sum_row, v_pre, v_count, v_accum);

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

  mul = combine_128_pair_hbd(_mm_load_si128((const __m128i *)neighbors_fst[0]),
                             _mm_load_si128((const __m128i *)neighbors_snd[0]));

  u_sum_row = _mm256_add_epi32(u_sum_row_2, u_sum_row_3);
  v_sum_row = _mm256_add_epi32(v_sum_row_2, v_sum_row_3);

  highbd_add_luma_dist_to_8_chroma_mod_avx2(y_dist, ss_x, ss_y, &u_sum_row,
                                            &v_sum_row);

  u_sum_row =
      highbd_average_8_avx2(u_sum_row, &mul, strength, rounding, &weight);
  v_sum_row =
      highbd_average_8_avx2(v_sum_row, &mul, strength, rounding, &weight);
  highbd_accumulate_and_store_8_avx2(u_sum_row, u_pre, u_count, u_accum);
  highbd_accumulate_and_store_8_avx2(v_sum_row, v_pre, v_count, v_accum);
}

static void vp9_highbd_apply_temporal_filter_chroma_avx2(
    const uint16_t *u_pre, const uint16_t *v_pre, int uv_pre_stride,
    unsigned int block_width, unsigned int block_height, int ss_x, int ss_y,
    int strength, const int *blk_fw, int use_whole_blk, uint32_t *u_accum,
    uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count,
    const uint32_t *y_dist, const uint32_t *u_dist, const uint32_t *v_dist) {
  const unsigned int uv_width = block_width >> ss_x;
  const unsigned int uv_height = block_height >> ss_y;
  unsigned int blk_col = 0, uv_blk_col = 0;
  const unsigned int uv_blk_col_step = 8, blk_col_step = 8 << ss_x;
  const unsigned int uv_mid_width = uv_width >> 1,
                     uv_last_width = uv_width - uv_blk_col_step;
  int top_weight = blk_fw[0];
  int bottom_weight = use_whole_blk ? blk_fw[0] : blk_fw[2];
  const uint32_t *const *neighbors_fst;
  const uint32_t *const *neighbors_snd;

  if (uv_width == 8) {
    assert(ss_x);
    if (ss_y) {
      neighbors_fst = HIGHBD_CHROMA_DOUBLE_SS_LEFT_COLUMN_NEIGHBORS;
      neighbors_snd = HIGHBD_CHROMA_DOUBLE_SS_RIGHT_COLUMN_NEIGHBORS;
    } else {
      neighbors_fst = HIGHBD_CHROMA_SINGLE_SS_LEFT_COLUMN_NEIGHBORS;
      neighbors_snd = HIGHBD_CHROMA_SINGLE_SS_RIGHT_COLUMN_NEIGHBORS;
    }
    vp9_highbd_apply_temporal_filter_chroma_8_avx2(
        u_pre, v_pre, uv_pre_stride, uv_height, ss_x, ss_y, strength, u_accum,
        u_count, v_accum, v_count, y_dist, u_dist, v_dist, neighbors_fst,
        neighbors_snd, top_weight, bottom_weight,
        use_whole_blk ? NULL : blk_fw);
    return;
  }

  if (ss_x && ss_y) {
    neighbors_fst = HIGHBD_CHROMA_DOUBLE_SS_LEFT_COLUMN_NEIGHBORS;
    neighbors_snd = HIGHBD_CHROMA_DOUBLE_SS_MIDDLE_COLUMN_NEIGHBORS;
  } else if (ss_x || ss_y) {
    neighbors_fst = HIGHBD_CHROMA_SINGLE_SS_LEFT_COLUMN_NEIGHBORS;
    neighbors_snd = HIGHBD_CHROMA_SINGLE_SS_MIDDLE_COLUMN_NEIGHBORS;
  } else {
    neighbors_fst = HIGHBD_CHROMA_NO_SS_LEFT_COLUMN_NEIGHBORS;
    neighbors_snd = HIGHBD_CHROMA_NO_SS_MIDDLE_COLUMN_NEIGHBORS;
  }

  vp9_highbd_apply_temporal_filter_chroma_8_avx2(
      u_pre + uv_blk_col, v_pre + uv_blk_col, uv_pre_stride, uv_height, ss_x,
      ss_y, strength, u_accum + uv_blk_col, u_count + uv_blk_col,
      v_accum + uv_blk_col, v_count + uv_blk_col, y_dist + blk_col,
      u_dist + uv_blk_col, v_dist + uv_blk_col, neighbors_fst, neighbors_snd,
      top_weight, bottom_weight, NULL);

  blk_col += blk_col_step;
  uv_blk_col += uv_blk_col_step;

  if (ss_x && ss_y) {
    neighbors_fst = HIGHBD_CHROMA_DOUBLE_SS_MIDDLE_COLUMN_NEIGHBORS;
  } else if (ss_x || ss_y) {
    neighbors_fst = HIGHBD_CHROMA_SINGLE_SS_MIDDLE_COLUMN_NEIGHBORS;
  } else {
    neighbors_fst = HIGHBD_CHROMA_NO_SS_MIDDLE_COLUMN_NEIGHBORS;
  }

  for (; uv_blk_col < uv_mid_width;
       blk_col += blk_col_step, uv_blk_col += uv_blk_col_step) {
    vp9_highbd_apply_temporal_filter_chroma_8_avx2(
        u_pre + uv_blk_col, v_pre + uv_blk_col, uv_pre_stride, uv_height, ss_x,
        ss_y, strength, u_accum + uv_blk_col, u_count + uv_blk_col,
        v_accum + uv_blk_col, v_count + uv_blk_col, y_dist + blk_col,
        u_dist + uv_blk_col, v_dist + uv_blk_col, neighbors_fst, neighbors_snd,
        top_weight, bottom_weight, NULL);
  }

  if (!use_whole_blk) {
    top_weight = blk_fw[1];
    bottom_weight = blk_fw[3];
  }

  for (; uv_blk_col < uv_last_width;
       blk_col += blk_col_step, uv_blk_col += uv_blk_col_step) {
    vp9_highbd_apply_temporal_filter_chroma_8_avx2(
        u_pre + uv_blk_col, v_pre + uv_blk_col, uv_pre_stride, uv_height, ss_x,
        ss_y, strength, u_accum + uv_blk_col, u_count + uv_blk_col,
        v_accum + uv_blk_col, v_count + uv_blk_col, y_dist + blk_col,
        u_dist + uv_blk_col, v_dist + uv_blk_col, neighbors_fst, neighbors_snd,
        top_weight, bottom_weight, NULL);
  }

  if (ss_x && ss_y) {
    neighbors_snd = HIGHBD_CHROMA_DOUBLE_SS_RIGHT_COLUMN_NEIGHBORS;
  } else if (ss_x || ss_y) {
    neighbors_snd = HIGHBD_CHROMA_SINGLE_SS_RIGHT_COLUMN_NEIGHBORS;
  } else {
    neighbors_snd = HIGHBD_CHROMA_NO_SS_RIGHT_COLUMN_NEIGHBORS;
  }

  vp9_highbd_apply_temporal_filter_chroma_8_avx2(
      u_pre + uv_blk_col, v_pre + uv_blk_col, uv_pre_stride, uv_height, ss_x,
      ss_y, strength, u_accum + uv_blk_col, u_count + uv_blk_col,
      v_accum + uv_blk_col, v_count + uv_blk_col, y_dist + blk_col,
      u_dist + uv_blk_col, v_dist + uv_blk_col, neighbors_fst, neighbors_snd,
      top_weight, bottom_weight, NULL);
}

void vp9_highbd_apply_temporal_filter_avx2(
    const uint16_t *y_src, int y_src_stride, const uint16_t *y_pre,
    int y_pre_stride, const uint16_t *u_src, const uint16_t *v_src,
    int uv_src_stride, const uint16_t *u_pre, const uint16_t *v_pre,
    int uv_pre_stride, unsigned int block_width, unsigned int block_height,
    int ss_x, int ss_y, int strength, const int *const blk_fw,
    int use_whole_blk, uint32_t *y_accum, uint16_t *y_count, uint32_t *u_accum,
    uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count) {
  const unsigned int chroma_height = block_height >> ss_y;
  const unsigned int chroma_width = block_width >> ss_x;

  DECLARE_ALIGNED(32, uint32_t, y_dist[BH * DIST_STRIDE]) = { 0 };
  DECLARE_ALIGNED(32, uint32_t, u_dist[BH * DIST_STRIDE]) = { 0 };
  DECLARE_ALIGNED(32, uint32_t, v_dist[BH * DIST_STRIDE]) = { 0 };

  uint32_t *y_dist_ptr = y_dist + 1;
  uint32_t *u_dist_ptr = u_dist + 1;
  uint32_t *v_dist_ptr = v_dist + 1;
  const uint16_t *y_src_ptr = y_src, *u_src_ptr = u_src, *v_src_ptr = v_src;
  const uint16_t *y_pre_ptr = y_pre, *u_pre_ptr = u_pre, *v_pre_ptr = v_pre;
  unsigned int row, blk_col;

  assert(block_width <= BW && "block width too large");
  assert(block_height <= BH && "block height too large");
  assert(block_width % 16 == 0 && "block width must be multiple of 16");
  assert(block_height % 2 == 0 && "block height must be even");
  assert((ss_x == 0 || ss_x == 1) && (ss_y == 0 || ss_y == 1) &&
         "invalid chroma subsampling");
  assert(strength >= 4 && strength <= 14 &&
         "invalid adjusted temporal filter strength");
  assert(blk_fw[0] >= 0 && "filter weight must be positive");
  assert(
      (use_whole_blk || (blk_fw[1] >= 0 && blk_fw[2] >= 0 && blk_fw[3] >= 0)) &&
      "subblock filter weight must be positive");
  assert(blk_fw[0] <= 2 && "sublock filter weight must be less than 2");
  assert(
      (use_whole_blk || (blk_fw[1] <= 2 && blk_fw[2] <= 2 && blk_fw[3] <= 2)) &&
      "subblock filter weight must be less than 2");

  for (row = 0; row < block_height; row++) {
    for (blk_col = 0; blk_col < block_width; blk_col += 8) {
      highbd_store_dist_8_avx2(y_src_ptr + blk_col, y_pre_ptr + blk_col,
                               y_dist_ptr + blk_col);
    }
    y_src_ptr += y_src_stride;
    y_pre_ptr += y_pre_stride;
    y_dist_ptr += DIST_STRIDE;
  }

  for (row = 0; row < chroma_height; row++) {
    for (blk_col = 0; blk_col < chroma_width; blk_col += 8) {
      highbd_store_dist_8_avx2(u_src_ptr + blk_col, u_pre_ptr + blk_col,
                               u_dist_ptr + blk_col);
      highbd_store_dist_8_avx2(v_src_ptr + blk_col, v_pre_ptr + blk_col,
                               v_dist_ptr + blk_col);
    }
    u_src_ptr += uv_src_stride;
    u_pre_ptr += uv_pre_stride;
    u_dist_ptr += DIST_STRIDE;
    v_src_ptr += uv_src_stride;
    v_pre_ptr += uv_pre_stride;
    v_dist_ptr += DIST_STRIDE;
  }

  y_dist_ptr = y_dist + 1;
  u_dist_ptr = u_dist + 1;
  v_dist_ptr = v_dist + 1;

  vp9_highbd_apply_temporal_filter_luma_avx2(
      y_pre, y_pre_stride, block_width, block_height, ss_x, ss_y, strength,
      blk_fw, use_whole_blk, y_accum, y_count, y_dist_ptr, u_dist_ptr,
      v_dist_ptr);

  vp9_highbd_apply_temporal_filter_chroma_avx2(
      u_pre, v_pre, uv_pre_stride, block_width, block_height, ss_x, ss_y,
      strength, blk_fw, use_whole_blk, u_accum, u_count, v_accum, v_count,
      y_dist_ptr, u_dist_ptr, v_dist_ptr);
}
