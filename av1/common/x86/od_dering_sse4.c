/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <smmintrin.h>
#include <emmintrin.h>
#include <tmmintrin.h>

#include "./av1_rtcd.h"
#include "av1/common/x86/od_dering_sse4.h"

/* partial A is a 16-bit vector of the form:
   [x8 x7 x6 x5 x4 x3 x2 x1] and partial B has the form:
   [0  y1 y2 y3 y4 y5 y6 y7].
   This function computes (x1^2+y1^2)*C1 + (x2^2+y2^2)*C2 + ...
   (x7^2+y2^7)*C7 + (x8^2+0^2)*C8 where the C1..C8 constants are in const1
   and const2. */
static INLINE __m128i fold_mul_and_sum(__m128i partiala, __m128i partialb,
                                       __m128i const1, __m128i const2) {
  __m128i tmp;
  /* Reverse partial B. */
  partialb = _mm_shuffle_epi8(
      partialb,
      _mm_set_epi8(15, 14, 1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12));
  /* Interleave the x and y values of identical indices and pair x8 with 0. */
  tmp = partiala;
  partiala = _mm_unpacklo_epi16(partiala, partialb);
  partialb = _mm_unpackhi_epi16(tmp, partialb);
  /* Square and add the corresponding x and y values. */
  partiala = _mm_madd_epi16(partiala, partiala);
  partialb = _mm_madd_epi16(partialb, partialb);
  /* Multiply by constant. */
  partiala = _mm_mullo_epi32(partiala, const1);
  partialb = _mm_mullo_epi32(partialb, const2);
  /* Sum all results. */
  partiala = _mm_add_epi32(partiala, partialb);
  return partiala;
}

static INLINE __m128i hsum4(__m128i x0, __m128i x1, __m128i x2, __m128i x3) {
  __m128i t0, t1, t2, t3;
  t0 = _mm_unpacklo_epi32(x0, x1);
  t1 = _mm_unpacklo_epi32(x2, x3);
  t2 = _mm_unpackhi_epi32(x0, x1);
  t3 = _mm_unpackhi_epi32(x2, x3);
  x0 = _mm_unpacklo_epi64(t0, t1);
  x1 = _mm_unpackhi_epi64(t0, t1);
  x2 = _mm_unpacklo_epi64(t2, t3);
  x3 = _mm_unpackhi_epi64(t2, t3);
  return _mm_add_epi32(_mm_add_epi32(x0, x1), _mm_add_epi32(x2, x3));
}

/* Horizontal sum of 8x16-bit unsigned values. */
static INLINE int32_t hsum_epi16(__m128i a) {
  a = _mm_madd_epi16(a, _mm_set1_epi16(1));
  a = _mm_hadd_epi32(a, a);
  a = _mm_hadd_epi32(a, a);
  return _mm_cvtsi128_si32(a);
}

/* Computes cost for directions 0, 5, 6 and 7. We can call this function again
   to compute the remaining directions. */
static INLINE __m128i compute_directions(__m128i lines[8],
                                         int32_t tmp_cost1[4]) {
  __m128i partial4a, partial4b, partial5a, partial5b, partial7a, partial7b;
  __m128i partial6;
  __m128i tmp;
  /* Partial sums for lines 0 and 1. */
  partial4a = _mm_slli_si128(lines[0], 14);
  partial4b = _mm_srli_si128(lines[0], 2);
  partial4a = _mm_add_epi16(partial4a, _mm_slli_si128(lines[1], 12));
  partial4b = _mm_add_epi16(partial4b, _mm_srli_si128(lines[1], 4));
  tmp = _mm_add_epi16(lines[0], lines[1]);
  partial5a = _mm_slli_si128(tmp, 10);
  partial5b = _mm_srli_si128(tmp, 6);
  partial7a = _mm_slli_si128(tmp, 4);
  partial7b = _mm_srli_si128(tmp, 12);
  partial6 = tmp;

  /* Partial sums for lines 2 and 3. */
  partial4a = _mm_add_epi16(partial4a, _mm_slli_si128(lines[2], 10));
  partial4b = _mm_add_epi16(partial4b, _mm_srli_si128(lines[2], 6));
  partial4a = _mm_add_epi16(partial4a, _mm_slli_si128(lines[3], 8));
  partial4b = _mm_add_epi16(partial4b, _mm_srli_si128(lines[3], 8));
  tmp = _mm_add_epi16(lines[2], lines[3]);
  partial5a = _mm_add_epi16(partial5a, _mm_slli_si128(tmp, 8));
  partial5b = _mm_add_epi16(partial5b, _mm_srli_si128(tmp, 8));
  partial7a = _mm_add_epi16(partial7a, _mm_slli_si128(tmp, 6));
  partial7b = _mm_add_epi16(partial7b, _mm_srli_si128(tmp, 10));
  partial6 = _mm_add_epi16(partial6, tmp);

  /* Partial sums for lines 4 and 5. */
  partial4a = _mm_add_epi16(partial4a, _mm_slli_si128(lines[4], 6));
  partial4b = _mm_add_epi16(partial4b, _mm_srli_si128(lines[4], 10));
  partial4a = _mm_add_epi16(partial4a, _mm_slli_si128(lines[5], 4));
  partial4b = _mm_add_epi16(partial4b, _mm_srli_si128(lines[5], 12));
  tmp = _mm_add_epi16(lines[4], lines[5]);
  partial5a = _mm_add_epi16(partial5a, _mm_slli_si128(tmp, 6));
  partial5b = _mm_add_epi16(partial5b, _mm_srli_si128(tmp, 10));
  partial7a = _mm_add_epi16(partial7a, _mm_slli_si128(tmp, 8));
  partial7b = _mm_add_epi16(partial7b, _mm_srli_si128(tmp, 8));
  partial6 = _mm_add_epi16(partial6, tmp);

  /* Partial sums for lines 6 and 7. */
  partial4a = _mm_add_epi16(partial4a, _mm_slli_si128(lines[6], 2));
  partial4b = _mm_add_epi16(partial4b, _mm_srli_si128(lines[6], 14));
  partial4a = _mm_add_epi16(partial4a, lines[7]);
  tmp = _mm_add_epi16(lines[6], lines[7]);
  partial5a = _mm_add_epi16(partial5a, _mm_slli_si128(tmp, 4));
  partial5b = _mm_add_epi16(partial5b, _mm_srli_si128(tmp, 12));
  partial7a = _mm_add_epi16(partial7a, _mm_slli_si128(tmp, 10));
  partial7b = _mm_add_epi16(partial7b, _mm_srli_si128(tmp, 6));
  partial6 = _mm_add_epi16(partial6, tmp);

  /* Compute costs in terms of partial sums. */
  partial4a =
      fold_mul_and_sum(partial4a, partial4b, _mm_set_epi32(210, 280, 420, 840),
                       _mm_set_epi32(105, 120, 140, 168));
  partial7a =
      fold_mul_and_sum(partial7a, partial7b, _mm_set_epi32(210, 420, 0, 0),
                       _mm_set_epi32(105, 105, 105, 140));
  partial5a =
      fold_mul_and_sum(partial5a, partial5b, _mm_set_epi32(210, 420, 0, 0),
                       _mm_set_epi32(105, 105, 105, 140));
  partial6 = _mm_madd_epi16(partial6, partial6);
  partial6 = _mm_mullo_epi32(partial6, _mm_set1_epi32(105));

  partial4a = hsum4(partial4a, partial5a, partial6, partial7a);
  _mm_storeu_si128((__m128i *)tmp_cost1, partial4a);
  return partial4a;
}

/* transpose and reverse the order of the lines -- equivalent to a 90-degree
   counter-clockwise rotation of the pixels. */
static INLINE void array_reverse_transpose_8x8(__m128i *in, __m128i *res) {
  const __m128i tr0_0 = _mm_unpacklo_epi16(in[0], in[1]);
  const __m128i tr0_1 = _mm_unpacklo_epi16(in[2], in[3]);
  const __m128i tr0_2 = _mm_unpackhi_epi16(in[0], in[1]);
  const __m128i tr0_3 = _mm_unpackhi_epi16(in[2], in[3]);
  const __m128i tr0_4 = _mm_unpacklo_epi16(in[4], in[5]);
  const __m128i tr0_5 = _mm_unpacklo_epi16(in[6], in[7]);
  const __m128i tr0_6 = _mm_unpackhi_epi16(in[4], in[5]);
  const __m128i tr0_7 = _mm_unpackhi_epi16(in[6], in[7]);

  const __m128i tr1_0 = _mm_unpacklo_epi32(tr0_0, tr0_1);
  const __m128i tr1_1 = _mm_unpacklo_epi32(tr0_4, tr0_5);
  const __m128i tr1_2 = _mm_unpackhi_epi32(tr0_0, tr0_1);
  const __m128i tr1_3 = _mm_unpackhi_epi32(tr0_4, tr0_5);
  const __m128i tr1_4 = _mm_unpacklo_epi32(tr0_2, tr0_3);
  const __m128i tr1_5 = _mm_unpacklo_epi32(tr0_6, tr0_7);
  const __m128i tr1_6 = _mm_unpackhi_epi32(tr0_2, tr0_3);
  const __m128i tr1_7 = _mm_unpackhi_epi32(tr0_6, tr0_7);

  res[7] = _mm_unpacklo_epi64(tr1_0, tr1_1);
  res[6] = _mm_unpackhi_epi64(tr1_0, tr1_1);
  res[5] = _mm_unpacklo_epi64(tr1_2, tr1_3);
  res[4] = _mm_unpackhi_epi64(tr1_2, tr1_3);
  res[3] = _mm_unpacklo_epi64(tr1_4, tr1_5);
  res[2] = _mm_unpackhi_epi64(tr1_4, tr1_5);
  res[1] = _mm_unpacklo_epi64(tr1_6, tr1_7);
  res[0] = _mm_unpackhi_epi64(tr1_6, tr1_7);
}

int od_dir_find8_sse4_1(const od_dering_in *img, int stride, int32_t *var,
                        int coeff_shift) {
  int i;
  int32_t cost[8];
  int32_t best_cost = 0;
  int best_dir = 0;
  __m128i lines[8];
  __m128i dir03, dir47;
  __m128i max;
  for (i = 0; i < 8; i++) {
    lines[i] = _mm_loadu_si128((__m128i *)&img[i * stride]);
    lines[i] = _mm_sub_epi16(_mm_srai_epi16(lines[i], coeff_shift),
                             _mm_set1_epi16(128));
  }

  /* Compute "mostly vertical" directions. */
  dir47 = compute_directions(lines, cost + 4);

  array_reverse_transpose_8x8(lines, lines);

  /* Compute "mostly horizontal" directions. */
  dir03 = compute_directions(lines, cost);

#if 1
  max = _mm_max_epi32(dir03, dir47);
  max = _mm_max_epi32(max, _mm_shuffle_epi32(max, _MM_SHUFFLE(1, 0, 3, 2)));
  max = _mm_max_epi32(max, _mm_shuffle_epi32(max, _MM_SHUFFLE(2, 3, 0, 1)));
  dir03 = _mm_and_si128(_mm_cmpeq_epi32(max, dir03),
                        _mm_setr_epi32(-1, -2, -3, -4));
  dir47 = _mm_and_si128(_mm_cmpeq_epi32(max, dir47),
                        _mm_setr_epi32(-5, -6, -7, -8));
  dir03 = _mm_max_epu32(dir03, dir47);
  dir03 = _mm_max_epu32(dir03, _mm_unpackhi_epi64(dir03, dir03));
  dir03 =
      _mm_max_epu32(dir03, _mm_shufflelo_epi16(dir03, _MM_SHUFFLE(1, 0, 3, 2)));
  dir03 = _mm_xor_si128(dir03, _mm_set1_epi32(0xFFFFFFFF));

  best_dir = _mm_cvtsi128_si32(dir03);
  best_cost = _mm_cvtsi128_si32(max);
#else
  for (i = 0; i < 8; i++) {
    if (cost[i] > best_cost) {
      best_cost = cost[i];
      best_dir = i;
    }
  }
#endif
  /* Difference between the optimal variance and the variance along the
     orthogonal direction. Again, the sum(x^2) terms cancel out. */
  *var = best_cost - cost[(best_dir + 4) & 7];
  /* We'd normally divide by 840, but dividing by 1024 is close enough
     for what we're going to do with this. */
  *var >>= 10;
  return best_dir;
}

static INLINE __m128i od_cmplt_abs_epi16(__m128i in, __m128i threshold) {
  return _mm_cmplt_epi16(_mm_abs_epi16(in), threshold);
}

int od_filter_dering_direction_4x4_sse4_1(int16_t *y, int ystride,
                                          const int16_t *in, int threshold,
                                          int dir) {
  int i;
  __m128i sum;
  __m128i p;
  __m128i cmp;
  __m128i row;
  __m128i res;
  __m128i tmp;
  __m128i thresh;
  __m128i total_abs;
  int off1, off2;
  off1 = OD_DIRECTION_OFFSETS_TABLE[dir][0];
  off2 = OD_DIRECTION_OFFSETS_TABLE[dir][1];
  total_abs = _mm_setzero_si128();
  thresh = _mm_set1_epi16(threshold);
  for (i = 0; i < 4; i += 2) {
    sum = _mm_set1_epi16(0);
    row = _mm_unpacklo_epi64(
        _mm_loadl_epi64((__m128i *)&in[i * OD_FILT_BSTRIDE]),
        _mm_loadl_epi64((__m128i *)&in[(i + 1) * OD_FILT_BSTRIDE]));

    /*p = in[i*OD_FILT_BSTRIDE + offset] - row*/
    tmp = _mm_unpacklo_epi64(
        _mm_loadl_epi64((__m128i *)&in[i * OD_FILT_BSTRIDE + off1]),
        _mm_loadl_epi64((__m128i *)&in[(i + 1) * OD_FILT_BSTRIDE + off1]));
    p = _mm_sub_epi16(tmp, row);
    /*if (abs(p) < thresh) sum += taps[k]*p*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_slli_epi16(p, 2);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);
    /*p = in[i*OD_FILT_BSTRIDE - offset] - row*/
    tmp = _mm_unpacklo_epi64(
        _mm_loadl_epi64((__m128i *)&in[i * OD_FILT_BSTRIDE - off1]),
        _mm_loadl_epi64((__m128i *)&in[(i + 1) * OD_FILT_BSTRIDE - off1]));
    p = _mm_sub_epi16(tmp, row);
    /*if (abs(p) < thresh) sum += taps[k]*p1*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_slli_epi16(p, 2);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);

    /*p = in[i*OD_FILT_BSTRIDE + offset] - row*/
    tmp = _mm_unpacklo_epi64(
        _mm_loadl_epi64((__m128i *)&in[i * OD_FILT_BSTRIDE + off2]),
        _mm_loadl_epi64((__m128i *)&in[(i + 1) * OD_FILT_BSTRIDE + off2]));
    p = _mm_sub_epi16(tmp, row);
    /*if (abs(p) < thresh) sum += taps[k]*p*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);
    /*p = in[i*OD_FILT_BSTRIDE - offset] - row*/
    tmp = _mm_unpacklo_epi64(
        _mm_loadl_epi64((__m128i *)&in[i * OD_FILT_BSTRIDE - off2]),
        _mm_loadl_epi64((__m128i *)&in[(i + 1) * OD_FILT_BSTRIDE - off2]));
    p = _mm_sub_epi16(tmp, row);
    /*if (abs(p) < thresh) sum += taps[k]*p1*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);

    /*res = row + ((sum + 8) >> 4)*/
    res = _mm_add_epi16(sum, _mm_set1_epi16(8));
    res = _mm_srai_epi16(res, 4);
    total_abs = _mm_add_epi16(total_abs, _mm_abs_epi16(res));
    res = _mm_add_epi16(row, res);
    _mm_storel_epi64((__m128i *)&y[i * ystride], res);
    _mm_storel_epi64((__m128i *)&y[(i + 1) * ystride],
                     _mm_unpackhi_epi64(res, res));
  }
  return (hsum_epi16(total_abs) + 2) >> 2;
}

int od_filter_dering_direction_8x8_sse4_1(int16_t *y, int ystride,
                                          const int16_t *in, int threshold,
                                          int dir) {
  int i;
  __m128i sum;
  __m128i p;
  __m128i cmp;
  __m128i row;
  __m128i res;
  __m128i thresh;
  __m128i total_abs;
  int off1, off2, off3;
  off1 = OD_DIRECTION_OFFSETS_TABLE[dir][0];
  off2 = OD_DIRECTION_OFFSETS_TABLE[dir][1];
  off3 = OD_DIRECTION_OFFSETS_TABLE[dir][2];
  total_abs = _mm_setzero_si128();
  thresh = _mm_set1_epi16(threshold);
  for (i = 0; i < 8; i++) {
    sum = _mm_set1_epi16(0);
    row = _mm_loadu_si128((__m128i *)&in[i * OD_FILT_BSTRIDE]);

    /*p = in[i*OD_FILT_BSTRIDE + offset] - row*/
    p = _mm_sub_epi16(
        _mm_loadu_si128((__m128i *)&in[i * OD_FILT_BSTRIDE + off1]), row);
    /*if (abs(p) < thresh) sum += taps[k]*p*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_add_epi16(p, _mm_slli_epi16(p, 1));
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);

    /*p = in[i*OD_FILT_BSTRIDE - offset] - row*/
    p = _mm_sub_epi16(
        _mm_loadu_si128((__m128i *)&in[i * OD_FILT_BSTRIDE - off1]), row);
    /*if (abs(p) < thresh) sum += taps[k]*p1*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_add_epi16(p, _mm_slli_epi16(p, 1));
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);

    /*p = in[i*OD_FILT_BSTRIDE + offset] - row*/
    p = _mm_sub_epi16(
        _mm_loadu_si128((__m128i *)&in[i * OD_FILT_BSTRIDE + off2]), row);
    /*if (abs(p) < thresh) sum += taps[k]*p*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_slli_epi16(p, 1);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);

    /*p = in[i*OD_FILT_BSTRIDE - offset] - row*/
    p = _mm_sub_epi16(
        _mm_loadu_si128((__m128i *)&in[i * OD_FILT_BSTRIDE - off2]), row);
    /*if (abs(p) < thresh) sum += taps[k]*p1*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_slli_epi16(p, 1);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);

    /*p = in[i*OD_FILT_BSTRIDE + offset] - row*/
    p = _mm_sub_epi16(
        _mm_loadu_si128((__m128i *)&in[i * OD_FILT_BSTRIDE + off3]), row);
    /*if (abs(p) < thresh) sum += taps[k]*p*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);

    /*p = in[i*OD_FILT_BSTRIDE - offset] - row*/
    p = _mm_sub_epi16(
        _mm_loadu_si128((__m128i *)&in[i * OD_FILT_BSTRIDE - off3]), row);
    /*if (abs(p) < thresh) sum += taps[k]*p1*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);

    /*res = row + ((sum + 8) >> 4)*/
    res = _mm_add_epi16(sum, _mm_set1_epi16(8));
    res = _mm_srai_epi16(res, 4);
    total_abs = _mm_add_epi16(total_abs, _mm_abs_epi16(res));
    res = _mm_add_epi16(row, res);
    _mm_storeu_si128((__m128i *)&y[i * ystride], res);
  }
  return (hsum_epi16(total_abs) + 8) >> 4;
}

void od_filter_dering_orthogonal_4x4_sse4_1(int16_t *y, int ystride,
                                            const int16_t *in, int threshold,
                                            int dir) {
  int i;
  int offset;
  __m128i res;
  __m128i p;
  __m128i cmp;
  __m128i row;
  __m128i sum;
  __m128i tmp;
  __m128i thresh;
  thresh = _mm_set1_epi16(threshold);
  if (dir > 0 && dir < 4)
    offset = OD_FILT_BSTRIDE;
  else
    offset = 1;
  for (i = 0; i < 4; i += 2) {
    sum = _mm_set1_epi16(0);
    row = _mm_unpacklo_epi64(
        _mm_loadl_epi64((__m128i *)&in[i * OD_FILT_BSTRIDE]),
        _mm_loadl_epi64((__m128i *)&in[(i + 1) * OD_FILT_BSTRIDE]));

    /*p = in[i*OD_FILT_BSTRIDE + k*offset] - row*/
    tmp = _mm_unpacklo_epi64(
        _mm_loadl_epi64((__m128i *)&in[i * OD_FILT_BSTRIDE + offset]),
        _mm_loadl_epi64((__m128i *)&in[(i + 1) * OD_FILT_BSTRIDE + offset]));
    p = _mm_sub_epi16(tmp, row);
    /*if (abs(p) < threshold) sum += p*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);
    /*p = in[i*OD_FILT_BSTRIDE - k*offset] - row*/
    tmp = _mm_unpacklo_epi64(
        _mm_loadl_epi64((__m128i *)&in[i * OD_FILT_BSTRIDE - offset]),
        _mm_loadl_epi64((__m128i *)&in[(i + 1) * OD_FILT_BSTRIDE - offset]));
    p = _mm_sub_epi16(tmp, row);
    /*if (abs(p) < threshold) sum += p*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);

    /*row + ((5*sum + 8) >> 4)*/
    res = _mm_mullo_epi16(sum, _mm_set1_epi16(5));
    res = _mm_add_epi16(res, _mm_set1_epi16(8));
    res = _mm_srai_epi16(res, 4);
    res = _mm_add_epi16(res, row);
    _mm_storel_epi64((__m128i *)&y[i * ystride], res);
    _mm_storel_epi64((__m128i *)&y[(i + 1) * ystride],
                     _mm_unpackhi_epi64(res, res));
  }
}

void od_filter_dering_orthogonal_8x8_sse4_1(int16_t *y, int ystride,
                                            const int16_t *in, int threshold,
                                            int dir) {
  int i;
  int offset;
  __m128i res;
  __m128i p;
  __m128i cmp;
  __m128i row;
  __m128i sum;
  __m128i thresh;
  thresh = _mm_set1_epi16(threshold);
  if (dir > 0 && dir < 4)
    offset = OD_FILT_BSTRIDE;
  else
    offset = 1;
  for (i = 0; i < 8; i++) {
    sum = _mm_set1_epi16(0);
    row = _mm_loadu_si128((__m128i *)&in[i * OD_FILT_BSTRIDE]);

    /*p = in[i*OD_FILT_BSTRIDE + k*offset] - row*/
    p = _mm_sub_epi16(
        _mm_loadu_si128((__m128i *)&in[i * OD_FILT_BSTRIDE + 1 * offset]), row);
    /*if (abs(p) < thresh) sum += p*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);
    /*p = in[i*OD_FILT_BSTRIDE - k*offset] - row*/
    p = _mm_sub_epi16(
        _mm_loadu_si128((__m128i *)&in[i * OD_FILT_BSTRIDE - 1 * offset]), row);
    /*if (abs(p) < threshold) sum += p*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);

    /*p = in[i*OD_FILT_BSTRIDE + k*offset] - row*/
    p = _mm_sub_epi16(
        _mm_loadu_si128((__m128i *)&in[i * OD_FILT_BSTRIDE + 2 * offset]), row);
    /*if (abs(p) < threshold) sum += p*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);
    /*p = in[i*OD_FILT_BSTRIDE - k*offset] - row*/
    p = _mm_sub_epi16(
        _mm_loadu_si128((__m128i *)&in[i * OD_FILT_BSTRIDE - 2 * offset]), row);
    /*if (abs(p) < threshold) sum += p*/
    cmp = od_cmplt_abs_epi16(p, thresh);
    p = _mm_and_si128(p, cmp);
    sum = _mm_add_epi16(sum, p);

    /*row + ((3*sum + 8) >> 4)*/
    res = _mm_mullo_epi16(sum, _mm_set1_epi16(3));
    res = _mm_add_epi16(res, _mm_set1_epi16(8));
    res = _mm_srai_epi16(res, 4);
    res = _mm_add_epi16(res, row);
    _mm_storeu_si128((__m128i *)&y[i * ystride], res);
  }
}
