/*
 *  Copyright (c) 2025 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <immintrin.h>

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"

static INLINE __m256i avg3_epu16_avx2(const __m256i *x, const __m256i *y,
                                      const __m256i *z) {
  const __m256i one = _mm256_set1_epi16(1);
  const __m256i a = _mm256_avg_epu16(*x, *z);
  const __m256i b =
      _mm256_subs_epu16(a, _mm256_and_si256(_mm256_xor_si256(*x, *z), one));
  return _mm256_avg_epu16(b, *y);
}

/*
 palignr in AVX2 operates in-lane
 hi: hi_hi | hi_lo
 mid hi_lo | lo_hi
 lo: lo_hi | lo_lo
*/
#define ALIGNR_256(res, hi, lo, i)                       \
  do {                                                   \
    __m256i _mid = _mm256_permute2x128_si256(hi, lo, 3); \
    res = _mm256_alignr_epi8(_mid, lo, (i));             \
  } while (0)

#define ALIGNR_256_HI(res, hi, lo, i)                    \
  do {                                                   \
    __m256i _mid = _mm256_permute2x128_si256(hi, lo, 3); \
    res = _mm256_alignr_epi8(hi, _mid, (i));             \
  } while (0)

void vpx_highbd_d63_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t stride,
                                         const uint16_t *above,
                                         const uint16_t *left, int bd) {
  (void)left;
  (void)bd;

  const __m256i A0 = _mm256_loadu_si256((const __m256i *)(above));  // 0..15
  const __m256i A1 =
      _mm256_loadu_si256((const __m256i *)(above + 16));  // 16..31

  // Replicate above[31] for right-edge extension
  const __m256i AR =
      _mm256_permute4x64_epi64(_mm256_shufflehi_epi16(A1, 0xff), 0xff);

  __m256i B0;
  ALIGNR_256(B0, A1, A0, 2);
  __m256i B1;
  ALIGNR_256(B1, AR, A1, 2);

  __m256i C0;
  ALIGNR_256(C0, A1, A0, 4);
  __m256i C1;
  ALIGNR_256(C1, AR, A1, 4);

  __m256i avg2_0 = _mm256_avg_epu16(A0, B0);
  __m256i avg2_1 = _mm256_avg_epu16(A1, B1);
  __m256i avg3_0 = avg3_epu16_avx2(&A0, &B0, &C0);
  __m256i avg3_1 = avg3_epu16_avx2(&A1, &B1, &C1);

  for (int i = 0; i < 30; i += 2) {
    _mm256_storeu_si256((__m256i *)dst, avg2_0);
    _mm256_storeu_si256((__m256i *)(dst + 16), avg2_1);
    dst += stride;

    _mm256_storeu_si256((__m256i *)dst, avg3_0);
    _mm256_storeu_si256((__m256i *)(dst + 16), avg3_1);
    dst += stride;

    ALIGNR_256(avg2_0, avg2_1, avg2_0, 2);
    ALIGNR_256(avg2_1, AR, avg2_1, 2);
    ALIGNR_256(avg3_0, avg3_1, avg3_0, 2);
    ALIGNR_256(avg3_1, AR, avg3_1, 2);
  }

  _mm256_storeu_si256((__m256i *)dst, avg2_0);
  _mm256_storeu_si256((__m256i *)(dst + 16), avg2_1);
  dst += stride;

  _mm256_storeu_si256((__m256i *)dst, avg3_0);
  _mm256_storeu_si256((__m256i *)(dst + 16), avg3_1);
}

static INLINE void d207_store_8x32_avx2(uint16_t **dst, const ptrdiff_t stride,
                                        const __m256i *ab, const __m256i *cd,
                                        const __m256i *ef) {
  _mm256_storeu_si256((__m256i *)*dst, *ab);
  _mm256_storeu_si256((__m256i *)(*dst + 16), *cd);
  *dst += stride;

  __m256i ab_shift;
  __m256i cd_shift;

  ALIGNR_256(ab_shift, *cd, *ab, 4);
  ALIGNR_256(cd_shift, *ef, *cd, 4);
  _mm256_storeu_si256((__m256i *)*dst, ab_shift);
  _mm256_storeu_si256((__m256i *)(*dst + 16), cd_shift);
  *dst += stride;

  ALIGNR_256(ab_shift, *cd, *ab, 8);
  ALIGNR_256(cd_shift, *ef, *cd, 8);
  _mm256_storeu_si256((__m256i *)*dst, ab_shift);
  _mm256_storeu_si256((__m256i *)(*dst + 16), cd_shift);
  *dst += stride;

  ALIGNR_256(ab_shift, *cd, *ab, 12);
  ALIGNR_256(cd_shift, *ef, *cd, 12);
  _mm256_storeu_si256((__m256i *)*dst, ab_shift);
  _mm256_storeu_si256((__m256i *)(*dst + 16), cd_shift);
  *dst += stride;

  ab_shift = _mm256_permute2x128_si256(*cd, *ab, 3);
  cd_shift = _mm256_permute2x128_si256(*ef, *cd, 3);
  _mm256_storeu_si256((__m256i *)*dst, ab_shift);
  _mm256_storeu_si256((__m256i *)(*dst + 16), cd_shift);
  *dst += stride;

  ALIGNR_256_HI(ab_shift, *cd, *ab, 4);
  ALIGNR_256_HI(cd_shift, *ef, *cd, 4);
  _mm256_storeu_si256((__m256i *)*dst, ab_shift);
  _mm256_storeu_si256((__m256i *)(*dst + 16), cd_shift);
  *dst += stride;

  ALIGNR_256_HI(ab_shift, *cd, *ab, 8);
  ALIGNR_256_HI(cd_shift, *ef, *cd, 8);
  _mm256_storeu_si256((__m256i *)*dst, ab_shift);
  _mm256_storeu_si256((__m256i *)(*dst + 16), cd_shift);
  *dst += stride;

  ALIGNR_256_HI(ab_shift, *cd, *ab, 12);
  ALIGNR_256_HI(cd_shift, *ef, *cd, 12);
  _mm256_storeu_si256((__m256i *)*dst, ab_shift);
  _mm256_storeu_si256((__m256i *)(*dst + 16), cd_shift);
  *dst += stride;
}

void vpx_highbd_d207_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t stride,
                                          const uint16_t *above,
                                          const uint16_t *left, int bd) {
  (void)above;
  (void)bd;

  const __m256i A0 = _mm256_loadu_si256((const __m256i *)left);  // 0..15
  const __m256i A1 =
      _mm256_loadu_si256((const __m256i *)(left + 16));  // 16..31

  const __m256i AR =
      _mm256_permute4x64_epi64(_mm256_shufflehi_epi16(A1, 0xff), 0xff);

  __m256i B0;
  ALIGNR_256(B0, A1, A0, 2);
  __m256i B1;
  ALIGNR_256(B1, AR, A1, 2);
  __m256i C0;
  ALIGNR_256(C0, A1, A0, 4);
  __m256i C1;
  ALIGNR_256(C1, AR, A1, 4);

  __m256i avg2_0 = _mm256_avg_epu16(A0, B0);
  __m256i avg2_1 = _mm256_avg_epu16(A1, B1);
  __m256i avg3_0 = avg3_epu16_avx2(&A0, &B0, &C0);
  __m256i avg3_1 = avg3_epu16_avx2(&A1, &B1, &C1);

  __m256i out_ac = _mm256_unpacklo_epi16(avg2_0, avg3_0);
  __m256i out_bd = _mm256_unpackhi_epi16(avg2_0, avg3_0);
  __m256i out_eg = _mm256_unpacklo_epi16(avg2_1, avg3_1);
  __m256i out_fh = _mm256_unpackhi_epi16(avg2_1, avg3_1);

  __m256i out_ab = _mm256_permute2x128_si256(out_ac, out_bd, 0x20);
  __m256i out_cd = _mm256_permute2x128_si256(out_ac, out_bd, 0x31);
  __m256i out_ef = _mm256_permute2x128_si256(out_eg, out_fh, 0x20);
  __m256i out_gh = _mm256_permute2x128_si256(out_eg, out_fh, 0x31);

  d207_store_8x32_avx2(&dst, stride, &out_ab, &out_cd, &out_ef);
  d207_store_8x32_avx2(&dst, stride, &out_cd, &out_ef, &out_gh);
  d207_store_8x32_avx2(&dst, stride, &out_ef, &out_gh, &AR);
  d207_store_8x32_avx2(&dst, stride, &out_gh, &AR, &AR);
}
