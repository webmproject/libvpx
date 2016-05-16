/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <smmintrin.h> /* SSE4.1 */

#include "./vp10_rtcd.h"
#include "./vpx_config.h"
#include "vp10/common/vp10_inv_txfm2d_cfg.h"


static INLINE void load_buffer_4x4(const int32_t *coeff, __m128i *in) {
  in[0] = _mm_loadu_si128((const __m128i *)(coeff + 0));
  in[1] = _mm_loadu_si128((const __m128i *)(coeff + 4));
  in[2] = _mm_loadu_si128((const __m128i *)(coeff + 8));
  in[3] = _mm_loadu_si128((const __m128i *)(coeff + 12));
}

static void idct4x4_sse4_1(__m128i *in, int bit) {
  const int32_t *cospi = cospi_arr[bit - cos_bit_min];
  const __m128i cospi32 = _mm_set1_epi32(cospi[32]);
  const __m128i cospi48 = _mm_set1_epi32(cospi[48]);
  const __m128i cospi16 = _mm_set1_epi32(cospi[16]);
  const __m128i cospim16 = _mm_set1_epi32(-cospi[16]);
  const __m128i rnding = _mm_set1_epi32(1 << (bit - 1));
  __m128i u0, u1, u2, u3;
  __m128i v0, v1, v2, v3, x, y;

  v0 = _mm_unpacklo_epi32(in[0], in[1]);
  v1 = _mm_unpackhi_epi32(in[0], in[1]);
  v2 = _mm_unpacklo_epi32(in[2], in[3]);
  v3 = _mm_unpackhi_epi32(in[2], in[3]);

  u0 = _mm_unpacklo_epi64(v0, v2);
  u1 = _mm_unpackhi_epi64(v0, v2);
  u2 = _mm_unpacklo_epi64(v1, v3);
  u3 = _mm_unpackhi_epi64(v1, v3);

  x = _mm_mullo_epi32(u0, cospi32);
  y = _mm_mullo_epi32(u2, cospi32);
  v0 = _mm_add_epi32(x, y);
  v0 = _mm_add_epi32(v0, rnding);
  v0 = _mm_srai_epi32(v0, bit);

  v1 = _mm_sub_epi32(x, y);
  v1 = _mm_add_epi32(v1, rnding);
  v1 = _mm_srai_epi32(v1, bit);

  x = _mm_mullo_epi32(u1, cospi48);
  y = _mm_mullo_epi32(u3, cospim16);
  v2 = _mm_add_epi32(x, y);
  v2 = _mm_add_epi32(v2, rnding);
  v2 = _mm_srai_epi32(v2, bit);

  x = _mm_mullo_epi32(u1, cospi16);
  y = _mm_mullo_epi32(u3, cospi48);
  v3 = _mm_add_epi32(x, y);
  v3 = _mm_add_epi32(v3, rnding);
  v3 = _mm_srai_epi32(v3, bit);

  in[0] = _mm_add_epi32(v0, v3);
  in[1] = _mm_add_epi32(v1, v2);
  in[2] = _mm_sub_epi32(v1, v2);
  in[3] = _mm_sub_epi32(v0, v3);
}

static void iadst4x4_sse4_1(__m128i *in, int bit) {
  const int32_t *cospi = cospi_arr[bit - cos_bit_min];
  const __m128i cospi32 = _mm_set1_epi32(cospi[32]);
  const __m128i cospi8 = _mm_set1_epi32(cospi[8]);
  const __m128i cospim8 = _mm_set1_epi32(-cospi[8]);
  const __m128i cospi40 = _mm_set1_epi32(cospi[40]);
  const __m128i cospim40 = _mm_set1_epi32(-cospi[40]);
  const __m128i cospi24 = _mm_set1_epi32(cospi[24]);
  const __m128i cospi56 = _mm_set1_epi32(cospi[56]);
  const __m128i rnding = _mm_set1_epi32(1 << (bit - 1));
  const __m128i zero = _mm_setzero_si128();
  __m128i u0, u1, u2, u3;
  __m128i v0, v1, v2, v3, x, y;

  v0 = _mm_unpacklo_epi32(in[0], in[1]);
  v1 = _mm_unpackhi_epi32(in[0], in[1]);
  v2 = _mm_unpacklo_epi32(in[2], in[3]);
  v3 = _mm_unpackhi_epi32(in[2], in[3]);

  u0 = _mm_unpacklo_epi64(v0, v2);
  u1 = _mm_unpackhi_epi64(v0, v2);
  u2 = _mm_unpacklo_epi64(v1, v3);
  u3 = _mm_unpackhi_epi64(v1, v3);

  // stage 0
  // stage 1
  u1 = _mm_sub_epi32(zero, u1);
  u3 = _mm_sub_epi32(zero, u3);

  // stage 2
  v0 = u0;
  v1 = u3;
  x = _mm_mullo_epi32(u1, cospi32);
  y = _mm_mullo_epi32(u2, cospi32);
  v2 = _mm_add_epi32(x, y);
  v2 = _mm_add_epi32(v2, rnding);
  v2 = _mm_srai_epi32(v2, bit);

  v3 = _mm_sub_epi32(x, y);
  v3 = _mm_add_epi32(v3, rnding);
  v3 = _mm_srai_epi32(v3, bit);

  // stage 3
  u0 = _mm_add_epi32(v0, v2);
  u1 = _mm_add_epi32(v1, v3);
  u2 = _mm_sub_epi32(v0, v2);
  u3 = _mm_sub_epi32(v1, v3);

  // stage 4
  x = _mm_mullo_epi32(u0, cospi8);
  y = _mm_mullo_epi32(u1, cospi56);
  in[3] = _mm_add_epi32(x, y);
  in[3] = _mm_add_epi32(in[3], rnding);
  in[3] = _mm_srai_epi32(in[3], bit);

  x = _mm_mullo_epi32(u0, cospi56);
  y = _mm_mullo_epi32(u1, cospim8);
  in[0] = _mm_add_epi32(x, y);
  in[0] = _mm_add_epi32(in[0], rnding);
  in[0] = _mm_srai_epi32(in[0], bit);

  x = _mm_mullo_epi32(u2, cospi40);
  y = _mm_mullo_epi32(u3, cospi24);
  in[1] = _mm_add_epi32(x, y);
  in[1] = _mm_add_epi32(in[1], rnding);
  in[1] = _mm_srai_epi32(in[1], bit);

  x = _mm_mullo_epi32(u2, cospi24);
  y = _mm_mullo_epi32(u3, cospim40);
  in[2] = _mm_add_epi32(x, y);
  in[2] = _mm_add_epi32(in[2], rnding);
  in[2] = _mm_srai_epi32(in[2], bit);
}

static INLINE void round_shift_4x4(__m128i *in, int shift) {
  __m128i rnding = _mm_set1_epi32(1 << (shift - 1));

  in[0] = _mm_add_epi32(in[0], rnding);
  in[1] = _mm_add_epi32(in[1], rnding);
  in[2] = _mm_add_epi32(in[2], rnding);
  in[3] = _mm_add_epi32(in[3], rnding);

  in[0] = _mm_srai_epi32(in[0], shift);
  in[1] = _mm_srai_epi32(in[1], shift);
  in[2] = _mm_srai_epi32(in[2], shift);
  in[3] = _mm_srai_epi32(in[3], shift);
}

static INLINE __m128i highbd_clamp_epi16(__m128i u, int bd) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i one = _mm_set1_epi16(1);
  const __m128i max = _mm_sub_epi16(_mm_slli_epi16(one, bd), one);
  __m128i clamped, mask;

  mask = _mm_cmpgt_epi16(u, max);
  clamped = _mm_andnot_si128(mask, u);
  mask = _mm_and_si128(mask, max);
  clamped = _mm_or_si128(mask, clamped);
  mask = _mm_cmpgt_epi16(clamped, zero);
  clamped = _mm_and_si128(clamped, mask);

  return clamped;
}

static void write_buffer_4x4(__m128i *in, uint16_t *output, int stride,
                             int flipud, int fliplr, int shift, int bd) {
  const __m128i zero = _mm_setzero_si128();
  __m128i u0, u1, u2, u3;
  __m128i v0, v1, v2, v3;

  round_shift_4x4(in, shift);

  v0 = _mm_loadl_epi64((__m128i const *)(output + 0 * stride));
  v1 = _mm_loadl_epi64((__m128i const *)(output + 1 * stride));
  v2 = _mm_loadl_epi64((__m128i const *)(output + 2 * stride));
  v3 = _mm_loadl_epi64((__m128i const *)(output + 3 * stride));

  v0 = _mm_unpacklo_epi16(v0, zero);
  v1 = _mm_unpacklo_epi16(v1, zero);
  v2 = _mm_unpacklo_epi16(v2, zero);
  v3 = _mm_unpacklo_epi16(v3, zero);

  u0 = _mm_add_epi32(in[0], v0);
  u1 = _mm_add_epi32(in[1], v1);
  u2 = _mm_add_epi32(in[2], v2);
  u3 = _mm_add_epi32(in[3], v3);

  v0 = _mm_packus_epi32(u0, u1);
  v2 = _mm_packus_epi32(u2, u3);

  u0 = highbd_clamp_epi16(v0, bd);
  u2 = highbd_clamp_epi16(v2, bd);

  v0 = _mm_unpacklo_epi64(u0, u0);
  v1 = _mm_unpackhi_epi64(u0, u0);
  v2 = _mm_unpacklo_epi64(u2, u2);
  v3 = _mm_unpackhi_epi64(u2, u2);

  _mm_storel_epi64((__m128i *)(output + 0 * stride), v0);
  _mm_storel_epi64((__m128i *)(output + 1 * stride), v1);
  _mm_storel_epi64((__m128i *)(output + 2 * stride), v2);
  _mm_storel_epi64((__m128i *)(output + 3 * stride), v3);

  (void) flipud;
  (void) fliplr;
}

void vp10_inv_txfm2d_add_4x4_sse4_1(const int32_t *coeff, uint16_t *output,
                                    int stride, int tx_type, int bd) {
  __m128i in[4];
  const TXFM_2D_CFG *cfg = NULL;

  switch (tx_type) {
    case DCT_DCT:
      cfg = &inv_txfm_2d_cfg_dct_dct_4;
      load_buffer_4x4(coeff, in);
      idct4x4_sse4_1(in, cfg->cos_bit_row[2]);
      idct4x4_sse4_1(in, cfg->cos_bit_row[2]);
      write_buffer_4x4(in, output, stride, 0, 0, -cfg->shift[1], bd);
      break;
    case ADST_DCT:
      cfg = &inv_txfm_2d_cfg_adst_dct_4;
      load_buffer_4x4(coeff, in);
      idct4x4_sse4_1(in, cfg->cos_bit_row[2]);
      iadst4x4_sse4_1(in, cfg->cos_bit_row[2]);
      write_buffer_4x4(in, output, stride, 0, 0, -cfg->shift[1], bd);
      break;
    case DCT_ADST:
      cfg = &inv_txfm_2d_cfg_dct_adst_4;
      load_buffer_4x4(coeff, in);
      iadst4x4_sse4_1(in, cfg->cos_bit_row[2]);
      idct4x4_sse4_1(in, cfg->cos_bit_row[2]);
      write_buffer_4x4(in, output, stride, 0, 0, -cfg->shift[1], bd);
      break;
    case ADST_ADST:
      cfg = &inv_txfm_2d_cfg_adst_adst_4;
      load_buffer_4x4(coeff, in);
      iadst4x4_sse4_1(in, cfg->cos_bit_row[2]);
      iadst4x4_sse4_1(in, cfg->cos_bit_row[2]);
      write_buffer_4x4(in, output, stride, 0, 0, -cfg->shift[1], bd);
      break;
    default:
      assert(0);
  }
}
