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
#include "vp10/common/vp10_fwd_txfm2d_cfg.h"
#include "vp10/common/vp10_txfm.h"
#include "vpx_dsp/txfm_common.h"
#include "vpx_ports/mem.h"

static INLINE void load_buffer_4x4(const int16_t *input, __m128i *in,
                                   int stride, int flipud, int fliplr,
                                   int shift) {
  if (!flipud) {
    in[0] = _mm_loadl_epi64((const __m128i *)(input + 0 * stride));
    in[1] = _mm_loadl_epi64((const __m128i *)(input + 1 * stride));
    in[2] = _mm_loadl_epi64((const __m128i *)(input + 2 * stride));
    in[3] = _mm_loadl_epi64((const __m128i *)(input + 3 * stride));
  } else {
    in[0] = _mm_loadl_epi64((const __m128i *)(input + 3 * stride));
    in[1] = _mm_loadl_epi64((const __m128i *)(input + 2 * stride));
    in[2] = _mm_loadl_epi64((const __m128i *)(input + 1 * stride));
    in[3] = _mm_loadl_epi64((const __m128i *)(input + 0 * stride));
  }

  if (fliplr) {
    in[0] = _mm_shufflelo_epi16(in[0], 0x1b);
    in[1] = _mm_shufflelo_epi16(in[1], 0x1b);
    in[2] = _mm_shufflelo_epi16(in[2], 0x1b);
    in[3] = _mm_shufflelo_epi16(in[3], 0x1b);
  }

  in[0] = _mm_cvtepi16_epi32(in[0]);
  in[1] = _mm_cvtepi16_epi32(in[1]);
  in[2] = _mm_cvtepi16_epi32(in[2]);
  in[3] = _mm_cvtepi16_epi32(in[3]);

  in[0] = _mm_slli_epi32(in[0], shift);
  in[1] = _mm_slli_epi32(in[1], shift);
  in[2] = _mm_slli_epi32(in[2], shift);
  in[3] = _mm_slli_epi32(in[3], shift);
}

// We only use stage-2 bit;
// shift[0] is used in load_buffer_4x4()
// shift[1] is used in txfm_func_col()
// shift[2] is used in txfm_func_row()
static void fdct4x4_sse4_1(__m128i *in, int bit) {
  const int32_t *cospi = cospi_arr[bit - cos_bit_min];
  const __m128i cospi32 = _mm_set1_epi32(cospi[32]);
  const __m128i cospi48 = _mm_set1_epi32(cospi[48]);
  const __m128i cospi16 = _mm_set1_epi32(cospi[16]);
  const __m128i rnding = _mm_set1_epi32(1 << (bit - 1));
  __m128i s0, s1, s2, s3;
  __m128i u0, u1, u2, u3;
  __m128i v0, v1, v2, v3;

  s0 = _mm_add_epi32(in[0], in[3]);
  s1 = _mm_add_epi32(in[1], in[2]);
  s2 = _mm_sub_epi32(in[1], in[2]);
  s3 = _mm_sub_epi32(in[0], in[3]);

  // btf_32_sse4_1_type0(cospi32, cospi32, s[01], u[02], bit);
  u0 = _mm_mullo_epi32(s0, cospi32);
  u1 = _mm_mullo_epi32(s1, cospi32);
  u2 = _mm_add_epi32(u0, u1);
  v0 = _mm_sub_epi32(u0, u1);

  u3 = _mm_add_epi32(u2, rnding);
  v1 = _mm_add_epi32(v0, rnding);

  u0 = _mm_srai_epi32(u3, bit);
  u2 = _mm_srai_epi32(v1, bit);

  // btf_32_sse4_1_type1(cospi48, cospi16, s[23], u[13], bit);
  v0 = _mm_mullo_epi32(s2, cospi48);
  v1 = _mm_mullo_epi32(s3, cospi16);
  v2 = _mm_add_epi32(v0, v1);

  v3 = _mm_add_epi32(v2, rnding);
  u1 = _mm_srai_epi32(v3, bit);

  v0 = _mm_mullo_epi32(s2, cospi16);
  v1 = _mm_mullo_epi32(s3, cospi48);
  v2 = _mm_sub_epi32(v1, v0);

  v3 = _mm_add_epi32(v2, rnding);
  u3 = _mm_srai_epi32(v3, bit);

  // Note: shift[1] and shift[2] are zeros

  // Transpose 4x4 32-bit
  v0 = _mm_unpacklo_epi32(u0, u1);
  v1 = _mm_unpackhi_epi32(u0, u1);
  v2 = _mm_unpacklo_epi32(u2, u3);
  v3 = _mm_unpackhi_epi32(u2, u3);

  in[0] = _mm_unpacklo_epi64(v0, v2);
  in[1] = _mm_unpackhi_epi64(v0, v2);
  in[2] = _mm_unpacklo_epi64(v1, v3);
  in[3] = _mm_unpackhi_epi64(v1, v3);
}

static INLINE void write_buffer_4x4(__m128i *res, tran_low_t *output) {
  _mm_store_si128((__m128i *)(output + 0 * 4), res[0]);
  _mm_store_si128((__m128i *)(output + 1 * 4), res[1]);
  _mm_store_si128((__m128i *)(output + 2 * 4), res[2]);
  _mm_store_si128((__m128i *)(output + 3 * 4), res[3]);
}

// Note:
//  We implement vp10_fwd_txfm2d_4x4(). This function is kept here since
//  vp10_highbd_fht4x4_c() is not removed yet
void vp10_highbd_fht4x4_sse4_1(const int16_t *input, tran_low_t *output,
                               int stride, int tx_type) {
  (void)input;
  (void)output;
  (void)stride;
  (void)tx_type;
  assert(0);
}

static void fadst4x4_sse4_1(__m128i *in, int bit) {
  const int32_t *cospi = cospi_arr[bit - cos_bit_min];
  const __m128i cospi8 = _mm_set1_epi32(cospi[8]);
  const __m128i cospi56 = _mm_set1_epi32(cospi[56]);
  const __m128i cospi40 = _mm_set1_epi32(cospi[40]);
  const __m128i cospi24 = _mm_set1_epi32(cospi[24]);
  const __m128i cospi32 = _mm_set1_epi32(cospi[32]);
  const __m128i rnding = _mm_set1_epi32(1 << (bit - 1));
  const __m128i kZero = _mm_setzero_si128();
  __m128i s0, s1, s2, s3;
  __m128i u0, u1, u2, u3;
  __m128i v0, v1, v2, v3;

  // stage 0
  // stage 1
  // stage 2
  u0 = _mm_mullo_epi32(in[3], cospi8);
  u1 = _mm_mullo_epi32(in[0], cospi56);
  u2 = _mm_add_epi32(u0, u1);
  s0 = _mm_add_epi32(u2, rnding);
  s0 = _mm_srai_epi32(s0, bit);

  v0 = _mm_mullo_epi32(in[3], cospi56);
  v1 = _mm_mullo_epi32(in[0], cospi8);
  v2 = _mm_sub_epi32(v0, v1);
  s1 = _mm_add_epi32(v2, rnding);
  s1 = _mm_srai_epi32(s1, bit);

  u0 = _mm_mullo_epi32(in[1], cospi40);
  u1 = _mm_mullo_epi32(in[2], cospi24);
  u2 = _mm_add_epi32(u0, u1);
  s2 = _mm_add_epi32(u2, rnding);
  s2 = _mm_srai_epi32(s2, bit);

  v0 = _mm_mullo_epi32(in[1], cospi24);
  v1 = _mm_mullo_epi32(in[2], cospi40);
  v2 = _mm_sub_epi32(v0, v1);
  s3 = _mm_add_epi32(v2, rnding);
  s3 = _mm_srai_epi32(s3, bit);

  // stage 3
  u0 = _mm_add_epi32(s0, s2);
  u2 = _mm_sub_epi32(s0, s2);
  u1 = _mm_add_epi32(s1, s3);
  u3 = _mm_sub_epi32(s1, s3);

  // stage 4
  v0 = _mm_mullo_epi32(u2, cospi32);
  v1 = _mm_mullo_epi32(u3, cospi32);
  v2 = _mm_add_epi32(v0, v1);
  s2 = _mm_add_epi32(v2, rnding);
  u2 = _mm_srai_epi32(s2, bit);

  v2 = _mm_sub_epi32(v0, v1);
  s3 = _mm_add_epi32(v2, rnding);
  u3 = _mm_srai_epi32(s3, bit);

  // u0, u1, u2, u3
  u2 = _mm_sub_epi32(kZero, u2);
  u1 = _mm_sub_epi32(kZero, u1);

  // u0, u2, u3, u1
  // Transpose 4x4 32-bit
  v0 = _mm_unpacklo_epi32(u0, u2);
  v1 = _mm_unpackhi_epi32(u0, u2);
  v2 = _mm_unpacklo_epi32(u3, u1);
  v3 = _mm_unpackhi_epi32(u3, u1);

  in[0] = _mm_unpacklo_epi64(v0, v2);
  in[1] = _mm_unpackhi_epi64(v0, v2);
  in[2] = _mm_unpacklo_epi64(v1, v3);
  in[3] = _mm_unpackhi_epi64(v1, v3);
}

void vp10_fwd_txfm2d_4x4_sse4_1(const int16_t *input, tran_low_t *coeff,
                                int input_stride, int tx_type,
                                const int bd) {
  __m128i in[4];
  const TXFM_2D_CFG *cfg = NULL;

  switch (tx_type) {
    case DCT_DCT:
      cfg = &fwd_txfm_2d_cfg_dct_dct_4;
      load_buffer_4x4(input, in, input_stride, 0, 0, cfg->shift[0]);
      fdct4x4_sse4_1(in, cfg->cos_bit_col[2]);
      fdct4x4_sse4_1(in, cfg->cos_bit_row[2]);
      write_buffer_4x4(in, coeff);
      break;
    case ADST_DCT:
      cfg = &fwd_txfm_2d_cfg_adst_dct_4;
      load_buffer_4x4(input, in, input_stride, 0, 0, cfg->shift[0]);
      fadst4x4_sse4_1(in, cfg->cos_bit_col[2]);
      fdct4x4_sse4_1(in, cfg->cos_bit_row[2]);
      write_buffer_4x4(in, coeff);
      break;
    case DCT_ADST:
      cfg = &fwd_txfm_2d_cfg_dct_adst_4;
      load_buffer_4x4(input, in, input_stride, 0, 0, cfg->shift[0]);
      fdct4x4_sse4_1(in, cfg->cos_bit_col[2]);
      fadst4x4_sse4_1(in, cfg->cos_bit_row[2]);
      write_buffer_4x4(in, coeff);
      break;
    case ADST_ADST:
      cfg = &fwd_txfm_2d_cfg_adst_adst_4;
      load_buffer_4x4(input, in, input_stride, 0, 0, cfg->shift[0]);
      fadst4x4_sse4_1(in, cfg->cos_bit_col[2]);
      fadst4x4_sse4_1(in, cfg->cos_bit_row[2]);
      write_buffer_4x4(in, coeff);
      break;
    default:
      assert(0);
  }
  (void)bd;
}
