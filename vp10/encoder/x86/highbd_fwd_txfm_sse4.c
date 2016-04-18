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

static INLINE void write_buffer_4x4(tran_low_t *output, __m128i *res) {
  _mm_store_si128((__m128i *)(output + 0 * 4), res[0]);
  _mm_store_si128((__m128i *)(output + 1 * 4), res[1]);
  _mm_store_si128((__m128i *)(output + 2 * 4), res[2]);
  _mm_store_si128((__m128i *)(output + 3 * 4), res[3]);
}

void vp10_highbd_fht4x4_sse4_1(const int16_t *input, tran_low_t *output,
                               int stride, int tx_type) {
  __m128i in[4];
  const TXFM_2D_CFG *cfg;
  int bit;

  switch (tx_type) {
    case DCT_DCT:
      cfg = &fwd_txfm_2d_cfg_dct_dct_4;
      load_buffer_4x4(input, in, stride, 0, 0, cfg->shift[0]);
      bit = cfg->cos_bit_col[2];
      fdct4x4_sse4_1(in, bit);
      bit = cfg->cos_bit_row[2];
      fdct4x4_sse4_1(in, bit);
      write_buffer_4x4(output, in);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_highbd_fht4x4_c(input, output, stride, tx_type);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      vp10_highbd_fht4x4_c(input, output, stride, tx_type);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      vp10_highbd_fht4x4_c(input, output, stride, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
  }
}
