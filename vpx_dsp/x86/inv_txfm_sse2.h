/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_DSP_X86_INV_TXFM_SSE2_H_
#define VPX_DSP_X86_INV_TXFM_SSE2_H_

#include <emmintrin.h>  // SSE2

#include "./vpx_config.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/inv_txfm.h"
#include "vpx_dsp/x86/transpose_sse2.h"
#include "vpx_dsp/x86/txfm_common_sse2.h"

static INLINE void idct8x8_12_transpose_16bit_4x8(const __m128i *const in,
                                                  __m128i *const out) {
  // Unpack 16 bit elements. Goes from:
  // in[0]: 30 31 32 33  00 01 02 03
  // in[1]: 20 21 22 23  10 11 12 13
  // in[2]: 40 41 42 43  70 71 72 73
  // in[3]: 50 51 52 53  60 61 62 63
  // to:
  // tr0_0: 00 10 01 11  02 12 03 13
  // tr0_1: 20 30 21 31  22 32 23 33
  // tr0_2: 40 50 41 51  42 52 43 53
  // tr0_3: 60 70 61 71  62 72 63 73
  const __m128i tr0_0 = _mm_unpackhi_epi16(in[0], in[1]);
  const __m128i tr0_1 = _mm_unpacklo_epi16(in[1], in[0]);
  const __m128i tr0_2 = _mm_unpacklo_epi16(in[2], in[3]);
  const __m128i tr0_3 = _mm_unpackhi_epi16(in[3], in[2]);

  // Unpack 32 bit elements resulting in:
  // tr1_0: 00 10 20 30  01 11 21 31
  // tr1_1: 02 12 22 32  03 13 23 33
  // tr1_2: 40 50 60 70  41 51 61 71
  // tr1_3: 42 52 62 72  43 53 63 73
  const __m128i tr1_0 = _mm_unpacklo_epi32(tr0_0, tr0_1);
  const __m128i tr1_1 = _mm_unpacklo_epi32(tr0_2, tr0_3);
  const __m128i tr1_2 = _mm_unpackhi_epi32(tr0_0, tr0_1);
  const __m128i tr1_3 = _mm_unpackhi_epi32(tr0_2, tr0_3);

  // Unpack 64 bit elements resulting in:
  // out[0]: 00 10 20 30  40 50 60 70
  // out[1]: 01 11 21 31  41 51 61 71
  // out[2]: 02 12 22 32  42 52 62 72
  // out[3]: 03 13 23 33  43 53 63 73
  out[0] = _mm_unpacklo_epi64(tr1_0, tr1_1);
  out[1] = _mm_unpackhi_epi64(tr1_0, tr1_1);
  out[2] = _mm_unpacklo_epi64(tr1_2, tr1_3);
  out[3] = _mm_unpackhi_epi64(tr1_2, tr1_3);
}

static INLINE __m128i dct_const_round_shift_sse2(const __m128i in) {
  const __m128i t = _mm_add_epi32(in, _mm_set1_epi32(DCT_CONST_ROUNDING));
  return _mm_srai_epi32(t, DCT_CONST_BITS);
}

static INLINE __m128i idct_madd_round_shift_sse2(const __m128i in,
                                                 const __m128i cospi) {
  const __m128i t = _mm_madd_epi16(in, cospi);
  return dct_const_round_shift_sse2(t);
}

// Calculate the dot product between in0/1 and x and wrap to short.
static INLINE __m128i idct_calc_wraplow_sse2(const __m128i in0,
                                             const __m128i in1,
                                             const __m128i x) {
  const __m128i t0 = idct_madd_round_shift_sse2(in0, x);
  const __m128i t1 = idct_madd_round_shift_sse2(in1, x);
  return _mm_packs_epi32(t0, t1);
}

// Multiply elements by constants and add them together.
static INLINE void butterfly(const __m128i in0, const __m128i in1, const int c0,
                             const int c1, __m128i *const res0,
                             __m128i *const res1) {
  const __m128i cst0 = pair_set_epi16(c0, -c1);
  const __m128i cst1 = pair_set_epi16(c1, c0);
  const __m128i lo = _mm_unpacklo_epi16(in0, in1);
  const __m128i hi = _mm_unpackhi_epi16(in0, in1);
  *res0 = idct_calc_wraplow_sse2(lo, hi, cst0);
  *res1 = idct_calc_wraplow_sse2(lo, hi, cst1);
}

// Functions to allow 8 bit optimisations to be used when profile 0 is used with
// highbitdepth enabled
static INLINE __m128i load_input_data4(const tran_low_t *data) {
#if CONFIG_VP9_HIGHBITDEPTH
  const __m128i zero = _mm_setzero_si128();
  const __m128i in = _mm_load_si128((const __m128i *)data);
  return _mm_packs_epi32(in, zero);
#else
  return _mm_loadl_epi64((const __m128i *)data);
#endif
}

static INLINE __m128i load_input_data8(const tran_low_t *data) {
#if CONFIG_VP9_HIGHBITDEPTH
  const __m128i in0 = _mm_load_si128((const __m128i *)data);
  const __m128i in1 = _mm_load_si128((const __m128i *)(data + 4));
  return _mm_packs_epi32(in0, in1);
#else
  return _mm_load_si128((const __m128i *)data);
#endif
}

static INLINE void load_buffer_8x8(const tran_low_t *const input,
                                   __m128i *const in) {
  in[0] = load_input_data8(input + 0 * 8);
  in[1] = load_input_data8(input + 1 * 8);
  in[2] = load_input_data8(input + 2 * 8);
  in[3] = load_input_data8(input + 3 * 8);
  in[4] = load_input_data8(input + 4 * 8);
  in[5] = load_input_data8(input + 5 * 8);
  in[6] = load_input_data8(input + 6 * 8);
  in[7] = load_input_data8(input + 7 * 8);
}

static INLINE void load_buffer_8x16(const tran_low_t *const input,
                                    __m128i *const in) {
  in[0] = load_input_data8(input + 0 * 16);
  in[1] = load_input_data8(input + 1 * 16);
  in[2] = load_input_data8(input + 2 * 16);
  in[3] = load_input_data8(input + 3 * 16);
  in[4] = load_input_data8(input + 4 * 16);
  in[5] = load_input_data8(input + 5 * 16);
  in[6] = load_input_data8(input + 6 * 16);
  in[7] = load_input_data8(input + 7 * 16);

  in[8] = load_input_data8(input + 8 * 16);
  in[9] = load_input_data8(input + 9 * 16);
  in[10] = load_input_data8(input + 10 * 16);
  in[11] = load_input_data8(input + 11 * 16);
  in[12] = load_input_data8(input + 12 * 16);
  in[13] = load_input_data8(input + 13 * 16);
  in[14] = load_input_data8(input + 14 * 16);
  in[15] = load_input_data8(input + 15 * 16);
}

static INLINE void recon_and_store(uint8_t *const dest, const __m128i in_x) {
  const __m128i zero = _mm_setzero_si128();
  __m128i d0 = _mm_loadl_epi64((__m128i *)(dest));
  d0 = _mm_unpacklo_epi8(d0, zero);
  d0 = _mm_add_epi16(in_x, d0);
  d0 = _mm_packus_epi16(d0, d0);
  _mm_storel_epi64((__m128i *)(dest), d0);
}

static INLINE void round_shift_8x8(const __m128i *const in,
                                   __m128i *const out) {
  const __m128i final_rounding = _mm_set1_epi16(1 << 4);

  out[0] = _mm_add_epi16(in[0], final_rounding);
  out[1] = _mm_add_epi16(in[1], final_rounding);
  out[2] = _mm_add_epi16(in[2], final_rounding);
  out[3] = _mm_add_epi16(in[3], final_rounding);
  out[4] = _mm_add_epi16(in[4], final_rounding);
  out[5] = _mm_add_epi16(in[5], final_rounding);
  out[6] = _mm_add_epi16(in[6], final_rounding);
  out[7] = _mm_add_epi16(in[7], final_rounding);

  out[0] = _mm_srai_epi16(out[0], 5);
  out[1] = _mm_srai_epi16(out[1], 5);
  out[2] = _mm_srai_epi16(out[2], 5);
  out[3] = _mm_srai_epi16(out[3], 5);
  out[4] = _mm_srai_epi16(out[4], 5);
  out[5] = _mm_srai_epi16(out[5], 5);
  out[6] = _mm_srai_epi16(out[6], 5);
  out[7] = _mm_srai_epi16(out[7], 5);
}

static INLINE void write_buffer_8x8(const __m128i *const in,
                                    uint8_t *const dest, const int stride) {
  __m128i t[8];

  round_shift_8x8(in, t);

  recon_and_store(dest + 0 * stride, t[0]);
  recon_and_store(dest + 1 * stride, t[1]);
  recon_and_store(dest + 2 * stride, t[2]);
  recon_and_store(dest + 3 * stride, t[3]);
  recon_and_store(dest + 4 * stride, t[4]);
  recon_and_store(dest + 5 * stride, t[5]);
  recon_and_store(dest + 6 * stride, t[6]);
  recon_and_store(dest + 7 * stride, t[7]);
}

static INLINE void write_buffer_8x16(uint8_t *const dest, __m128i *const in,
                                     const int stride) {
  const __m128i final_rounding = _mm_set1_epi16(1 << 5);
  // Final rounding and shift
  in[0] = _mm_adds_epi16(in[0], final_rounding);
  in[1] = _mm_adds_epi16(in[1], final_rounding);
  in[2] = _mm_adds_epi16(in[2], final_rounding);
  in[3] = _mm_adds_epi16(in[3], final_rounding);
  in[4] = _mm_adds_epi16(in[4], final_rounding);
  in[5] = _mm_adds_epi16(in[5], final_rounding);
  in[6] = _mm_adds_epi16(in[6], final_rounding);
  in[7] = _mm_adds_epi16(in[7], final_rounding);
  in[8] = _mm_adds_epi16(in[8], final_rounding);
  in[9] = _mm_adds_epi16(in[9], final_rounding);
  in[10] = _mm_adds_epi16(in[10], final_rounding);
  in[11] = _mm_adds_epi16(in[11], final_rounding);
  in[12] = _mm_adds_epi16(in[12], final_rounding);
  in[13] = _mm_adds_epi16(in[13], final_rounding);
  in[14] = _mm_adds_epi16(in[14], final_rounding);
  in[15] = _mm_adds_epi16(in[15], final_rounding);

  in[0] = _mm_srai_epi16(in[0], 6);
  in[1] = _mm_srai_epi16(in[1], 6);
  in[2] = _mm_srai_epi16(in[2], 6);
  in[3] = _mm_srai_epi16(in[3], 6);
  in[4] = _mm_srai_epi16(in[4], 6);
  in[5] = _mm_srai_epi16(in[5], 6);
  in[6] = _mm_srai_epi16(in[6], 6);
  in[7] = _mm_srai_epi16(in[7], 6);
  in[8] = _mm_srai_epi16(in[8], 6);
  in[9] = _mm_srai_epi16(in[9], 6);
  in[10] = _mm_srai_epi16(in[10], 6);
  in[11] = _mm_srai_epi16(in[11], 6);
  in[12] = _mm_srai_epi16(in[12], 6);
  in[13] = _mm_srai_epi16(in[13], 6);
  in[14] = _mm_srai_epi16(in[14], 6);
  in[15] = _mm_srai_epi16(in[15], 6);

  recon_and_store(dest + 0 * stride, in[0]);
  recon_and_store(dest + 1 * stride, in[1]);
  recon_and_store(dest + 2 * stride, in[2]);
  recon_and_store(dest + 3 * stride, in[3]);
  recon_and_store(dest + 4 * stride, in[4]);
  recon_and_store(dest + 5 * stride, in[5]);
  recon_and_store(dest + 6 * stride, in[6]);
  recon_and_store(dest + 7 * stride, in[7]);
  recon_and_store(dest + 8 * stride, in[8]);
  recon_and_store(dest + 9 * stride, in[9]);
  recon_and_store(dest + 10 * stride, in[10]);
  recon_and_store(dest + 11 * stride, in[11]);
  recon_and_store(dest + 12 * stride, in[12]);
  recon_and_store(dest + 13 * stride, in[13]);
  recon_and_store(dest + 14 * stride, in[14]);
  recon_and_store(dest + 15 * stride, in[15]);
}

static INLINE void recon_and_store4x4_sse2(const __m128i *const in,
                                           uint8_t *const dest,
                                           const int stride) {
  const __m128i zero = _mm_setzero_si128();
  __m128i d[2];

  // Reconstruction and Store
  d[0] = _mm_cvtsi32_si128(*(const int *)(dest));
  d[1] = _mm_cvtsi32_si128(*(const int *)(dest + stride * 3));
  d[0] = _mm_unpacklo_epi32(d[0],
                            _mm_cvtsi32_si128(*(const int *)(dest + stride)));
  d[1] = _mm_unpacklo_epi32(
      _mm_cvtsi32_si128(*(const int *)(dest + stride * 2)), d[1]);
  d[0] = _mm_unpacklo_epi8(d[0], zero);
  d[1] = _mm_unpacklo_epi8(d[1], zero);
  d[0] = _mm_add_epi16(d[0], in[0]);
  d[1] = _mm_add_epi16(d[1], in[1]);
  d[0] = _mm_packus_epi16(d[0], d[1]);

  *(int *)dest = _mm_cvtsi128_si32(d[0]);
  d[0] = _mm_srli_si128(d[0], 4);
  *(int *)(dest + stride) = _mm_cvtsi128_si32(d[0]);
  d[0] = _mm_srli_si128(d[0], 4);
  *(int *)(dest + stride * 2) = _mm_cvtsi128_si32(d[0]);
  d[0] = _mm_srli_si128(d[0], 4);
  *(int *)(dest + stride * 3) = _mm_cvtsi128_si32(d[0]);
}

static INLINE void store_buffer_8x32(__m128i *in, uint8_t *dst, int stride) {
  const __m128i final_rounding = _mm_set1_epi16(1 << 5);
  int j = 0;
  while (j < 32) {
    in[j] = _mm_adds_epi16(in[j], final_rounding);
    in[j + 1] = _mm_adds_epi16(in[j + 1], final_rounding);

    in[j] = _mm_srai_epi16(in[j], 6);
    in[j + 1] = _mm_srai_epi16(in[j + 1], 6);

    recon_and_store(dst, in[j]);
    dst += stride;
    recon_and_store(dst, in[j + 1]);
    dst += stride;
    j += 2;
  }
}

// Only do addition and subtraction butterfly, size = 16, 32
static INLINE void add_sub_butterfly(const __m128i *in, __m128i *out,
                                     int size) {
  int i = 0;
  const int num = size >> 1;
  const int bound = size - 1;
  while (i < num) {
    out[i] = _mm_add_epi16(in[i], in[bound - i]);
    out[bound - i] = _mm_sub_epi16(in[i], in[bound - i]);
    i++;
  }
}

static INLINE void idct8(const __m128i *const in /*in[8]*/,
                         __m128i *const out /*out[8]*/) {
  __m128i step1[8], step2[8];

  // stage 1
  butterfly(in[1], in[7], (int)cospi_28_64, (int)cospi_4_64, &step1[4],
            &step1[7]);
  butterfly(in[5], in[3], (int)cospi_12_64, (int)cospi_20_64, &step1[5],
            &step1[6]);

  // stage 2
  butterfly(in[0], in[4], (int)cospi_16_64, (int)cospi_16_64, &step2[1],
            &step2[0]);
  butterfly(in[2], in[6], (int)cospi_24_64, (int)cospi_8_64, &step2[2],
            &step2[3]);

  step2[4] = _mm_add_epi16(step1[4], step1[5]);
  step2[5] = _mm_sub_epi16(step1[4], step1[5]);
  step2[6] = _mm_sub_epi16(step1[7], step1[6]);
  step2[7] = _mm_add_epi16(step1[7], step1[6]);

  // stage 3
  step1[0] = _mm_add_epi16(step2[0], step2[3]);
  step1[1] = _mm_add_epi16(step2[1], step2[2]);
  step1[2] = _mm_sub_epi16(step2[1], step2[2]);
  step1[3] = _mm_sub_epi16(step2[0], step2[3]);
  butterfly(step2[6], step2[5], (int)cospi_16_64, (int)cospi_16_64, &step1[5],
            &step1[6]);

  // stage 4
  out[0] = _mm_add_epi16(step1[0], step2[7]);
  out[1] = _mm_add_epi16(step1[1], step1[6]);
  out[2] = _mm_add_epi16(step1[2], step1[5]);
  out[3] = _mm_add_epi16(step1[3], step2[4]);
  out[4] = _mm_sub_epi16(step1[3], step2[4]);
  out[5] = _mm_sub_epi16(step1[2], step1[5]);
  out[6] = _mm_sub_epi16(step1[1], step1[6]);
  out[7] = _mm_sub_epi16(step1[0], step2[7]);
}

static INLINE void idct8x8_12_add_kernel_sse2(__m128i *const io /*io[8]*/) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i cp_16_16 = pair_set_epi16(cospi_16_64, cospi_16_64);
  const __m128i cp_16_n16 = pair_set_epi16(cospi_16_64, -cospi_16_64);
  __m128i step1[8], step2[8], tmp[4];

  transpose_16bit_4x4(io, io);
  // io[0]: 00 10 20 30  01 11 21 31
  // io[1]: 02 12 22 32  03 13 23 33

  // stage 1
  {
    const __m128i cp_28_n4 = pair_set_epi16(cospi_28_64, -cospi_4_64);
    const __m128i cp_4_28 = pair_set_epi16(cospi_4_64, cospi_28_64);
    const __m128i cp_n20_12 = pair_set_epi16(-cospi_20_64, cospi_12_64);
    const __m128i cp_12_20 = pair_set_epi16(cospi_12_64, cospi_20_64);
    const __m128i lo_1 = _mm_unpackhi_epi16(io[0], zero);
    const __m128i lo_3 = _mm_unpackhi_epi16(io[1], zero);
    step1[4] = idct_calc_wraplow_sse2(cp_28_n4, cp_4_28, lo_1);    // step1 4&7
    step1[5] = idct_calc_wraplow_sse2(cp_n20_12, cp_12_20, lo_3);  // step1 5&6
  }

  // stage 2
  {
    const __m128i cp_24_n8 = pair_set_epi16(cospi_24_64, -cospi_8_64);
    const __m128i cp_8_24 = pair_set_epi16(cospi_8_64, cospi_24_64);
    const __m128i lo_0 = _mm_unpacklo_epi16(io[0], zero);
    const __m128i lo_2 = _mm_unpacklo_epi16(io[1], zero);
    const __m128i t = idct_madd_round_shift_sse2(cp_16_16, lo_0);
    step2[0] = _mm_packs_epi32(t, t);                            // step2 0&1
    step2[2] = idct_calc_wraplow_sse2(cp_8_24, cp_24_n8, lo_2);  // step2 3&2
    step2[4] = _mm_add_epi16(step1[4], step1[5]);                // step2 4&7
    step2[5] = _mm_sub_epi16(step1[4], step1[5]);                // step2 5&6
    step2[6] = _mm_unpackhi_epi64(step2[5], zero);               // step2 6
  }

  // stage 3
  {
    const __m128i lo_65 = _mm_unpacklo_epi16(step2[6], step2[5]);
    tmp[0] = _mm_add_epi16(step2[0], step2[2]);                     // step1 0&1
    tmp[1] = _mm_sub_epi16(step2[0], step2[2]);                     // step1 3&2
    step1[2] = _mm_unpackhi_epi64(tmp[1], tmp[0]);                  // step1 2&1
    step1[3] = _mm_unpacklo_epi64(tmp[1], tmp[0]);                  // step1 3&0
    step1[5] = idct_calc_wraplow_sse2(cp_16_n16, cp_16_16, lo_65);  // step1 5&6
  }

  // stage 4
  tmp[0] = _mm_add_epi16(step1[3], step2[4]);  // output 3&0
  tmp[1] = _mm_add_epi16(step1[2], step1[5]);  // output 2&1
  tmp[2] = _mm_sub_epi16(step1[3], step2[4]);  // output 4&7
  tmp[3] = _mm_sub_epi16(step1[2], step1[5]);  // output 5&6

  idct8x8_12_transpose_16bit_4x8(tmp, io);
  io[4] = io[5] = io[6] = io[7] = zero;

  idct8(io, io);
}

static INLINE void idct16_8col(__m128i *const io /*io[16]*/) {
  __m128i step1[16], step2[16];

  // stage 2
  butterfly(io[1], io[15], (int)cospi_30_64, (int)cospi_2_64, &step2[8],
            &step2[15]);
  butterfly(io[9], io[7], (int)cospi_14_64, (int)cospi_18_64, &step2[9],
            &step2[14]);
  butterfly(io[5], io[11], (int)cospi_22_64, (int)cospi_10_64, &step2[10],
            &step2[13]);
  butterfly(io[13], io[3], (int)cospi_6_64, (int)cospi_26_64, &step2[11],
            &step2[12]);

  // stage 3
  butterfly(io[2], io[14], (int)cospi_28_64, (int)cospi_4_64, &step1[4],
            &step1[7]);
  butterfly(io[10], io[6], (int)cospi_12_64, (int)cospi_20_64, &step1[5],
            &step1[6]);
  step1[8] = _mm_add_epi16(step2[8], step2[9]);
  step1[9] = _mm_sub_epi16(step2[8], step2[9]);
  step1[10] = _mm_sub_epi16(step2[11], step2[10]);
  step1[11] = _mm_add_epi16(step2[10], step2[11]);
  step1[12] = _mm_add_epi16(step2[12], step2[13]);
  step1[13] = _mm_sub_epi16(step2[12], step2[13]);
  step1[14] = _mm_sub_epi16(step2[15], step2[14]);
  step1[15] = _mm_add_epi16(step2[14], step2[15]);

  // stage 4
  butterfly(io[0], io[8], (int)cospi_16_64, (int)cospi_16_64, &step2[1],
            &step2[0]);
  butterfly(io[4], io[12], (int)cospi_24_64, (int)cospi_8_64, &step2[2],
            &step2[3]);
  butterfly(step1[14], step1[9], (int)cospi_24_64, (int)cospi_8_64, &step2[9],
            &step2[14]);
  butterfly(step1[10], step1[13], -(int)cospi_8_64, -(int)cospi_24_64,
            &step2[13], &step2[10]);
  step2[5] = _mm_sub_epi16(step1[4], step1[5]);
  step1[4] = _mm_add_epi16(step1[4], step1[5]);
  step2[6] = _mm_sub_epi16(step1[7], step1[6]);
  step1[7] = _mm_add_epi16(step1[6], step1[7]);
  step2[8] = step1[8];
  step2[11] = step1[11];
  step2[12] = step1[12];
  step2[15] = step1[15];

  // stage 5
  step1[0] = _mm_add_epi16(step2[0], step2[3]);
  step1[1] = _mm_add_epi16(step2[1], step2[2]);
  step1[2] = _mm_sub_epi16(step2[1], step2[2]);
  step1[3] = _mm_sub_epi16(step2[0], step2[3]);
  butterfly(step2[6], step2[5], (int)cospi_16_64, (int)cospi_16_64, &step1[5],
            &step1[6]);
  step1[8] = _mm_add_epi16(step2[8], step2[11]);
  step1[9] = _mm_add_epi16(step2[9], step2[10]);
  step1[10] = _mm_sub_epi16(step2[9], step2[10]);
  step1[11] = _mm_sub_epi16(step2[8], step2[11]);
  step1[12] = _mm_sub_epi16(step2[15], step2[12]);
  step1[13] = _mm_sub_epi16(step2[14], step2[13]);
  step1[14] = _mm_add_epi16(step2[14], step2[13]);
  step1[15] = _mm_add_epi16(step2[15], step2[12]);

  // stage 6
  step2[0] = _mm_add_epi16(step1[0], step1[7]);
  step2[1] = _mm_add_epi16(step1[1], step1[6]);
  step2[2] = _mm_add_epi16(step1[2], step1[5]);
  step2[3] = _mm_add_epi16(step1[3], step1[4]);
  step2[4] = _mm_sub_epi16(step1[3], step1[4]);
  step2[5] = _mm_sub_epi16(step1[2], step1[5]);
  step2[6] = _mm_sub_epi16(step1[1], step1[6]);
  step2[7] = _mm_sub_epi16(step1[0], step1[7]);
  butterfly(step1[13], step1[10], (int)cospi_16_64, (int)cospi_16_64,
            &step2[10], &step2[13]);
  butterfly(step1[12], step1[11], (int)cospi_16_64, (int)cospi_16_64,
            &step2[11], &step2[12]);

  // stage 7
  io[0] = _mm_add_epi16(step2[0], step1[15]);
  io[1] = _mm_add_epi16(step2[1], step1[14]);
  io[2] = _mm_add_epi16(step2[2], step2[13]);
  io[3] = _mm_add_epi16(step2[3], step2[12]);
  io[4] = _mm_add_epi16(step2[4], step2[11]);
  io[5] = _mm_add_epi16(step2[5], step2[10]);
  io[6] = _mm_add_epi16(step2[6], step1[9]);
  io[7] = _mm_add_epi16(step2[7], step1[8]);
  io[8] = _mm_sub_epi16(step2[7], step1[8]);
  io[9] = _mm_sub_epi16(step2[6], step1[9]);
  io[10] = _mm_sub_epi16(step2[5], step2[10]);
  io[11] = _mm_sub_epi16(step2[4], step2[11]);
  io[12] = _mm_sub_epi16(step2[3], step2[12]);
  io[13] = _mm_sub_epi16(step2[2], step2[13]);
  io[14] = _mm_sub_epi16(step2[1], step1[14]);
  io[15] = _mm_sub_epi16(step2[0], step1[15]);
}

static INLINE void idct16x16_10_pass1(const __m128i *const input /*input[4]*/,
                                      __m128i *const output /*output[16]*/) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i k__cospi_p16_p16 = pair_set_epi16(cospi_16_64, cospi_16_64);
  const __m128i k__cospi_m16_p16 = pair_set_epi16(-cospi_16_64, cospi_16_64);
  __m128i step1[16], step2[16];

  transpose_16bit_4x4(input, output);

  // stage 2
  {
    const __m128i k__cospi_p30_m02 = pair_set_epi16(cospi_30_64, -cospi_2_64);
    const __m128i k__cospi_p02_p30 = pair_set_epi16(cospi_2_64, cospi_30_64);
    const __m128i k__cospi_p06_m26 = pair_set_epi16(cospi_6_64, -cospi_26_64);
    const __m128i k__cospi_p26_p06 = pair_set_epi16(cospi_26_64, cospi_6_64);
    const __m128i lo_1_15 = _mm_unpackhi_epi16(output[0], zero);
    const __m128i lo_13_3 = _mm_unpackhi_epi16(zero, output[1]);
    step2[8] = idct_calc_wraplow_sse2(k__cospi_p30_m02, k__cospi_p02_p30,
                                      lo_1_15);  // step2 8&15
    step2[11] = idct_calc_wraplow_sse2(k__cospi_p06_m26, k__cospi_p26_p06,
                                       lo_13_3);  // step2 11&12
  }

  // stage 3
  {
    const __m128i k__cospi_p28_m04 = pair_set_epi16(cospi_28_64, -cospi_4_64);
    const __m128i k__cospi_p04_p28 = pair_set_epi16(cospi_4_64, cospi_28_64);
    const __m128i lo_2_14 = _mm_unpacklo_epi16(output[1], zero);
    step1[4] = idct_calc_wraplow_sse2(k__cospi_p28_m04, k__cospi_p04_p28,
                                      lo_2_14);  // step1 4&7
    step1[13] = _mm_unpackhi_epi64(step2[11], zero);
    step1[14] = _mm_unpackhi_epi64(step2[8], zero);
  }

  // stage 4
  {
    const __m128i k__cospi_m08_p24 = pair_set_epi16(-cospi_8_64, cospi_24_64);
    const __m128i k__cospi_p24_p08 = pair_set_epi16(cospi_24_64, cospi_8_64);
    const __m128i k__cospi_m24_m08 = pair_set_epi16(-cospi_24_64, -cospi_8_64);
    const __m128i lo_0_8 = _mm_unpacklo_epi16(output[0], zero);
    const __m128i lo_9_14 = _mm_unpacklo_epi16(step2[8], step1[14]);
    const __m128i lo_10_13 = _mm_unpacklo_epi16(step2[11], step1[13]);
    const __m128i t = idct_madd_round_shift_sse2(lo_0_8, k__cospi_p16_p16);
    step1[0] = _mm_packs_epi32(t, t);  // step2 0&1
    step2[9] = idct_calc_wraplow_sse2(k__cospi_m08_p24, k__cospi_p24_p08,
                                      lo_9_14);  // step2 9&14
    step2[10] = idct_calc_wraplow_sse2(k__cospi_m24_m08, k__cospi_m08_p24,
                                       lo_10_13);  // step2 10&13
    step2[6] = _mm_unpackhi_epi64(step1[4], zero);
  }

  // stage 5
  {
    const __m128i lo_5_6 = _mm_unpacklo_epi16(step1[4], step2[6]);
    step1[6] = idct_calc_wraplow_sse2(k__cospi_p16_p16, k__cospi_m16_p16,
                                      lo_5_6);  // step1 6&5
    step1[8] = _mm_add_epi16(step2[8], step2[11]);
    step1[9] = _mm_add_epi16(step2[9], step2[10]);
    step1[10] = _mm_sub_epi16(step2[9], step2[10]);
    step1[11] = _mm_sub_epi16(step2[8], step2[11]);
    step1[12] = _mm_unpackhi_epi64(step1[11], zero);
    step1[13] = _mm_unpackhi_epi64(step1[10], zero);
    step1[14] = _mm_unpackhi_epi64(step1[9], zero);
    step1[15] = _mm_unpackhi_epi64(step1[8], zero);
  }

  // stage 6
  {
    const __m128i lo_10_13 = _mm_unpacklo_epi16(step1[10], step1[13]);
    const __m128i lo_11_12 = _mm_unpacklo_epi16(step1[11], step1[12]);
    step2[10] = idct_calc_wraplow_sse2(k__cospi_m16_p16, k__cospi_p16_p16,
                                       lo_10_13);  // step2 10&13
    step2[11] = idct_calc_wraplow_sse2(k__cospi_m16_p16, k__cospi_p16_p16,
                                       lo_11_12);  // step2 11&12
    step2[13] = _mm_unpackhi_epi64(step2[10], zero);
    step2[12] = _mm_unpackhi_epi64(step2[11], zero);
    step2[3] = _mm_add_epi16(step1[0], step1[4]);
    step2[1] = _mm_add_epi16(step1[0], step1[6]);
    step2[6] = _mm_sub_epi16(step1[0], step1[6]);
    step2[4] = _mm_sub_epi16(step1[0], step1[4]);
    step2[0] = _mm_unpackhi_epi64(step2[3], zero);
    step2[2] = _mm_unpackhi_epi64(step2[1], zero);
    step2[5] = _mm_unpackhi_epi64(step2[6], zero);
    step2[7] = _mm_unpackhi_epi64(step2[4], zero);
  }

  // stage 7. Left 8x16 only.
  output[0] = _mm_add_epi16(step2[0], step1[15]);
  output[1] = _mm_add_epi16(step2[1], step1[14]);
  output[2] = _mm_add_epi16(step2[2], step2[13]);
  output[3] = _mm_add_epi16(step2[3], step2[12]);
  output[4] = _mm_add_epi16(step2[4], step2[11]);
  output[5] = _mm_add_epi16(step2[5], step2[10]);
  output[6] = _mm_add_epi16(step2[6], step1[9]);
  output[7] = _mm_add_epi16(step2[7], step1[8]);
  output[8] = _mm_sub_epi16(step2[7], step1[8]);
  output[9] = _mm_sub_epi16(step2[6], step1[9]);
  output[10] = _mm_sub_epi16(step2[5], step2[10]);
  output[11] = _mm_sub_epi16(step2[4], step2[11]);
  output[12] = _mm_sub_epi16(step2[3], step2[12]);
  output[13] = _mm_sub_epi16(step2[2], step2[13]);
  output[14] = _mm_sub_epi16(step2[1], step1[14]);
  output[15] = _mm_sub_epi16(step2[0], step1[15]);
}

static INLINE void idct16x16_10_pass2(__m128i *const l /*l[8]*/,
                                      __m128i *const io /*io[16]*/) {
  const __m128i zero = _mm_setzero_si128();
  __m128i step1[16], step2[16];

  transpose_16bit_4x8(l, io);

  // stage 2
  butterfly(io[1], zero, (int)cospi_30_64, (int)cospi_2_64, &step2[8],
            &step2[15]);
  butterfly(zero, io[3], (int)cospi_6_64, (int)cospi_26_64, &step2[11],
            &step2[12]);

  // stage 3
  butterfly(io[2], zero, (int)cospi_28_64, (int)cospi_4_64, &step1[4],
            &step1[7]);

  // stage 4
  butterfly(io[0], zero, (int)cospi_16_64, (int)cospi_16_64, &step1[1],
            &step1[0]);
  butterfly(step2[15], step2[8], (int)cospi_24_64, (int)cospi_8_64, &step2[9],
            &step2[14]);
  butterfly(step2[11], step2[12], -(int)cospi_8_64, -(int)cospi_24_64,
            &step2[13], &step2[10]);

  // stage 5
  butterfly(step1[7], step1[4], (int)cospi_16_64, (int)cospi_16_64, &step1[5],
            &step1[6]);
  step1[8] = _mm_add_epi16(step2[8], step2[11]);
  step1[9] = _mm_add_epi16(step2[9], step2[10]);
  step1[10] = _mm_sub_epi16(step2[9], step2[10]);
  step1[11] = _mm_sub_epi16(step2[8], step2[11]);
  step1[12] = _mm_sub_epi16(step2[15], step2[12]);
  step1[13] = _mm_sub_epi16(step2[14], step2[13]);
  step1[14] = _mm_add_epi16(step2[14], step2[13]);
  step1[15] = _mm_add_epi16(step2[15], step2[12]);

  // stage 6
  step2[0] = _mm_add_epi16(step1[0], step1[7]);
  step2[1] = _mm_add_epi16(step1[1], step1[6]);
  step2[2] = _mm_add_epi16(step1[1], step1[5]);
  step2[3] = _mm_add_epi16(step1[0], step1[4]);
  step2[4] = _mm_sub_epi16(step1[0], step1[4]);
  step2[5] = _mm_sub_epi16(step1[1], step1[5]);
  step2[6] = _mm_sub_epi16(step1[1], step1[6]);
  step2[7] = _mm_sub_epi16(step1[0], step1[7]);
  butterfly(step1[13], step1[10], (int)cospi_16_64, (int)cospi_16_64,
            &step2[10], &step2[13]);
  butterfly(step1[12], step1[11], (int)cospi_16_64, (int)cospi_16_64,
            &step2[11], &step2[12]);

  // stage 7
  io[0] = _mm_add_epi16(step2[0], step1[15]);
  io[1] = _mm_add_epi16(step2[1], step1[14]);
  io[2] = _mm_add_epi16(step2[2], step2[13]);
  io[3] = _mm_add_epi16(step2[3], step2[12]);
  io[4] = _mm_add_epi16(step2[4], step2[11]);
  io[5] = _mm_add_epi16(step2[5], step2[10]);
  io[6] = _mm_add_epi16(step2[6], step1[9]);
  io[7] = _mm_add_epi16(step2[7], step1[8]);
  io[8] = _mm_sub_epi16(step2[7], step1[8]);
  io[9] = _mm_sub_epi16(step2[6], step1[9]);
  io[10] = _mm_sub_epi16(step2[5], step2[10]);
  io[11] = _mm_sub_epi16(step2[4], step2[11]);
  io[12] = _mm_sub_epi16(step2[3], step2[12]);
  io[13] = _mm_sub_epi16(step2[2], step2[13]);
  io[14] = _mm_sub_epi16(step2[1], step1[14]);
  io[15] = _mm_sub_epi16(step2[0], step1[15]);
}

void idct4_sse2(__m128i *in);
void idct8_sse2(__m128i *in);
void idct16_sse2(__m128i *in0, __m128i *in1);
void iadst4_sse2(__m128i *in);
void iadst8_sse2(__m128i *in);
void iadst16_sse2(__m128i *in0, __m128i *in1);

#endif  // VPX_DSP_X86_INV_TXFM_SSE2_H_
