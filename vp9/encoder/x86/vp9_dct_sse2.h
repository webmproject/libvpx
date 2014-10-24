/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_X86_VP9_DCT_SSE2_H_
#define VP9_ENCODER_X86_VP9_DCT_SSE2_H_

#ifdef __cplusplus
extern "C" {
#endif

#define pair_set_epi32(a, b) \
  _mm_set_epi32(b, a, b, a)

void vp9_fdct4x4_sse2(const int16_t *input, tran_low_t *output, int stride);
void vp9_fdct8x8_sse2(const int16_t *input, tran_low_t *output, int stride);
void vp9_fdct16x16_sse2(const int16_t *input, tran_low_t *output, int stride);
void vp9_highbd_fdct4x4_sse2(const int16_t *input, tran_low_t *output,
                            int stride);
void vp9_highbd_fdct8x8_sse2(const int16_t *input, tran_low_t *output,
                            int stride);
void vp9_highbd_fdct16x16_sse2(const int16_t *input, tran_low_t *output,
                            int stride);

static INLINE __m128i k_madd_epi32(__m128i a, __m128i b) {
  __m128i buf0, buf1;
  buf0 = _mm_mul_epu32(a, b);
  a = _mm_srli_epi64(a, 32);
  b = _mm_srli_epi64(b, 32);
  buf1 = _mm_mul_epu32(a, b);
  return _mm_add_epi64(buf0, buf1);
}

static INLINE __m128i k_packs_epi64(__m128i a, __m128i b) {
  __m128i buf0 = _mm_shuffle_epi32(a, _MM_SHUFFLE(0, 0, 2, 0));
  __m128i buf1 = _mm_shuffle_epi32(b, _MM_SHUFFLE(0, 0, 2, 0));
  return _mm_unpacklo_epi64(buf0, buf1);
}

static INLINE int check_epi16_overflow_x2(__m128i reg0, __m128i reg1) {
  const __m128i max_overflow = _mm_set1_epi16(0x7fff);
  const __m128i min_overflow = _mm_set1_epi16(0x8000);
  __m128i cmp0 = _mm_or_si128(_mm_cmpeq_epi16(reg0, max_overflow),
                              _mm_cmpeq_epi16(reg0, min_overflow));
  __m128i cmp1 = _mm_or_si128(_mm_cmpeq_epi16(reg1, max_overflow),
                              _mm_cmpeq_epi16(reg1, min_overflow));
  cmp0 = _mm_or_si128(cmp0, cmp1);
  return _mm_movemask_epi8(cmp0);
}

static INLINE int check_epi16_overflow_x4(__m128i reg0, __m128i reg1,
                                   __m128i reg2, __m128i reg3) {
  const __m128i max_overflow = _mm_set1_epi16(0x7fff);
  const __m128i min_overflow = _mm_set1_epi16(0x8000);
  __m128i cmp0 = _mm_or_si128(_mm_cmpeq_epi16(reg0, max_overflow),
                              _mm_cmpeq_epi16(reg0, min_overflow));
  __m128i cmp1 = _mm_or_si128(_mm_cmpeq_epi16(reg1, max_overflow),
                              _mm_cmpeq_epi16(reg1, min_overflow));
  __m128i cmp2 = _mm_or_si128(_mm_cmpeq_epi16(reg2, max_overflow),
                              _mm_cmpeq_epi16(reg2, min_overflow));
  __m128i cmp3 = _mm_or_si128(_mm_cmpeq_epi16(reg3, max_overflow),
                              _mm_cmpeq_epi16(reg3, min_overflow));
  cmp0 = _mm_or_si128(_mm_or_si128(cmp0, cmp1), _mm_or_si128(cmp2, cmp3));
  return _mm_movemask_epi8(cmp0);
}

static INLINE int check_epi16_overflow_x8(__m128i reg0, __m128i reg1,
                                   __m128i reg2, __m128i reg3, __m128i reg4,
                                   __m128i reg5, __m128i reg6, __m128i reg7) {
  int res0, res1;
  res0 = check_epi16_overflow_x4(reg0, reg1, reg2, reg3);
  res1 = check_epi16_overflow_x4(reg4, reg5, reg6, reg7);
  return res0 + res1;
}

static INLINE int check_epi16_overflow_x12(__m128i reg0, __m128i reg1,
                                   __m128i reg2, __m128i reg3, __m128i reg4,
                                   __m128i reg5, __m128i reg6, __m128i reg7,
                                   __m128i reg8, __m128i reg9, __m128i reg10,
                                   __m128i reg11) {
  int res0, res1;
  res0 = check_epi16_overflow_x4(reg0, reg1, reg2, reg3);
  res1 = check_epi16_overflow_x4(reg4, reg5, reg6, reg7);
  if (!res0)
    res0 = check_epi16_overflow_x4(reg8, reg9, reg10, reg11);
  return res0 + res1;
}

static INLINE int check_epi16_overflow_x16(__m128i reg0,  __m128i reg1,
                                    __m128i reg2, __m128i reg3,  __m128i reg4,
                                    __m128i reg5, __m128i reg6,  __m128i reg7,
                                    __m128i reg8, __m128i reg9,  __m128i reg10,
                                    __m128i reg11, __m128i reg12, __m128i reg13,
                                    __m128i reg14, __m128i reg15) {
  int res0, res1;
  res0 = check_epi16_overflow_x4(reg0, reg1, reg2, reg3);
  res1 = check_epi16_overflow_x4(reg4, reg5, reg6, reg7);
  if (!res0) {
    res0 = check_epi16_overflow_x4(reg8, reg9, reg10, reg11);
    if (!res1)
      res1 = check_epi16_overflow_x4(reg12, reg13, reg14, reg15);
  }
  return res0 + res1;
}

static INLINE int check_epi16_overflow_x32(__m128i reg0,  __m128i reg1,
                                __m128i reg2, __m128i reg3,  __m128i reg4,
                                __m128i reg5, __m128i reg6,  __m128i reg7,
                                __m128i reg8, __m128i reg9,  __m128i reg10,
                                __m128i reg11, __m128i reg12, __m128i reg13,
                                __m128i reg14, __m128i reg15, __m128i reg16,
                                __m128i reg17, __m128i reg18, __m128i reg19,
                                __m128i reg20, __m128i reg21, __m128i reg22,
                                __m128i reg23, __m128i reg24, __m128i reg25,
                                __m128i reg26, __m128i reg27, __m128i reg28,
                                __m128i reg29, __m128i reg30, __m128i reg31) {
  int res0, res1;
  res0 = check_epi16_overflow_x4(reg0, reg1, reg2, reg3);
  res1 = check_epi16_overflow_x4(reg4, reg5, reg6, reg7);
  if (!res0) {
    res0 = check_epi16_overflow_x4(reg8, reg9, reg10, reg11);
    if (!res1) {
      res1 = check_epi16_overflow_x4(reg12, reg13, reg14, reg15);
      if (!res0) {
        res0 = check_epi16_overflow_x4(reg16, reg17, reg18, reg19);
        if (!res1) {
          res1 = check_epi16_overflow_x4(reg20, reg21, reg22, reg23);
          if (!res0) {
            res0 = check_epi16_overflow_x4(reg24, reg25, reg26, reg27);
            if (!res1)
              res1 = check_epi16_overflow_x4(reg28, reg29, reg30, reg31);
          }
        }
      }
    }
  }
  return res0 + res1;
}

static INLINE int k_check_epi32_overflow_4(__m128i reg0, __m128i reg1,
                 __m128i reg2, __m128i reg3, const __m128i* zero) {
  __m128i minus_one = _mm_set1_epi32(-1);
  // Check for overflows
  __m128i reg0_shifted = _mm_slli_epi64(reg0, 1);
  __m128i reg1_shifted = _mm_slli_epi64(reg1, 1);
  __m128i reg2_shifted = _mm_slli_epi64(reg2, 1);
  __m128i reg3_shifted = _mm_slli_epi64(reg3, 1);
  __m128i reg0_top_dwords = _mm_shuffle_epi32(
                                reg0_shifted, _MM_SHUFFLE(0, 0, 3, 1));
  __m128i reg1_top_dwords = _mm_shuffle_epi32(
                                reg1_shifted, _MM_SHUFFLE(0, 0, 3, 1));
  __m128i reg2_top_dwords = _mm_shuffle_epi32(
                                reg2_shifted, _MM_SHUFFLE(0, 0, 3, 1));
  __m128i reg3_top_dwords = _mm_shuffle_epi32(
                                reg3_shifted, _MM_SHUFFLE(0, 0, 3, 1));
  __m128i top_dwords_01 = _mm_unpacklo_epi64(reg0_top_dwords, reg1_top_dwords);
  __m128i top_dwords_23 = _mm_unpacklo_epi64(reg2_top_dwords, reg3_top_dwords);
  __m128i valid_positve_01 = _mm_cmpeq_epi32(top_dwords_01, *zero);
  __m128i valid_positve_23 = _mm_cmpeq_epi32(top_dwords_23, *zero);
  __m128i valid_negative_01 = _mm_cmpeq_epi32(top_dwords_01, minus_one);
  __m128i valid_negative_23 = _mm_cmpeq_epi32(top_dwords_23, minus_one);
  int overflow_01 = _mm_movemask_epi8(
            _mm_cmpeq_epi32(valid_positve_01, valid_negative_01));
  int overflow_23 = _mm_movemask_epi8(
            _mm_cmpeq_epi32(valid_positve_23, valid_negative_23));
  return (overflow_01 + overflow_23);
}

static INLINE int k_check_epi32_overflow_8(__m128i reg0, __m128i reg1,
                 __m128i reg2, __m128i reg3, __m128i reg4, __m128i reg5,
                 __m128i reg6, __m128i reg7, const __m128i* zero) {
  int overflow = k_check_epi32_overflow_4(reg0, reg1, reg2, reg3, zero);
  if (!overflow) {
    overflow = k_check_epi32_overflow_4(reg4, reg5, reg6, reg7, zero);
  }
  return overflow;
}

static INLINE int k_check_epi32_overflow_16(__m128i reg0, __m128i reg1,
                 __m128i reg2, __m128i reg3, __m128i reg4, __m128i reg5,
                 __m128i reg6, __m128i reg7, __m128i reg8, __m128i reg9,
                 __m128i reg10, __m128i reg11, __m128i reg12, __m128i reg13,
                 __m128i reg14, __m128i reg15, const __m128i* zero) {
  int overflow = k_check_epi32_overflow_4(reg0, reg1, reg2, reg3, zero);
  if (!overflow) {
    overflow = k_check_epi32_overflow_4(reg4, reg5, reg6, reg7, zero);
    if (!overflow) {
      overflow = k_check_epi32_overflow_4(reg8, reg9, reg10, reg11, zero);
      if (!overflow) {
        overflow = k_check_epi32_overflow_4(reg12, reg13, reg14, reg15, zero);
      }
    }
  }
  return overflow;
}

static INLINE int k_check_epi32_overflow_32(__m128i reg0, __m128i reg1,
                 __m128i reg2, __m128i reg3, __m128i reg4, __m128i reg5,
                 __m128i reg6, __m128i reg7, __m128i reg8, __m128i reg9,
                 __m128i reg10, __m128i reg11, __m128i reg12, __m128i reg13,
                 __m128i reg14, __m128i reg15, __m128i reg16, __m128i reg17,
                 __m128i reg18, __m128i reg19, __m128i reg20, __m128i reg21,
                 __m128i reg22, __m128i reg23, __m128i reg24, __m128i reg25,
                 __m128i reg26, __m128i reg27, __m128i reg28, __m128i reg29,
                 __m128i reg30, __m128i reg31, const __m128i* zero) {
  int overflow = k_check_epi32_overflow_4(reg0, reg1, reg2, reg3, zero);
  if (!overflow) {
    overflow = k_check_epi32_overflow_4(reg4, reg5, reg6, reg7, zero);
    if (!overflow) {
      overflow = k_check_epi32_overflow_4(reg8, reg9, reg10, reg11, zero);
      if (!overflow) {
        overflow = k_check_epi32_overflow_4(reg12, reg13, reg14, reg15, zero);
        if (!overflow) {
          overflow = k_check_epi32_overflow_4(reg16, reg17, reg18, reg19, zero);
          if (!overflow) {
            overflow = k_check_epi32_overflow_4(reg20, reg21,
                                                reg22, reg23, zero);
            if (!overflow) {
              overflow = k_check_epi32_overflow_4(reg24, reg25,
                                                  reg26, reg27, zero);
              if (!overflow) {
                overflow = k_check_epi32_overflow_4(reg28, reg29,
                                                    reg30, reg31, zero);
              }
            }
          }
        }
      }
    }
  }
  return overflow;
}

static INLINE void store_output(const __m128i output, tran_low_t* dst_ptr) {
#if CONFIG_VP9_HIGHBITDEPTH
  const __m128i zero = _mm_setzero_si128();
  const __m128i sign_bits = _mm_cmplt_epi16(output, zero);
  __m128i out0 = _mm_unpacklo_epi16(output, sign_bits);
  __m128i out1 = _mm_unpackhi_epi16(output, sign_bits);
  _mm_store_si128((__m128i *)(dst_ptr), out0);
  _mm_store_si128((__m128i *)(dst_ptr + 4), out1);
#else
  _mm_store_si128((__m128i *)(dst_ptr), output);
#endif
}

static INLINE void storeu_output(const __m128i output, tran_low_t* dst_ptr) {
#if CONFIG_VP9_HIGHBITDEPTH
  const __m128i zero = _mm_setzero_si128();
  const __m128i sign_bits = _mm_cmplt_epi16(output, zero);
  __m128i out0 = _mm_unpacklo_epi16(output, sign_bits);
  __m128i out1 = _mm_unpackhi_epi16(output, sign_bits);
  _mm_storeu_si128((__m128i *)(dst_ptr), out0);
  _mm_storeu_si128((__m128i *)(dst_ptr + 4), out1);
#else
  _mm_storeu_si128((__m128i *)(dst_ptr), output);
#endif
}


static INLINE __m128i mult_round_shift(const __m128i in0, const __m128i in1,
                                       const __m128i multiplier,
                                       const __m128i rounding,
                                       const int shift) {
  const __m128i u0 = _mm_madd_epi16(in0, multiplier);
  const __m128i u1 = _mm_madd_epi16(in1, multiplier);
  const __m128i v0 = _mm_add_epi32(u0, rounding);
  const __m128i v1 = _mm_add_epi32(u1, rounding);
  const __m128i w0 = _mm_srai_epi32(v0, shift);
  const __m128i w1 = _mm_srai_epi32(v1, shift);
  return _mm_packs_epi32(w0, w1);
}

static INLINE void transpose_and_output8x8(
                                const __m128i in00, const __m128i in01,
                                const __m128i in02, const __m128i in03,
                                const __m128i in04, const __m128i in05,
                                const __m128i in06, const __m128i in07,
                                const int pass, int16_t* out0_ptr,
                                tran_low_t* out1_ptr) {
  // 00 01 02 03 04 05 06 07
  // 10 11 12 13 14 15 16 17
  // 20 21 22 23 24 25 26 27
  // 30 31 32 33 34 35 36 37
  // 40 41 42 43 44 45 46 47
  // 50 51 52 53 54 55 56 57
  // 60 61 62 63 64 65 66 67
  // 70 71 72 73 74 75 76 77
  const __m128i tr0_0 = _mm_unpacklo_epi16(in00, in01);
  const __m128i tr0_1 = _mm_unpacklo_epi16(in02, in03);
  const __m128i tr0_2 = _mm_unpackhi_epi16(in00, in01);
  const __m128i tr0_3 = _mm_unpackhi_epi16(in02, in03);
  const __m128i tr0_4 = _mm_unpacklo_epi16(in04, in05);
  const __m128i tr0_5 = _mm_unpacklo_epi16(in06, in07);
  const __m128i tr0_6 = _mm_unpackhi_epi16(in04, in05);
  const __m128i tr0_7 = _mm_unpackhi_epi16(in06, in07);
  // 00 10 01 11 02 12 03 13
  // 20 30 21 31 22 32 23 33
  // 04 14 05 15 06 16 07 17
  // 24 34 25 35 26 36 27 37
  // 40 50 41 51 42 52 43 53
  // 60 70 61 71 62 72 63 73
  // 54 54 55 55 56 56 57 57
  // 64 74 65 75 66 76 67 77
  const __m128i tr1_0 = _mm_unpacklo_epi32(tr0_0, tr0_1);
  const __m128i tr1_1 = _mm_unpacklo_epi32(tr0_2, tr0_3);
  const __m128i tr1_2 = _mm_unpackhi_epi32(tr0_0, tr0_1);
  const __m128i tr1_3 = _mm_unpackhi_epi32(tr0_2, tr0_3);
  const __m128i tr1_4 = _mm_unpacklo_epi32(tr0_4, tr0_5);
  const __m128i tr1_5 = _mm_unpacklo_epi32(tr0_6, tr0_7);
  const __m128i tr1_6 = _mm_unpackhi_epi32(tr0_4, tr0_5);
  const __m128i tr1_7 = _mm_unpackhi_epi32(tr0_6, tr0_7);
  // 00 10 20 30 01 11 21 31
  // 40 50 60 70 41 51 61 71
  // 02 12 22 32 03 13 23 33
  // 42 52 62 72 43 53 63 73
  // 04 14 24 34 05 15 21 36
  // 44 54 64 74 45 55 61 76
  // 06 16 26 36 07 17 27 37
  // 46 56 66 76 47 57 67 77
  const __m128i tr2_0 = _mm_unpacklo_epi64(tr1_0, tr1_4);
  const __m128i tr2_1 = _mm_unpackhi_epi64(tr1_0, tr1_4);
  const __m128i tr2_2 = _mm_unpacklo_epi64(tr1_2, tr1_6);
  const __m128i tr2_3 = _mm_unpackhi_epi64(tr1_2, tr1_6);
  const __m128i tr2_4 = _mm_unpacklo_epi64(tr1_1, tr1_5);
  const __m128i tr2_5 = _mm_unpackhi_epi64(tr1_1, tr1_5);
  const __m128i tr2_6 = _mm_unpacklo_epi64(tr1_3, tr1_7);
  const __m128i tr2_7 = _mm_unpackhi_epi64(tr1_3, tr1_7);
  // 00 10 20 30 40 50 60 70
  // 01 11 21 31 41 51 61 71
  // 02 12 22 32 42 52 62 72
  // 03 13 23 33 43 53 63 73
  // 04 14 24 34 44 54 64 74
  // 05 15 25 35 45 55 65 75
  // 06 16 26 36 46 56 66 76
  // 07 17 27 37 47 57 67 77
  if (pass == 0) {
    _mm_storeu_si128((__m128i*)(out0_ptr + 0 * 16), tr2_0);
    _mm_storeu_si128((__m128i*)(out0_ptr + 1 * 16), tr2_1);
    _mm_storeu_si128((__m128i*)(out0_ptr + 2 * 16), tr2_2);
    _mm_storeu_si128((__m128i*)(out0_ptr + 3 * 16), tr2_3);
    _mm_storeu_si128((__m128i*)(out0_ptr + 4 * 16), tr2_4);
    _mm_storeu_si128((__m128i*)(out0_ptr + 5 * 16), tr2_5);
    _mm_storeu_si128((__m128i*)(out0_ptr + 6 * 16), tr2_6);
    _mm_storeu_si128((__m128i*)(out0_ptr + 7 * 16), tr2_7);
  } else {
    storeu_output(tr2_0, (out1_ptr + 0 * 16));
    storeu_output(tr2_1, (out1_ptr + 1 * 16));
    storeu_output(tr2_2, (out1_ptr + 2 * 16));
    storeu_output(tr2_3, (out1_ptr + 3 * 16));
    storeu_output(tr2_4, (out1_ptr + 4 * 16));
    storeu_output(tr2_5, (out1_ptr + 5 * 16));
    storeu_output(tr2_6, (out1_ptr + 6 * 16));
    storeu_output(tr2_7, (out1_ptr + 7 * 16));
  }
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_X86_VP9_DCT_SSE2_H_
