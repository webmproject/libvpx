#ifndef VP10_TXMF1D_SSE2_H_
#define VP10_TXMF1D_SSE2_H_

#include <emmintrin.h>
#include "vp10/common/vp10_txfm.h"

#ifdef __cplusplus
extern "C" {
#endif

void vp10_fdct4_new_sse2(const __m128i* input, __m128i* output,
                         const int8_t* cos_bit, const int8_t* stage_range);
void vp10_fdct8_new_sse2(const __m128i* input, __m128i* output,
                         const int8_t* cos_bit, const int8_t* stage_range);
void vp10_fdct16_new_sse2(const __m128i* input, __m128i* output,
                          const int8_t* cos_bit, const int8_t* stage_range);
void vp10_fdct32_new_sse2(const __m128i* input, __m128i* output,
                          const int8_t* cos_bit, const int8_t* stage_range);
void vp10_fdct64_new_sse2(const __m128i* input, __m128i* output,
                          const int8_t* cos_bit, const int8_t* stage_range);

void vp10_fadst4_new_sse2(const __m128i* input, __m128i* output,
                          const int8_t* cos_bit, const int8_t* stage_range);
void vp10_fadst8_new_sse2(const __m128i* input, __m128i* output,
                          const int8_t* cos_bit, const int8_t* stage_range);
void vp10_fadst16_new_sse2(const __m128i* input, __m128i* output,
                           const int8_t* cos_bit, const int8_t* stage_range);
void vp10_fadst32_new_sse2(const __m128i* input, __m128i* output,
                           const int8_t* cos_bit, const int8_t* stage_range);

void vp10_idct4_new_sse2(const __m128i* input, __m128i* output,
                         const int8_t* cos_bit, const int8_t* stage_range);
void vp10_idct8_new_sse2(const __m128i* input, __m128i* output,
                         const int8_t* cos_bit, const int8_t* stage_range);
void vp10_idct16_new_sse2(const __m128i* input, __m128i* output,
                          const int8_t* cos_bit, const int8_t* stage_range);
void vp10_idct32_new_sse2(const __m128i* input, __m128i* output,
                          const int8_t* cos_bit, const int8_t* stage_range);
void vp10_idct64_new_sse2(const __m128i* input, __m128i* output,
                          const int8_t* cos_bit, const int8_t* stage_range);

void vp10_iadst4_new_sse2(const __m128i* input, __m128i* output,
                          const int8_t* cos_bit, const int8_t* stage_range);
void vp10_iadst8_new_sse2(const __m128i* input, __m128i* output,
                          const int8_t* cos_bit, const int8_t* stage_range);
void vp10_iadst16_new_sse2(const __m128i* input, __m128i* output,
                           const int8_t* cos_bit, const int8_t* stage_range);
void vp10_iadst32_new_sse2(const __m128i* input, __m128i* output,
                           const int8_t* cos_bit, const int8_t* stage_range);

static INLINE void transpose_32_4x4(int stride, const __m128i* input,
                                    __m128i* output) {
  __m128i temp0 = _mm_unpacklo_epi32(input[0 * stride], input[2 * stride]);
  __m128i temp1 = _mm_unpackhi_epi32(input[0 * stride], input[2 * stride]);
  __m128i temp2 = _mm_unpacklo_epi32(input[1 * stride], input[3 * stride]);
  __m128i temp3 = _mm_unpackhi_epi32(input[1 * stride], input[3 * stride]);

  output[0 * stride] = _mm_unpacklo_epi32(temp0, temp2);
  output[1 * stride] = _mm_unpackhi_epi32(temp0, temp2);
  output[2 * stride] = _mm_unpacklo_epi32(temp1, temp3);
  output[3 * stride] = _mm_unpackhi_epi32(temp1, temp3);
}

// the entire input block can be represent by a grid of 4x4 blocks
// each 4x4 blocks can be represent by 4 vertical __m128i
// we first transpose each 4x4 block internally
// than transpose the grid
static INLINE void transpose_32(int txfm_size, const __m128i* input,
                                __m128i* output) {
  const int num_per_128 = 4;
  const int row_size = txfm_size;
  const int col_size = txfm_size / num_per_128;
  int r, c;

  // transpose each 4x4 block internally
  for (r = 0; r < row_size; r += 4) {
    for (c = 0; c < col_size; c++) {
      transpose_32_4x4(col_size, &input[r * col_size + c],
                       &output[c * 4 * col_size + r / 4]);
    }
  }
}

#define mullo_epi32(a, b)                                                   \
  ({                                                                        \
    __m128i tmp1 = _mm_mul_epu32(a, b);                                       \
    __m128i tmp2 = _mm_mul_epu32(_mm_srli_si128(a, 4), _mm_srli_si128(b, 4)); \
    _mm_unpacklo_epi32(_mm_shuffle_epi32(tmp1, _MM_SHUFFLE(0, 0, 2, 0)),    \
                       _mm_shuffle_epi32(tmp2, _MM_SHUFFLE(0, 0, 2, 0)));   \
  })

#define round_shift_32_simple_sse2(input, bit)          \
  ({                                                    \
    __m128i round = _mm_set1_epi32((1 << (bit - 1)) - 1); \
    __m128i tmp1 = _mm_add_epi32(input, round);           \
    _mm_srai_epi32(tmp1, bit);                          \
  })

#define round_shift_32_sse2(vec, bit)             \
  ({                                              \
    __m128i sign, tmp, round;                       \
    sign = _mm_srai_epi32(vec, 31);               \
    tmp = _mm_add_epi32(vec, sign);               \
    tmp = _mm_xor_si128(tmp, sign);               \
    round = _mm_set1_epi32((1 << (bit - 1)) - 1); \
    tmp = _mm_add_epi32(tmp, round);              \
    tmp = _mm_srli_epi32(tmp, bit);               \
    tmp = _mm_xor_si128(tmp, sign);               \
    _mm_sub_epi32(tmp, sign);                     \
  })

#define round_shift_array_32_sse2(input, output, size, bit) \
  ({                                                        \
    if (bit > 0) {                                          \
      int i;                                                \
      for (i = 0; i < size; i++) {                          \
        output[i] = round_shift_32_sse2(input[i], bit);     \
      }                                                     \
    } else {                                                \
      int i;                                                \
      for (i = 0; i < size; i++) {                          \
        output[i] = _mm_slli_epi32(input[i], -bit);         \
      }                                                     \
    }                                                       \
  })

// out0 = in0*w0 + in1*w1
// out1 = -in1*w0 + in0*w1
#define btf_32_sse2_type0(w0, w1, in0, in1, out0, out1, bit) \
  ({                                                         \
    __m128i ww0, ww1, in0_w0, in1_w1, in0_w1, in1_w0;          \
    ww0 = _mm_set1_epi32(w0);                                \
    ww1 = _mm_set1_epi32(w1);                                \
    in0_w0 = mullo_epi32(in0, ww0);                          \
    in1_w1 = mullo_epi32(in1, ww1);                          \
    out0 = _mm_add_epi32(in0_w0, in1_w1);                    \
    out0 = round_shift_32_sse2(out0, bit);                   \
    in0_w1 = mullo_epi32(in0, ww1);                          \
    in1_w0 = mullo_epi32(in1, ww0);                          \
    out1 = _mm_sub_epi32(in0_w1, in1_w0);                    \
    out1 = round_shift_32_sse2(out1, bit);                   \
  })

// out0 = in0*w0 + in1*w1
// out1 = in1*w0 - in0*w1
#define btf_32_sse2_type1(w0, w1, in0, in1, out0, out1, bit) \
  ({                                                         \
    __m128i ww0, ww1, in0_w0, in1_w1, in0_w1, in1_w0;          \
    ww0 = _mm_set1_epi32(w0);                                \
    ww1 = _mm_set1_epi32(w1);                                \
    in0_w0 = mullo_epi32(in0, ww0);                          \
    in1_w1 = mullo_epi32(in1, ww1);                          \
    out0 = _mm_add_epi32(in0_w0, in1_w1);                    \
    out0 = round_shift_32_sse2(out0, bit);                   \
    in0_w1 = mullo_epi32(in0, ww1);                          \
    in1_w0 = mullo_epi32(in1, ww0);                          \
    out1 = _mm_sub_epi32(in1_w0, in0_w1);                    \
    out1 = round_shift_32_sse2(out1, bit);                   \
  })

#ifdef __cplusplus
}
#endif

#endif  // VP10_TXMF1D_SSE2_H_
