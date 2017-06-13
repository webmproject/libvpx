/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "vpx_dsp/ppc/types_vsx.h"

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/inv_txfm.h"

static int16x8_t cospi2_v = { 16305, 16305, 16305, 16305,
                              16305, 16305, 16305, 16305 };
static int16x8_t cospi4_v = { 16069, 16069, 16069, 16069,
                              16069, 16069, 16069, 16069 };
static int16x8_t cospi6_v = { 15679, 15679, 15679, 15679,
                              15679, 15679, 15679, 15679 };
static int16x8_t cospi8_v = { 15137, 15137, 15137, 15137,
                              15137, 15137, 15137, 15137 };
static int16x8_t cospi10_v = { 14449, 14449, 14449, 14449,
                               14449, 14449, 14449, 14449 };
static int16x8_t cospi12_v = { 13623, 13623, 13623, 13623,
                               13623, 13623, 13623, 13623 };
static int16x8_t cospi14_v = { 12665, 12665, 12665, 12665,
                               12665, 12665, 12665, 12665 };
static int16x8_t cospi16_v = { 11585, 11585, 11585, 11585,
                               11585, 11585, 11585, 11585 };
static int16x8_t cospi18_v = { 10394, 10394, 10394, 10394,
                               10394, 10394, 10394, 10394 };
static int16x8_t cospi20_v = { 9102, 9102, 9102, 9102, 9102, 9102, 9102, 9102 };
static int16x8_t cospi22_v = { 7723, 7723, 7723, 7723, 7723, 7723, 7723, 7723 };
static int16x8_t cospi24_v = { 6270, 6270, 6270, 6270, 6270, 6270, 6270, 6270 };
static int16x8_t cospi24_mv = { -6270, -6270, -6270, -6270,
                                -6270, -6270, -6270, -6270 };
static int16x8_t cospi26_v = { 4756, 4756, 4756, 4756, 4756, 4756, 4756, 4756 };
static int16x8_t cospi28_v = { 3196, 3196, 3196, 3196, 3196, 3196, 3196, 3196 };
static int16x8_t cospi30_v = { 1606, 1606, 1606, 1606, 1606, 1606, 1606, 1606 };

#define ROUND_SHIFT_INIT                                               \
  const int32x4_t shift = vec_sl(vec_splat_s32(1), vec_splat_u32(13)); \
  const uint32x4_t shift14 = vec_splat_u32(14);

#define DCT_CONST_ROUND_SHIFT(vec) vec = vec_sra(vec_add(vec, shift), shift14);

#define PIXEL_ADD_INIT               \
  int16x8_t add8 = vec_splat_s16(8); \
  uint16x8_t shift4 = vec_splat_u16(4);

#define PIXEL_ADD4(out, in) out = vec_sra(vec_add(in, add8), shift4);

#define IDCT4(in0, in1, out0, out1)                                           \
  t0 = vec_add(in0, in1);                                                     \
  t1 = vec_sub(in0, in1);                                                     \
  tmp16_0 = vec_mergeh(t0, t1);                                               \
  temp1 = vec_sra(vec_add(vec_mule(tmp16_0, cospi16_v), shift), shift14);     \
  temp2 = vec_sra(vec_add(vec_mulo(tmp16_0, cospi16_v), shift), shift14);     \
                                                                              \
  tmp16_0 = vec_mergel(in0, in1);                                             \
  temp3 = vec_sub(vec_mule(tmp16_0, cospi24_v), vec_mulo(tmp16_0, cospi8_v)); \
  DCT_CONST_ROUND_SHIFT(temp3);                                               \
  temp4 = vec_add(vec_mule(tmp16_0, cospi8_v), vec_mulo(tmp16_0, cospi24_v)); \
  DCT_CONST_ROUND_SHIFT(temp4);                                               \
                                                                              \
  step0 = vec_packs(temp1, temp2);                                            \
  step1 = vec_packs(temp4, temp3);                                            \
  out0 = vec_add(step0, step1);                                               \
  out1 = vec_sub(step0, step1);                                               \
  out1 = vec_perm(out1, out1, mask0);

void vpx_idct4x4_16_add_vsx(const tran_low_t *input, uint8_t *dest,
                            int stride) {
  int32x4_t temp1, temp2, temp3, temp4;
  int16x8_t step0, step1, tmp16_0, tmp16_1, t_out0, t_out1;
  uint8x16_t mask0 = { 0x8, 0x9, 0xA, 0xB, 0xC, 0xD, 0xE, 0xF,
                       0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7 };
  uint8x16_t mask1 = { 0x0,  0x1,  0x2,  0x3,  0x4,  0x5,  0x6,  0x7,
                       0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17 };
  int16x8_t v0 = vec_vsx_ld(0, input);
  int16x8_t v1 = vec_vsx_ld(16, input);
  int16x8_t t0 = vec_mergeh(v0, v1);
  int16x8_t t1 = vec_mergel(v0, v1);

  uint8x16_t dest0 = vec_vsx_ld(0, dest);
  uint8x16_t dest1 = vec_vsx_ld(stride, dest);
  uint8x16_t dest2 = vec_vsx_ld(2 * stride, dest);
  uint8x16_t dest3 = vec_vsx_ld(3 * stride, dest);
  uint8x16_t zerov = vec_splat_u8(0);
  int16x8_t d_u0 = (int16x8_t)vec_mergeh(dest0, zerov);
  int16x8_t d_u1 = (int16x8_t)vec_mergeh(dest1, zerov);
  int16x8_t d_u2 = (int16x8_t)vec_mergeh(dest2, zerov);
  int16x8_t d_u3 = (int16x8_t)vec_mergeh(dest3, zerov);
  uint8x16_t output_v;
  uint8_t tmp_dest[16];
  ROUND_SHIFT_INIT
  PIXEL_ADD_INIT;

  v0 = vec_mergeh(t0, t1);
  v1 = vec_mergel(t0, t1);

  IDCT4(v0, v1, t_out0, t_out1);
  // transpose
  t0 = vec_mergeh(t_out0, t_out1);
  t1 = vec_mergel(t_out0, t_out1);
  v0 = vec_mergeh(t0, t1);
  v1 = vec_mergel(t0, t1);
  IDCT4(v0, v1, t_out0, t_out1);

  PIXEL_ADD4(v0, t_out0);
  PIXEL_ADD4(v1, t_out1);
  tmp16_0 = vec_add(vec_perm(d_u0, d_u1, mask1), v0);
  tmp16_1 = vec_add(vec_perm(d_u2, d_u3, mask1), v1);
  output_v = vec_packsu(tmp16_0, tmp16_1);

  vec_vsx_st(output_v, 0, tmp_dest);
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) dest[j * stride + i] = tmp_dest[j * 4 + i];
}

#define TRANSPOSE8x8(in0, in1, in2, in3, in4, in5, in6, in7, out0, out1, out2, \
                     out3, out4, out5, out6, out7)                             \
  out0 = vec_mergeh(in0, in1);                                                 \
  out1 = vec_mergel(in0, in1);                                                 \
  out2 = vec_mergeh(in2, in3);                                                 \
  out3 = vec_mergel(in2, in3);                                                 \
  out4 = vec_mergeh(in4, in5);                                                 \
  out5 = vec_mergel(in4, in5);                                                 \
  out6 = vec_mergeh(in6, in7);                                                 \
  out7 = vec_mergel(in6, in7);                                                 \
  in0 = (int16x8_t)vec_mergeh((int32x4_t)out0, (int32x4_t)out2);               \
  in1 = (int16x8_t)vec_mergel((int32x4_t)out0, (int32x4_t)out2);               \
  in2 = (int16x8_t)vec_mergeh((int32x4_t)out1, (int32x4_t)out3);               \
  in3 = (int16x8_t)vec_mergel((int32x4_t)out1, (int32x4_t)out3);               \
  in4 = (int16x8_t)vec_mergeh((int32x4_t)out4, (int32x4_t)out6);               \
  in5 = (int16x8_t)vec_mergel((int32x4_t)out4, (int32x4_t)out6);               \
  in6 = (int16x8_t)vec_mergeh((int32x4_t)out5, (int32x4_t)out7);               \
  in7 = (int16x8_t)vec_mergel((int32x4_t)out5, (int32x4_t)out7);               \
  out0 = vec_perm(in0, in4, tr8_mask0);                                        \
  out1 = vec_perm(in0, in4, tr8_mask1);                                        \
  out2 = vec_perm(in1, in5, tr8_mask0);                                        \
  out3 = vec_perm(in1, in5, tr8_mask1);                                        \
  out4 = vec_perm(in2, in6, tr8_mask0);                                        \
  out5 = vec_perm(in2, in6, tr8_mask1);                                        \
  out6 = vec_perm(in3, in7, tr8_mask0);                                        \
  out7 = vec_perm(in3, in7, tr8_mask1);

#define STEP8_0(inpt0, inpt1, outpt0, outpt1, cospi0, cospi1)             \
  tmp16_0 = vec_mergeh(inpt0, inpt1);                                     \
  tmp16_1 = vec_mergel(inpt0, inpt1);                                     \
  temp10 = vec_sub(vec_mule(tmp16_0, cospi0), vec_mulo(tmp16_0, cospi1)); \
  temp11 = vec_sub(vec_mule(tmp16_1, cospi0), vec_mulo(tmp16_1, cospi1)); \
  DCT_CONST_ROUND_SHIFT(temp10);                                          \
  DCT_CONST_ROUND_SHIFT(temp11);                                          \
  outpt0 = vec_packs(temp10, temp11);                                     \
  temp10 = vec_add(vec_mule(tmp16_0, cospi1), vec_mulo(tmp16_0, cospi0)); \
  temp11 = vec_add(vec_mule(tmp16_1, cospi1), vec_mulo(tmp16_1, cospi0)); \
  DCT_CONST_ROUND_SHIFT(temp10);                                          \
  DCT_CONST_ROUND_SHIFT(temp11);                                          \
  outpt1 = vec_packs(temp10, temp11);

#define STEP8_1(inpt0, inpt1, outpt0, outpt1, cospi) \
  tmp16_2 = vec_sub(inpt0, inpt1);                   \
  tmp16_3 = vec_add(inpt0, inpt1);                   \
  tmp16_0 = vec_mergeh(tmp16_2, tmp16_3);            \
  tmp16_1 = vec_mergel(tmp16_2, tmp16_3);            \
  temp10 = vec_mule(tmp16_0, cospi);                 \
  temp11 = vec_mule(tmp16_1, cospi);                 \
  DCT_CONST_ROUND_SHIFT(temp10);                     \
  DCT_CONST_ROUND_SHIFT(temp11);                     \
  outpt0 = vec_packs(temp10, temp11);                \
  temp10 = vec_mulo(tmp16_0, cospi);                 \
  temp11 = vec_mulo(tmp16_1, cospi);                 \
  DCT_CONST_ROUND_SHIFT(temp10);                     \
  DCT_CONST_ROUND_SHIFT(temp11);                     \
  outpt1 = vec_packs(temp10, temp11);

#define IDCT8(in0, in1, in2, in3, in4, in5, in6, in7)    \
  /* stage 1 */                                          \
  step0 = in0;                                           \
  step2 = in4;                                           \
  step1 = in2;                                           \
  step3 = in6;                                           \
                                                         \
  STEP8_0(in1, in7, step4, step7, cospi28_v, cospi4_v);  \
  STEP8_0(in5, in3, step5, step6, cospi12_v, cospi20_v); \
                                                         \
  /* stage 2 */                                          \
  STEP8_1(step0, step2, in1, in0, cospi16_v);            \
  STEP8_0(step1, step3, in2, in3, cospi24_v, cospi8_v);  \
  in4 = vec_add(step4, step5);                           \
  in5 = vec_sub(step4, step5);                           \
  in6 = vec_sub(step7, step6);                           \
  in7 = vec_add(step6, step7);                           \
                                                         \
  /* stage 3 */                                          \
  step0 = vec_add(in0, in3);                             \
  step1 = vec_add(in1, in2);                             \
  step2 = vec_sub(in1, in2);                             \
  step3 = vec_sub(in0, in3);                             \
  step4 = in4;                                           \
  STEP8_1(in6, in5, step5, step6, cospi16_v);            \
  step7 = in7;                                           \
                                                         \
  /* stage 4 */                                          \
  in0 = vec_add(step0, step7);                           \
  in1 = vec_add(step1, step6);                           \
  in2 = vec_add(step2, step5);                           \
  in3 = vec_add(step3, step4);                           \
  in4 = vec_sub(step3, step4);                           \
  in5 = vec_sub(step2, step5);                           \
  in6 = vec_sub(step1, step6);                           \
  in7 = vec_sub(step0, step7);

#define PIXEL_ADD(in, out, add, shiftx) \
  out = vec_add(vec_sra(vec_add(in, add), shiftx), out);

static uint8x16_t tr8_mask0 = {
  0x0,  0x1,  0x2,  0x3,  0x4,  0x5,  0x6,  0x7,
  0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17
};
static uint8x16_t tr8_mask1 = {
  0x8,  0x9,  0xA,  0xB,  0xC,  0xD,  0xE,  0xF,
  0x18, 0x19, 0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F
};
void vpx_idct8x8_64_add_vsx(const tran_low_t *input, uint8_t *dest,
                            int stride) {
  int32x4_t temp10, temp11;
  int16x8_t step0, step1, step2, step3, step4, step5, step6, step7;
  int16x8_t tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp16_0, tmp16_1,
      tmp16_2, tmp16_3;
  int16x8_t src0 = vec_vsx_ld(0, input);
  int16x8_t src1 = vec_vsx_ld(16, input);
  int16x8_t src2 = vec_vsx_ld(2 * 16, input);
  int16x8_t src3 = vec_vsx_ld(3 * 16, input);
  int16x8_t src4 = vec_vsx_ld(4 * 16, input);
  int16x8_t src5 = vec_vsx_ld(5 * 16, input);
  int16x8_t src6 = vec_vsx_ld(6 * 16, input);
  int16x8_t src7 = vec_vsx_ld(7 * 16, input);
  uint8x16_t dest0 = vec_vsx_ld(0, dest);
  uint8x16_t dest1 = vec_vsx_ld(stride, dest);
  uint8x16_t dest2 = vec_vsx_ld(2 * stride, dest);
  uint8x16_t dest3 = vec_vsx_ld(3 * stride, dest);
  uint8x16_t dest4 = vec_vsx_ld(4 * stride, dest);
  uint8x16_t dest5 = vec_vsx_ld(5 * stride, dest);
  uint8x16_t dest6 = vec_vsx_ld(6 * stride, dest);
  uint8x16_t dest7 = vec_vsx_ld(7 * stride, dest);
  uint8x16_t zerov = vec_splat_u8(0);
  int16x8_t d_u0 = (int16x8_t)vec_mergeh(dest0, zerov);
  int16x8_t d_u1 = (int16x8_t)vec_mergeh(dest1, zerov);
  int16x8_t d_u2 = (int16x8_t)vec_mergeh(dest2, zerov);
  int16x8_t d_u3 = (int16x8_t)vec_mergeh(dest3, zerov);
  int16x8_t d_u4 = (int16x8_t)vec_mergeh(dest4, zerov);
  int16x8_t d_u5 = (int16x8_t)vec_mergeh(dest5, zerov);
  int16x8_t d_u6 = (int16x8_t)vec_mergeh(dest6, zerov);
  int16x8_t d_u7 = (int16x8_t)vec_mergeh(dest7, zerov);
  int16x8_t add = vec_sl(vec_splat_s16(8), vec_splat_u16(1));
  uint16x8_t shift5 = vec_splat_u16(5);
  uint8x16_t output0, output1, output2, output3;
  ROUND_SHIFT_INIT;

  TRANSPOSE8x8(src0, src1, src2, src3, src4, src5, src6, src7, tmp0, tmp1, tmp2,
               tmp3, tmp4, tmp5, tmp6, tmp7);

  IDCT8(tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7);
  TRANSPOSE8x8(tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, src0, src1, src2,
               src3, src4, src5, src6, src7);
  IDCT8(src0, src1, src2, src3, src4, src5, src6, src7);
  PIXEL_ADD(src0, d_u0, add, shift5);
  PIXEL_ADD(src1, d_u1, add, shift5);
  PIXEL_ADD(src2, d_u2, add, shift5);
  PIXEL_ADD(src3, d_u3, add, shift5);
  PIXEL_ADD(src4, d_u4, add, shift5);
  PIXEL_ADD(src5, d_u5, add, shift5);
  PIXEL_ADD(src6, d_u6, add, shift5);
  PIXEL_ADD(src7, d_u7, add, shift5);
  output0 = vec_packsu(d_u0, d_u1);
  output1 = vec_packsu(d_u2, d_u3);
  output2 = vec_packsu(d_u4, d_u5);
  output3 = vec_packsu(d_u6, d_u7);

  vec_vsx_st(xxpermdi(output0, dest0, 1), 0, dest);
  vec_vsx_st(xxpermdi(output0, dest1, 3), stride, dest);
  vec_vsx_st(xxpermdi(output1, dest2, 1), 2 * stride, dest);
  vec_vsx_st(xxpermdi(output1, dest3, 3), 3 * stride, dest);
  vec_vsx_st(xxpermdi(output2, dest4, 1), 4 * stride, dest);
  vec_vsx_st(xxpermdi(output2, dest5, 3), 5 * stride, dest);
  vec_vsx_st(xxpermdi(output3, dest6, 1), 6 * stride, dest);
  vec_vsx_st(xxpermdi(output3, dest7, 3), 7 * stride, dest);
}

#define LOAD_INPUT16(source, offset, step, in0, in1, in2, in3, in4, in5, in6, \
                     in7, in8, in9, inA, inB, inC, inD, inE, inF)             \
  in0 = vec_vsx_ld(offset, source);                                           \
  in1 = vec_vsx_ld(step + offset, source);                                    \
  in2 = vec_vsx_ld(2 * step + offset, source);                                \
  in3 = vec_vsx_ld(3 * step + offset, source);                                \
  in4 = vec_vsx_ld(4 * step + offset, source);                                \
  in5 = vec_vsx_ld(5 * step + offset, source);                                \
  in6 = vec_vsx_ld(6 * step + offset, source);                                \
  in7 = vec_vsx_ld(7 * step + offset, source);                                \
  in8 = vec_vsx_ld(8 * step + offset, source);                                \
  in9 = vec_vsx_ld(9 * step + offset, source);                                \
  inA = vec_vsx_ld(10 * step + offset, source);                               \
  inB = vec_vsx_ld(11 * step + offset, source);                               \
  inC = vec_vsx_ld(12 * step + offset, source);                               \
  inD = vec_vsx_ld(13 * step + offset, source);                               \
  inE = vec_vsx_ld(14 * step + offset, source);                               \
  inF = vec_vsx_ld(15 * step + offset, source);

#define STEP16_1(inpt0, inpt1, outpt0, outpt1, cospi) \
  tmp16_0 = vec_mergeh(inpt0, inpt1);                 \
  tmp16_1 = vec_mergel(inpt0, inpt1);                 \
  temp10 = vec_mule(tmp16_0, cospi);                  \
  temp11 = vec_mule(tmp16_1, cospi);                  \
  temp20 = vec_mulo(tmp16_0, cospi);                  \
  temp21 = vec_mulo(tmp16_1, cospi);                  \
  temp30 = vec_sub(temp10, temp20);                   \
  temp10 = vec_add(temp10, temp20);                   \
  temp20 = vec_sub(temp11, temp21);                   \
  temp21 = vec_add(temp11, temp21);                   \
  DCT_CONST_ROUND_SHIFT(temp30);                      \
  DCT_CONST_ROUND_SHIFT(temp20);                      \
  outpt0 = vec_packs(temp30, temp20);                 \
  DCT_CONST_ROUND_SHIFT(temp10);                      \
  DCT_CONST_ROUND_SHIFT(temp21);                      \
  outpt1 = vec_packs(temp10, temp21);

#define IDCT16(in0, in1, in2, in3, in4, in5, in6, in7, in8, in9, inA, inB,     \
               inC, inD, inE, inF, out0, out1, out2, out3, out4, out5, out6,   \
               out7, out8, out9, outA, outB, outC, outD, outE, outF)           \
  /* stage 1 */                                                                \
  /* out0 = in0; */                                                            \
  out1 = in8;                                                                  \
  out2 = in4;                                                                  \
  out3 = inC;                                                                  \
  out4 = in2;                                                                  \
  out5 = inA;                                                                  \
  out6 = in6;                                                                  \
  out7 = inE;                                                                  \
  out8 = in1;                                                                  \
  out9 = in9;                                                                  \
  outA = in5;                                                                  \
  outB = inD;                                                                  \
  outC = in3;                                                                  \
  outD = inB;                                                                  \
  outE = in7;                                                                  \
  outF = inF;                                                                  \
                                                                               \
  /* stage 2 */                                                                \
  /* in0 = out0; */                                                            \
  in1 = out1;                                                                  \
  in2 = out2;                                                                  \
  in3 = out3;                                                                  \
  in4 = out4;                                                                  \
  in5 = out5;                                                                  \
  in6 = out6;                                                                  \
  in7 = out7;                                                                  \
                                                                               \
  STEP8_0(out8, outF, in8, inF, cospi30_v, cospi2_v);                          \
  STEP8_0(out9, outE, in9, inE, cospi14_v, cospi18_v);                         \
  STEP8_0(outA, outD, inA, inD, cospi22_v, cospi10_v);                         \
  STEP8_0(outB, outC, inB, inC, cospi6_v, cospi26_v);                          \
                                                                               \
  /* stage 3 */                                                                \
  out0 = in0;                                                                  \
  out1 = in1;                                                                  \
  out2 = in2;                                                                  \
  out3 = in3;                                                                  \
                                                                               \
  STEP8_0(in4, in7, out4, out7, cospi28_v, cospi4_v);                          \
  STEP8_0(in5, in6, out5, out6, cospi12_v, cospi20_v);                         \
                                                                               \
  out8 = vec_add(in8, in9);                                                    \
  out9 = vec_sub(in8, in9);                                                    \
  outA = vec_sub(inB, inA);                                                    \
  outB = vec_add(inA, inB);                                                    \
  outC = vec_add(inC, inD);                                                    \
  outD = vec_sub(inC, inD);                                                    \
  outE = vec_sub(inF, inE);                                                    \
  outF = vec_add(inE, inF);                                                    \
                                                                               \
  /* stage 4 */                                                                \
  STEP16_1(out0, out1, in1, in0, cospi16_v);                                   \
  STEP8_0(out2, out3, in2, in3, cospi24_v, cospi8_v);                          \
  in4 = vec_add(out4, out5);                                                   \
  in5 = vec_sub(out4, out5);                                                   \
  in6 = vec_sub(out7, out6);                                                   \
  in7 = vec_add(out6, out7);                                                   \
                                                                               \
  in8 = out8;                                                                  \
  inF = outF;                                                                  \
  tmp16_0 = vec_mergeh(out9, outE);                                            \
  tmp16_1 = vec_mergel(out9, outE);                                            \
  temp10 = vec_sub(vec_mulo(tmp16_0, cospi24_v), vec_mule(tmp16_0, cospi8_v)); \
  temp11 = vec_sub(vec_mulo(tmp16_1, cospi24_v), vec_mule(tmp16_1, cospi8_v)); \
  DCT_CONST_ROUND_SHIFT(temp10);                                               \
  DCT_CONST_ROUND_SHIFT(temp11);                                               \
  in9 = vec_packs(temp10, temp11);                                             \
  temp10 = vec_add(vec_mule(tmp16_0, cospi24_v), vec_mulo(tmp16_0, cospi8_v)); \
  temp11 = vec_add(vec_mule(tmp16_1, cospi24_v), vec_mulo(tmp16_1, cospi8_v)); \
  DCT_CONST_ROUND_SHIFT(temp10);                                               \
  DCT_CONST_ROUND_SHIFT(temp11);                                               \
  inE = vec_packs(temp10, temp11);                                             \
                                                                               \
  tmp16_0 = vec_mergeh(outA, outD);                                            \
  tmp16_1 = vec_mergel(outA, outD);                                            \
  temp10 =                                                                     \
      vec_sub(vec_mule(tmp16_0, cospi24_mv), vec_mulo(tmp16_0, cospi8_v));     \
  temp11 =                                                                     \
      vec_sub(vec_mule(tmp16_1, cospi24_mv), vec_mulo(tmp16_1, cospi8_v));     \
  DCT_CONST_ROUND_SHIFT(temp10);                                               \
  DCT_CONST_ROUND_SHIFT(temp11);                                               \
  inA = vec_packs(temp10, temp11);                                             \
  temp10 = vec_sub(vec_mulo(tmp16_0, cospi24_v), vec_mule(tmp16_0, cospi8_v)); \
  temp11 = vec_sub(vec_mulo(tmp16_1, cospi24_v), vec_mule(tmp16_1, cospi8_v)); \
  DCT_CONST_ROUND_SHIFT(temp10);                                               \
  DCT_CONST_ROUND_SHIFT(temp11);                                               \
  inD = vec_packs(temp10, temp11);                                             \
                                                                               \
  inB = outB;                                                                  \
  inC = outC;                                                                  \
                                                                               \
  /* stage 5 */                                                                \
  out0 = vec_add(in0, in3);                                                    \
  out1 = vec_add(in1, in2);                                                    \
  out2 = vec_sub(in1, in2);                                                    \
  out3 = vec_sub(in0, in3);                                                    \
  out4 = in4;                                                                  \
  STEP16_1(in6, in5, out5, out6, cospi16_v);                                   \
  out7 = in7;                                                                  \
                                                                               \
  out8 = vec_add(in8, inB);                                                    \
  out9 = vec_add(in9, inA);                                                    \
  outA = vec_sub(in9, inA);                                                    \
  outB = vec_sub(in8, inB);                                                    \
  outC = vec_sub(inF, inC);                                                    \
  outD = vec_sub(inE, inD);                                                    \
  outE = vec_add(inD, inE);                                                    \
  outF = vec_add(inC, inF);                                                    \
                                                                               \
  /* stage 6 */                                                                \
  in0 = vec_add(out0, out7);                                                   \
  in1 = vec_add(out1, out6);                                                   \
  in2 = vec_add(out2, out5);                                                   \
  in3 = vec_add(out3, out4);                                                   \
  in4 = vec_sub(out3, out4);                                                   \
  in5 = vec_sub(out2, out5);                                                   \
  in6 = vec_sub(out1, out6);                                                   \
  in7 = vec_sub(out0, out7);                                                   \
  in8 = out8;                                                                  \
  in9 = out9;                                                                  \
  STEP16_1(outD, outA, inA, inD, cospi16_v);                                   \
  STEP16_1(outC, outB, inB, inC, cospi16_v);                                   \
  inE = outE;                                                                  \
  inF = outF;                                                                  \
                                                                               \
  /* stage 7 */                                                                \
  out0 = vec_add(in0, inF);                                                    \
  out1 = vec_add(in1, inE);                                                    \
  out2 = vec_add(in2, inD);                                                    \
  out3 = vec_add(in3, inC);                                                    \
  out4 = vec_add(in4, inB);                                                    \
  out5 = vec_add(in5, inA);                                                    \
  out6 = vec_add(in6, in9);                                                    \
  out7 = vec_add(in7, in8);                                                    \
  out8 = vec_sub(in7, in8);                                                    \
  out9 = vec_sub(in6, in9);                                                    \
  outA = vec_sub(in5, inA);                                                    \
  outB = vec_sub(in4, inB);                                                    \
  outC = vec_sub(in3, inC);                                                    \
  outD = vec_sub(in2, inD);                                                    \
  outE = vec_sub(in1, inE);                                                    \
  outF = vec_sub(in0, inF);

#define PIXEL_ADD_STORE16(in0, in1, dst, offset) \
  d_uh = (int16x8_t)vec_mergeh(dst, zerov);      \
  d_ul = (int16x8_t)vec_mergel(dst, zerov);      \
  PIXEL_ADD(in0, d_uh, add, shift6);             \
  PIXEL_ADD(in1, d_ul, add, shift6);             \
  vec_vsx_st(vec_packsu(d_uh, d_ul), offset, dest);

void vpx_idct16x16_256_add_vsx(const tran_low_t *input, uint8_t *dest,
                               int stride) {
  int32x4_t temp10, temp11, temp20, temp21, temp30;
  int16x8_t src00, src01, src02, src03, src04, src05, src06, src07, src10,
      src11, src12, src13, src14, src15, src16, src17;
  int16x8_t src20, src21, src22, src23, src24, src25, src26, src27, src30,
      src31, src32, src33, src34, src35, src36, src37;
  int16x8_t tmp00, tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp10,
      tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17, tmp16_0, tmp16_1;
  int16x8_t tmp20, tmp21, tmp22, tmp23, tmp24, tmp25, tmp26, tmp27, tmp30,
      tmp31, tmp32, tmp33, tmp34, tmp35, tmp36, tmp37;
  uint8x16_t dest0, dest1, dest2, dest3, dest4, dest5, dest6, dest7, dest8,
      dest9, destA, destB, destC, destD, destE, destF;
  int16x8_t d_uh, d_ul;
  int16x8_t add = vec_sl(vec_splat_s16(8), vec_splat_u16(2));
  uint16x8_t shift6 = vec_splat_u16(6);
  uint8x16_t zerov = vec_splat_u8(0);
  ROUND_SHIFT_INIT;

  // transform rows
  // load and transform the upper half of 16x16 matrix
  LOAD_INPUT16(input, 0, 16, src00, src10, src01, src11, src02, src12, src03,
               src13, src04, src14, src05, src15, src06, src16, src07, src17);
  TRANSPOSE8x8(src00, src01, src02, src03, src04, src05, src06, src07, tmp00,
               tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07);
  TRANSPOSE8x8(src10, src11, src12, src13, src14, src15, src16, src17, tmp10,
               tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17);
  IDCT16(tmp00, tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp10, tmp11,
         tmp12, tmp13, tmp14, tmp15, tmp16, tmp17, src00, src01, src02, src03,
         src04, src05, src06, src07, src10, src11, src12, src13, src14, src15,
         src16, src17);
  TRANSPOSE8x8(src00, src01, src02, src03, src04, src05, src06, src07, tmp00,
               tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07);
  TRANSPOSE8x8(src10, src11, src12, src13, src14, src15, src16, src17, tmp10,
               tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17);

  // load and transform the lower half of 16x16 matrix
  LOAD_INPUT16(input, 16 * 16, 16, src20, src30, src21, src31, src22, src32,
               src23, src33, src24, src34, src25, src35, src26, src36, src27,
               src37);
  TRANSPOSE8x8(src20, src21, src22, src23, src24, src25, src26, src27, tmp20,
               tmp21, tmp22, tmp23, tmp24, tmp25, tmp26, tmp27);
  TRANSPOSE8x8(src30, src31, src32, src33, src34, src35, src36, src37, tmp30,
               tmp31, tmp32, tmp33, tmp34, tmp35, tmp36, tmp37);
  IDCT16(tmp20, tmp21, tmp22, tmp23, tmp24, tmp25, tmp26, tmp27, tmp30, tmp31,
         tmp32, tmp33, tmp34, tmp35, tmp36, tmp37, src20, src21, src22, src23,
         src24, src25, src26, src27, src30, src31, src32, src33, src34, src35,
         src36, src37);
  TRANSPOSE8x8(src20, src21, src22, src23, src24, src25, src26, src27, tmp20,
               tmp21, tmp22, tmp23, tmp24, tmp25, tmp26, tmp27);
  TRANSPOSE8x8(src30, src31, src32, src33, src34, src35, src36, src37, tmp30,
               tmp31, tmp32, tmp33, tmp34, tmp35, tmp36, tmp37);

  // transform columns
  // left half first
  IDCT16(tmp00, tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp20, tmp21,
         tmp22, tmp23, tmp24, tmp25, tmp26, tmp27, src00, src01, src02, src03,
         src04, src05, src06, src07, src20, src21, src22, src23, src24, src25,
         src26, src27);
  // right half
  IDCT16(tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17, tmp30, tmp31,
         tmp32, tmp33, tmp34, tmp35, tmp36, tmp37, src10, src11, src12, src13,
         src14, src15, src16, src17, src30, src31, src32, src33, src34, src35,
         src36, src37);

  // load dest
  LOAD_INPUT16(dest, 0, stride, dest0, dest1, dest2, dest3, dest4, dest5, dest6,
               dest7, dest8, dest9, destA, destB, destC, destD, destE, destF);

  PIXEL_ADD_STORE16(src00, src10, dest0, 0);
  PIXEL_ADD_STORE16(src01, src11, dest1, stride);
  PIXEL_ADD_STORE16(src02, src12, dest2, 2 * stride);
  PIXEL_ADD_STORE16(src03, src13, dest3, 3 * stride);
  PIXEL_ADD_STORE16(src04, src14, dest4, 4 * stride);
  PIXEL_ADD_STORE16(src05, src15, dest5, 5 * stride);
  PIXEL_ADD_STORE16(src06, src16, dest6, 6 * stride);
  PIXEL_ADD_STORE16(src07, src17, dest7, 7 * stride);

  PIXEL_ADD_STORE16(src20, src30, dest8, 8 * stride);
  PIXEL_ADD_STORE16(src21, src31, dest9, 9 * stride);
  PIXEL_ADD_STORE16(src22, src32, destA, 10 * stride);
  PIXEL_ADD_STORE16(src23, src33, destB, 11 * stride);
  PIXEL_ADD_STORE16(src24, src34, destC, 12 * stride);
  PIXEL_ADD_STORE16(src25, src35, destD, 13 * stride);
  PIXEL_ADD_STORE16(src26, src36, destE, 14 * stride);
  PIXEL_ADD_STORE16(src27, src37, destF, 15 * stride);
}
