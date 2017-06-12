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

static int16x8_t cospi4_v = { 16069, 16069, 16069, 16069,
                              16069, 16069, 16069, 16069 };
static int16x8_t cospi8_v = { 15137, 15137, 15137, 15137,
                              15137, 15137, 15137, 15137 };
static int16x8_t cospi12_v = { 13623, 13623, 13623, 13623,
                               13623, 13623, 13623, 13623 };
static int16x8_t cospi16_v = { 11585, 11585, 11585, 11585,
                               11585, 11585, 11585, 11585 };
static int16x8_t cospi20_v = { 9102, 9102, 9102, 9102, 9102, 9102, 9102, 9102 };
static int16x8_t cospi24_v = { 6270, 6270, 6270, 6270, 6270, 6270, 6270, 6270 };
static int16x8_t cospi28_v = { 3196, 3196, 3196, 3196, 3196, 3196, 3196, 3196 };

#define ROUND_SHIFT_INIT                                               \
  const int32x4_t shift = vec_sl(vec_splat_s32(1), vec_splat_u32(13)); \
  const uint32x4_t shift14 = vec_splat_u32(14);

#define DCT_CONST_ROUND_SHIFT(vec) vec = vec_sra(vec_add(vec, shift), shift14);

#define PIXEL_ADD_INIT               \
  int16x8_t add8 = vec_splat_s16(8); \
  uint16x8_t shift4 = vec_splat_u16(4);

#define PIXEL_ADD(out, in) out = vec_sra(vec_add(in, add8), shift4);

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

  PIXEL_ADD(v0, t_out0);
  PIXEL_ADD(v1, t_out1);
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

#define PIXEL_ADD8(in, out) \
  out = vec_add(vec_sra(vec_add(in, add), shift5), out);

void vpx_idct8x8_64_add_vsx(const tran_low_t *input, uint8_t *dest,
                            int stride) {
  int32x4_t temp10, temp11;
  int16x8_t step0, step1, step2, step3, step4, step5, step6, step7;
  uint8x16_t tr8_mask0 = { 0x0,  0x1,  0x2,  0x3,  0x4,  0x5,  0x6,  0x7,
                           0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17 };
  uint8x16_t tr8_mask1 = { 0x8,  0x9,  0xA,  0xB,  0xC,  0xD,  0xE,  0xF,
                           0x18, 0x19, 0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F };
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
  PIXEL_ADD8(src0, d_u0);
  PIXEL_ADD8(src1, d_u1);
  PIXEL_ADD8(src2, d_u2);
  PIXEL_ADD8(src3, d_u3);
  PIXEL_ADD8(src4, d_u4);
  PIXEL_ADD8(src5, d_u5);
  PIXEL_ADD8(src6, d_u6);
  PIXEL_ADD8(src7, d_u7);
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
