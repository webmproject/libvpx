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

static int16x8_t cospi8_v = { 15137, 15137, 15137, 15137,
                              15137, 15137, 15137, 15137 };
static int16x8_t cospi16_v = { 11585, 11585, 11585, 11585,
                               11585, 11585, 11585, 11585 };
static int16x8_t cospi24_v = { 6270, 6270, 6270, 6270, 6270, 6270, 6270, 6270 };

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
