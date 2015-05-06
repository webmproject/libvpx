/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>

#include "vpx_ports/mem.h"
#include "vp9/common/vp9_idct.h"
#include "vp9/common/mips/msa/vp9_macros_msa.h"

#define SET_COSPI_PAIR(c0_h, c1_h) ({  \
  v8i16 out0, r0_m, r1_m;              \
                                       \
  r0_m = __msa_fill_h(c0_h);           \
  r1_m = __msa_fill_h(c1_h);           \
  out0 = __msa_ilvev_h(r1_m, r0_m);    \
                                       \
  out0;                                \
})

#define DOTP_CONST_PAIR(reg0, reg1, const0, const1, out0, out1) {  \
  v8i16 k0_m = __msa_fill_h(const0);                               \
  v8i16 s0_m, s1_m, s2_m, s3_m;                                    \
                                                                   \
  s0_m = __msa_fill_h(const1);                                     \
  k0_m = __msa_ilvev_h(s0_m, k0_m);                                \
                                                                   \
  s0_m = __msa_ilvl_h(-reg1, reg0);                                \
  s1_m = __msa_ilvr_h(-reg1, reg0);                                \
  s2_m = __msa_ilvl_h(reg0, reg1);                                 \
  s3_m = __msa_ilvr_h(reg0, reg1);                                 \
  s1_m = (v8i16)__msa_dotp_s_w(s1_m, k0_m);                        \
  s0_m = (v8i16)__msa_dotp_s_w(s0_m, k0_m);                        \
  s1_m = (v8i16)__msa_srari_w((v4i32)s1_m, DCT_CONST_BITS);        \
  s0_m = (v8i16)__msa_srari_w((v4i32)s0_m, DCT_CONST_BITS);        \
  out0 = __msa_pckev_h(s0_m, s1_m);                                \
                                                                   \
  s1_m = (v8i16)__msa_dotp_s_w(s3_m, k0_m);                        \
  s0_m = (v8i16)__msa_dotp_s_w(s2_m, k0_m);                        \
  s1_m = (v8i16)__msa_srari_w((v4i32)s1_m, DCT_CONST_BITS);        \
  s0_m = (v8i16)__msa_srari_w((v4i32)s0_m, DCT_CONST_BITS);        \
  out1 = __msa_pckev_h(s0_m, s1_m);                                \
}

#define VP9_MADD_SHORT(m0, m1, c0, c1, res0, res1) {      \
  v4i32 madd0_m, madd1_m, madd2_m, madd3_m;               \
  v8i16 madd_s0_m, madd_s1_m;                             \
                                                          \
  ILV_H_LR_SH(m0, m1, madd_s1_m, madd_s0_m);              \
                                                          \
  DOTP_S_W_4VECS_SW(madd_s0_m, c0, madd_s1_m, c0,         \
                    madd_s0_m, c1, madd_s1_m, c1,         \
                    madd0_m, madd1_m, madd2_m, madd3_m);  \
                                                          \
  SRARI_W_4VECS_SW(madd0_m, madd1_m, madd2_m, madd3_m,    \
                   madd0_m, madd1_m, madd2_m, madd3_m,    \
                   DCT_CONST_BITS);                       \
                                                          \
  PCKEV_H_2VECS_SH(madd1_m, madd0_m, madd3_m, madd2_m,    \
                   res0, res1);                           \
}

#define VP9_MADD_BF(inp0, inp1, inp2, inp3,                   \
                    cst0, cst1, cst2, cst3,                   \
                    out0, out1, out2, out3) {                 \
  v8i16 madd_s0_m, madd_s1_m, madd_s2_m, madd_s3_m;           \
  v4i32 tmp0_m, tmp1_m, tmp2_m, tmp3_m;                       \
  v4i32 m4_m, m5_m;                                           \
                                                              \
  ILV_H_LRLR_SH(inp0, inp1, inp2, inp3,                       \
                madd_s1_m, madd_s0_m, madd_s3_m, madd_s2_m);  \
                                                              \
  DOTP_S_W_4VECS_SW(madd_s0_m, cst0, madd_s1_m, cst0,         \
                    madd_s2_m, cst2, madd_s3_m, cst2,         \
                    tmp0_m, tmp1_m, tmp2_m, tmp3_m);          \
                                                              \
  m4_m = tmp0_m + tmp2_m;                                     \
  m5_m = tmp1_m + tmp3_m;                                     \
  tmp3_m = tmp1_m - tmp3_m;                                   \
  tmp2_m = tmp0_m - tmp2_m;                                   \
                                                              \
  SRARI_W_4VECS_SW(m4_m, m5_m, tmp2_m, tmp3_m,                \
                   m4_m, m5_m, tmp2_m, tmp3_m,                \
                   DCT_CONST_BITS);                           \
                                                              \
  PCKEV_H_2VECS_SH(m5_m, m4_m, tmp3_m, tmp2_m, out0, out1);   \
                                                              \
  DOTP_S_W_4VECS_SW(madd_s0_m, cst1, madd_s1_m, cst1,         \
                    madd_s2_m, cst3, madd_s3_m, cst3,         \
                    tmp0_m, tmp1_m, tmp2_m, tmp3_m);          \
                                                              \
  m4_m = tmp0_m + tmp2_m;                                     \
  m5_m = tmp1_m + tmp3_m;                                     \
  tmp3_m = tmp1_m - tmp3_m;                                   \
  tmp2_m = tmp0_m - tmp2_m;                                   \
                                                              \
  SRARI_W_4VECS_SW(m4_m, m5_m, tmp2_m, tmp3_m,                \
                   m4_m, m5_m, tmp2_m, tmp3_m,                \
                   DCT_CONST_BITS);                           \
                                                              \
  PCKEV_H_2VECS_SH(m5_m, m4_m, tmp3_m, tmp2_m, out2, out3);   \
}

#define TRANSPOSE8x8_H1(in0, in1, in2, in3,                   \
                        in4, in5, in6, in7,                   \
                        out0, out1, out2, out3,               \
                        out4, out5, out6, out7) {             \
  v8i16 loc0_m, loc1_m;                                       \
  v8i16 tmp0_m, tmp1_m, tmp2_m, tmp3_m;                       \
  v8i16 tmp4_m, tmp5_m, tmp6_m, tmp7_m;                       \
                                                              \
  loc0_m = __msa_ilvr_h((in6), (in4));                        \
  loc1_m = __msa_ilvr_h((in7), (in5));                        \
  tmp0_m = __msa_ilvr_h(loc1_m, loc0_m);                      \
  tmp1_m = __msa_ilvl_h(loc1_m, loc0_m);                      \
                                                              \
  loc0_m = __msa_ilvl_h((in6), (in4));                        \
  loc1_m = __msa_ilvl_h((in7), (in5));                        \
  tmp2_m = __msa_ilvr_h(loc1_m, loc0_m);                      \
  tmp3_m = __msa_ilvl_h(loc1_m, loc0_m);                      \
                                                              \
  loc0_m = __msa_ilvr_h((in2), (in0));                        \
  loc1_m = __msa_ilvr_h((in3), (in1));                        \
  tmp4_m = __msa_ilvr_h(loc1_m, loc0_m);                      \
  tmp5_m = __msa_ilvl_h(loc1_m, loc0_m);                      \
                                                              \
  loc0_m = __msa_ilvl_h((in2), (in0));                        \
  loc1_m = __msa_ilvl_h((in3), (in1));                        \
  tmp6_m = __msa_ilvr_h(loc1_m, loc0_m);                      \
  tmp7_m = __msa_ilvl_h(loc1_m, loc0_m);                      \
                                                              \
  out0 = (v8i16)__msa_pckev_d((v2i64)tmp0_m, (v2i64)tmp4_m);  \
  out1 = (v8i16)__msa_pckod_d((v2i64)tmp0_m, (v2i64)tmp4_m);  \
  out2 = (v8i16)__msa_pckev_d((v2i64)tmp1_m, (v2i64)tmp5_m);  \
  out3 = (v8i16)__msa_pckod_d((v2i64)tmp1_m, (v2i64)tmp5_m);  \
  out4 = (v8i16)__msa_pckev_d((v2i64)tmp2_m, (v2i64)tmp6_m);  \
  out5 = (v8i16)__msa_pckod_d((v2i64)tmp2_m, (v2i64)tmp6_m);  \
  out6 = (v8i16)__msa_pckev_d((v2i64)tmp3_m, (v2i64)tmp7_m);  \
  out7 = (v8i16)__msa_pckod_d((v2i64)tmp3_m, (v2i64)tmp7_m);  \
}

#define VP9_IADST8x16_1D(r0, r1, r2, r3, r4, r5, r6, r7,                  \
                         r8, r9, r10, r11, r12, r13, r14, r15,            \
                         out0, out1, out2, out3, out4, out5, out6, out7,  \
                         out8, out9, out10, out11,                        \
                         out12, out13, out14, out15) {                    \
  v8i16 g0_m, g1_m, g2_m, g3_m, g4_m, g5_m, g6_m, g7_m;                   \
  v8i16 g8_m, g9_m, g10_m, g11_m, g12_m, g13_m, g14_m, g15_m;             \
  v8i16 h0_m, h1_m, h2_m, h3_m, h4_m, h5_m, h6_m, h7_m;                   \
  v8i16 h8_m, h9_m, h10_m, h11_m;                                         \
  v8i16 k0_m, k1_m, k2_m, k3_m;                                           \
                                                                          \
  /* stage 1 */                                                           \
  k0_m = SET_COSPI_PAIR(cospi_1_64, cospi_31_64);                         \
  k1_m = SET_COSPI_PAIR(cospi_31_64, -cospi_1_64);                        \
  k2_m = SET_COSPI_PAIR(cospi_17_64, cospi_15_64);                        \
  k3_m = SET_COSPI_PAIR(cospi_15_64, -cospi_17_64);                       \
  VP9_MADD_BF(r15, r0, r7, r8, k0_m, k1_m, k2_m, k3_m,                    \
              g0_m, g1_m, g2_m, g3_m);                                    \
                                                                          \
  k0_m = SET_COSPI_PAIR(cospi_5_64, cospi_27_64);                         \
  k1_m = SET_COSPI_PAIR(cospi_27_64, -cospi_5_64);                        \
  k2_m = SET_COSPI_PAIR(cospi_21_64, cospi_11_64);                        \
  k3_m = SET_COSPI_PAIR(cospi_11_64, -cospi_21_64);                       \
  VP9_MADD_BF(r13, r2, r5, r10, k0_m, k1_m, k2_m, k3_m,                   \
              g4_m, g5_m, g6_m, g7_m);                                    \
                                                                          \
  k0_m = SET_COSPI_PAIR(cospi_9_64, cospi_23_64);                         \
  k1_m = SET_COSPI_PAIR(cospi_23_64, -cospi_9_64);                        \
  k2_m = SET_COSPI_PAIR(cospi_25_64, cospi_7_64);                         \
  k3_m = SET_COSPI_PAIR(cospi_7_64, -cospi_25_64);                        \
  VP9_MADD_BF(r11, r4, r3, r12, k0_m, k1_m, k2_m, k3_m,                   \
              g8_m, g9_m, g10_m, g11_m);                                  \
                                                                          \
  k0_m = SET_COSPI_PAIR(cospi_13_64, cospi_19_64);                        \
  k1_m = SET_COSPI_PAIR(cospi_19_64, -cospi_13_64);                       \
  k2_m = SET_COSPI_PAIR(cospi_29_64, cospi_3_64);                         \
  k3_m = SET_COSPI_PAIR(cospi_3_64, -cospi_29_64);                        \
  VP9_MADD_BF(r9, r6, r1, r14, k0_m, k1_m, k2_m, k3_m,                    \
              g12_m, g13_m, g14_m, g15_m);                                \
                                                                          \
  /* stage 2 */                                                           \
  k0_m = SET_COSPI_PAIR(cospi_4_64, cospi_28_64);                         \
  k1_m = SET_COSPI_PAIR(cospi_28_64, -cospi_4_64);                        \
  k2_m = SET_COSPI_PAIR(-cospi_28_64, cospi_4_64);                        \
  VP9_MADD_BF(g1_m, g3_m, g9_m, g11_m, k0_m, k1_m, k2_m, k0_m,            \
              h0_m, h1_m, h2_m, h3_m);                                    \
                                                                          \
  k0_m = SET_COSPI_PAIR(cospi_12_64, cospi_20_64);                        \
  k1_m = SET_COSPI_PAIR(-cospi_20_64, cospi_12_64);                       \
  k2_m = SET_COSPI_PAIR(cospi_20_64, -cospi_12_64);                       \
  VP9_MADD_BF(g7_m, g5_m, g15_m, g13_m, k0_m, k1_m, k2_m, k0_m,           \
              h4_m, h5_m, h6_m, h7_m);                                    \
                                                                          \
  BUTTERFLY_4(h0_m, h2_m, h6_m, h4_m, out8, out9, out11, out10);          \
                                                                          \
  BUTTERFLY_8(g0_m, g2_m, g4_m, g6_m, g14_m, g12_m, g10_m, g8_m,          \
              h8_m, h9_m, h10_m, h11_m, h6_m, h4_m, h2_m, h0_m);          \
                                                                          \
  /* stage 3 */                                                           \
  BUTTERFLY_4(h8_m, h9_m, h11_m, h10_m, out0, out1, h11_m, h10_m);        \
                                                                          \
  k0_m = SET_COSPI_PAIR(cospi_8_64, cospi_24_64);                         \
  k1_m = SET_COSPI_PAIR(cospi_24_64, -cospi_8_64);                        \
  k2_m = SET_COSPI_PAIR(-cospi_24_64, cospi_8_64);                        \
  VP9_MADD_BF(h0_m, h2_m, h4_m, h6_m, k0_m, k1_m, k2_m, k0_m,             \
              out4, out6, out5, out7);                                    \
  VP9_MADD_BF(h1_m, h3_m, h5_m, h7_m, k0_m, k1_m, k2_m, k0_m,             \
              out12, out14, out13, out15);                                \
                                                                          \
  /* stage 4 */                                                           \
  k0_m = SET_COSPI_PAIR(cospi_16_64, cospi_16_64);                        \
  k1_m = SET_COSPI_PAIR(-cospi_16_64, -cospi_16_64);                      \
  k2_m = SET_COSPI_PAIR(cospi_16_64, -cospi_16_64);                       \
  k3_m = SET_COSPI_PAIR(-cospi_16_64, cospi_16_64);                       \
  VP9_MADD_SHORT(h10_m, h11_m, k1_m, k2_m, out2, out3);                   \
  VP9_MADD_SHORT(out6, out7, k0_m, k3_m, out6, out7);                     \
  VP9_MADD_SHORT(out10, out11, k0_m, k3_m, out10, out11);                 \
  VP9_MADD_SHORT(out14, out15, k1_m, k2_m, out14, out15);                 \
}

#define VP9_ADDBLK_CLIP_AND_STORE_8_BYTES_4(dest, dest_stride,     \
                                            in0, in1, in2, in3) {  \
  uint64_t out0_m, out1_m, out2_m, out3_m;                         \
  v8i16 res0_m, res1_m, res2_m, res3_m;                            \
  v16u8 dest0_m, dest1_m, dest2_m, dest3_m;                        \
  v16i8 tmp0_m, tmp1_m;                                            \
  v16i8 zero_m = { 0 };                                            \
  uint8_t *dst_m = (uint8_t *)(dest);                              \
                                                                   \
  LOAD_4VECS_UB(dst_m, (dest_stride),                              \
                dest0_m, dest1_m, dest2_m, dest3_m);               \
                                                                   \
  res0_m = (v8i16)__msa_ilvr_b(zero_m, (v16i8)dest0_m);            \
  res1_m = (v8i16)__msa_ilvr_b(zero_m, (v16i8)dest1_m);            \
  res2_m = (v8i16)__msa_ilvr_b(zero_m, (v16i8)dest2_m);            \
  res3_m = (v8i16)__msa_ilvr_b(zero_m, (v16i8)dest3_m);            \
                                                                   \
  res0_m += (v8i16)(in0);                                          \
  res1_m += (v8i16)(in1);                                          \
  res2_m += (v8i16)(in2);                                          \
  res3_m += (v8i16)(in3);                                          \
                                                                   \
  res0_m = CLIP_UNSIGNED_CHAR_H(res0_m);                           \
  res1_m = CLIP_UNSIGNED_CHAR_H(res1_m);                           \
  res2_m = CLIP_UNSIGNED_CHAR_H(res2_m);                           \
  res3_m = CLIP_UNSIGNED_CHAR_H(res3_m);                           \
                                                                   \
  tmp0_m = __msa_pckev_b((v16i8)res1_m, (v16i8)res0_m);            \
  tmp1_m = __msa_pckev_b((v16i8)res3_m, (v16i8)res2_m);            \
                                                                   \
  out0_m = __msa_copy_u_d((v2i64)tmp0_m, 0);                       \
  out1_m = __msa_copy_u_d((v2i64)tmp0_m, 1);                       \
  out2_m = __msa_copy_u_d((v2i64)tmp1_m, 0);                       \
  out3_m = __msa_copy_u_d((v2i64)tmp1_m, 1);                       \
                                                                   \
  STORE_DWORD(dst_m, out0_m);                                      \
  dst_m += (dest_stride);                                          \
  STORE_DWORD(dst_m, out1_m);                                      \
  dst_m += (dest_stride);                                          \
  STORE_DWORD(dst_m, out2_m);                                      \
  dst_m += (dest_stride);                                          \
  STORE_DWORD(dst_m, out3_m);                                      \
}

void vp9_idct16_1d_rows_msa(const int16_t *input, int16_t *output) {
  v8i16 loc0, loc1, loc2, loc3;
  v8i16 reg0, reg2, reg4, reg6, reg8, reg10, reg12, reg14;
  v8i16 reg3, reg13, reg11, reg5, reg7, reg9, reg1, reg15;
  v8i16 tmp5, tmp6, tmp7;

  /* load left top 8x8 */
  LOAD_8VECS_SH(input, 16, reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7);

  /* load right top 8x8 */
  LOAD_8VECS_SH((input + 8), 16,
                reg8, reg9, reg10, reg11, reg12, reg13, reg14, reg15);

  /* transpose block */
  TRANSPOSE8x8_H1(reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7,
                  reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7);

  /* transpose block */
  TRANSPOSE8x8_H1(reg8, reg9, reg10, reg11, reg12, reg13, reg14, reg15,
                  reg8, reg9, reg10, reg11, reg12, reg13, reg14, reg15);

  DOTP_CONST_PAIR(reg2, reg14, cospi_28_64, cospi_4_64, reg2, reg14);
  DOTP_CONST_PAIR(reg10, reg6, cospi_12_64, cospi_20_64, reg10, reg6);

  loc0 = reg2 + reg10;
  reg2 = reg2 - reg10;
  loc1 = reg14 + reg6;
  reg14 = reg14 - reg6;

  DOTP_CONST_PAIR(reg14, reg2, cospi_16_64, cospi_16_64, loc2, loc3);
  DOTP_CONST_PAIR(reg0, reg8, cospi_16_64, cospi_16_64, reg0, reg8);
  DOTP_CONST_PAIR(reg4, reg12, cospi_24_64, cospi_8_64, reg4, reg12);

  reg14 = reg8 - reg12;
  reg2 = reg8 + reg12;
  reg10 = reg0 - reg4;
  reg6 = reg0 + reg4;

  reg0 = reg2 - loc1;
  reg2 = reg2 + loc1;
  reg12 = reg14 - loc0;
  reg14 = reg14 + loc0;
  reg4 = reg6 - loc3;
  reg6 = reg6 + loc3;
  reg8 = reg10 - loc2;
  reg10 = reg10 + loc2;

  /* stage 2 */
  DOTP_CONST_PAIR(reg1, reg15, cospi_30_64, cospi_2_64, reg1, reg15);
  DOTP_CONST_PAIR(reg9, reg7, cospi_14_64, cospi_18_64, loc2, loc3);

  reg9 = reg1 - loc2;
  reg1 = reg1 + loc2;
  reg7 = reg15 - loc3;
  reg15 = reg15 + loc3;

  DOTP_CONST_PAIR(reg5, reg11, cospi_22_64, cospi_10_64, reg5, reg11);
  DOTP_CONST_PAIR(reg13, reg3, cospi_6_64, cospi_26_64, loc0, loc1);

  reg13 = loc0 + reg5;
  reg5 = loc0 - reg5;
  reg3 = loc1 + reg11;
  reg11 = loc1 - reg11;

  loc1 = reg15 + reg3;
  reg3 = reg15 - reg3;
  loc2 = reg2 + loc1;
  reg15 = reg2 - loc1;

  loc1 = reg1 + reg13;
  reg13 = reg1 - reg13;
  loc0 = reg0 + loc1;
  loc1 = reg0 - loc1;
  tmp6 = loc0;
  tmp7 = loc1;
  reg0 = loc2;

  DOTP_CONST_PAIR(reg7, reg9, cospi_24_64, cospi_8_64, reg7, reg9);
  DOTP_CONST_PAIR((-reg5), (-reg11), cospi_8_64, cospi_24_64, reg5, reg11);

  loc0 = reg9 + reg5;
  reg5 = reg9 - reg5;
  reg2 = reg6 + loc0;
  reg1 = reg6 - loc0;

  loc0 = reg7 + reg11;
  reg11 = reg7 - reg11;
  loc1 = reg4 + loc0;
  loc2 = reg4 - loc0;
  tmp5 = loc1;

  DOTP_CONST_PAIR(reg5, reg11, cospi_16_64, cospi_16_64, reg5, reg11);

  loc0 = reg8 + reg5;
  loc1 = reg8 - reg5;
  reg4 = reg10 + reg11;
  reg9 = reg10 - reg11;
  reg10 = loc0;
  reg11 = loc1;

  DOTP_CONST_PAIR(reg3, reg13, cospi_16_64, cospi_16_64, reg3, reg13);

  reg8 = reg12 + reg3;
  reg5 = reg12 - reg3;
  reg6 = reg14 + reg13;
  reg7 = reg14 - reg13;
  reg13 = loc2;

  /* Transpose and store the output */
  reg12 = tmp5;
  reg14 = tmp6;
  reg3 = tmp7;

  /* transpose block */
  TRANSPOSE8x8_H1(reg0, reg2, reg4, reg6, reg8, reg10, reg12, reg14,
                  reg0, reg2, reg4, reg6, reg8, reg10, reg12, reg14);

  STORE_8VECS_SH(output, 16, reg0, reg2, reg4, reg6, reg8, reg10, reg12, reg14);

  /* transpose block */
  TRANSPOSE8x8_H1(reg3, reg13, reg11, reg5, reg7, reg9, reg1, reg15,
                  reg3, reg13, reg11, reg5, reg7, reg9, reg1, reg15);

  STORE_8VECS_SH((output + 8), 16,
                 reg3, reg13, reg11, reg5, reg7, reg9, reg1, reg15);
}

void vp9_idct16_1d_columns_addblk_msa(int16_t *input, uint8_t *dest,
                                      int32_t dest_stride) {
  v8i16 loc0, loc1, loc2, loc3;
  v8i16 reg0, reg2, reg4, reg6, reg8, reg10, reg12, reg14;
  v8i16 reg3, reg13, reg11, reg5, reg7, reg9, reg1, reg15;
  v8i16 tmp5, tmp6, tmp7;

  /* load up 8x8 */
  LOAD_8VECS_SH(input, 16, reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7);

  /* load bottom 8x8 */
  LOAD_8VECS_SH((input + 8 * 16), 16,
                reg8, reg9, reg10, reg11, reg12, reg13, reg14, reg15);

  DOTP_CONST_PAIR(reg2, reg14, cospi_28_64, cospi_4_64, reg2, reg14);
  DOTP_CONST_PAIR(reg10, reg6, cospi_12_64, cospi_20_64, reg10, reg6);

  loc0 = reg2 + reg10;
  reg2 = reg2 - reg10;
  loc1 = reg14 + reg6;
  reg14 = reg14 - reg6;

  DOTP_CONST_PAIR(reg14, reg2, cospi_16_64, cospi_16_64, loc2, loc3);
  DOTP_CONST_PAIR(reg0, reg8, cospi_16_64, cospi_16_64, reg0, reg8);
  DOTP_CONST_PAIR(reg4, reg12, cospi_24_64, cospi_8_64, reg4, reg12);

  reg14 = reg8 - reg12;
  reg2 = reg8 + reg12;
  reg10 = reg0 - reg4;
  reg6 = reg0 + reg4;

  reg0 = reg2 - loc1;
  reg2 = reg2 + loc1;
  reg12 = reg14 - loc0;
  reg14 = reg14 + loc0;
  reg4 = reg6 - loc3;
  reg6 = reg6 + loc3;
  reg8 = reg10 - loc2;
  reg10 = reg10 + loc2;

  /* stage 2 */
  DOTP_CONST_PAIR(reg1, reg15, cospi_30_64, cospi_2_64, reg1, reg15);
  DOTP_CONST_PAIR(reg9, reg7, cospi_14_64, cospi_18_64, loc2, loc3);

  reg9 = reg1 - loc2;
  reg1 = reg1 + loc2;
  reg7 = reg15 - loc3;
  reg15 = reg15 + loc3;

  DOTP_CONST_PAIR(reg5, reg11, cospi_22_64, cospi_10_64, reg5, reg11);
  DOTP_CONST_PAIR(reg13, reg3, cospi_6_64, cospi_26_64, loc0, loc1);

  reg13 = loc0 + reg5;
  reg5 = loc0 - reg5;
  reg3 = loc1 + reg11;
  reg11 = loc1 - reg11;

  loc1 = reg15 + reg3;
  reg3 = reg15 - reg3;
  loc2 = reg2 + loc1;
  reg15 = reg2 - loc1;

  loc1 = reg1 + reg13;
  reg13 = reg1 - reg13;
  loc0 = reg0 + loc1;
  loc1 = reg0 - loc1;
  tmp6 = loc0;
  tmp7 = loc1;
  reg0 = loc2;

  DOTP_CONST_PAIR(reg7, reg9, cospi_24_64, cospi_8_64, reg7, reg9);
  DOTP_CONST_PAIR((-reg5), (-reg11), cospi_8_64, cospi_24_64, reg5, reg11);

  loc0 = reg9 + reg5;
  reg5 = reg9 - reg5;
  reg2 = reg6 + loc0;
  reg1 = reg6 - loc0;

  loc0 = reg7 + reg11;
  reg11 = reg7 - reg11;
  loc1 = reg4 + loc0;
  loc2 = reg4 - loc0;
  tmp5 = loc1;

  DOTP_CONST_PAIR(reg5, reg11, cospi_16_64, cospi_16_64, reg5, reg11);

  loc0 = reg8 + reg5;
  loc1 = reg8 - reg5;
  reg4 = reg10 + reg11;
  reg9 = reg10 - reg11;
  reg10 = loc0;
  reg11 = loc1;

  DOTP_CONST_PAIR(reg3, reg13, cospi_16_64, cospi_16_64, reg3, reg13);

  reg8 = reg12 + reg3;
  reg5 = reg12 - reg3;
  reg6 = reg14 + reg13;
  reg7 = reg14 - reg13;
  reg13 = loc2;

  /* Transpose and store the output */
  reg12 = tmp5;
  reg14 = tmp6;
  reg3 = tmp7;

  SRARI_H_4VECS_SH(reg0, reg2, reg4, reg6, reg0, reg2, reg4, reg6, 6);
  VP9_ADDBLK_CLIP_AND_STORE_8_BYTES_4(dest, dest_stride,
                                      reg0, reg2, reg4, reg6);
  SRARI_H_4VECS_SH(reg8, reg10, reg12, reg14, reg8, reg10, reg12, reg14, 6);
  VP9_ADDBLK_CLIP_AND_STORE_8_BYTES_4((dest + (4 * dest_stride)),
                                      dest_stride, reg8, reg10, reg12, reg14);
  SRARI_H_4VECS_SH(reg3, reg13, reg11, reg5, reg3, reg13, reg11, reg5, 6);
  VP9_ADDBLK_CLIP_AND_STORE_8_BYTES_4((dest + (8 * dest_stride)),
                                      dest_stride, reg3, reg13, reg11, reg5);
  SRARI_H_4VECS_SH(reg7, reg9, reg1, reg15, reg7, reg9, reg1, reg15, 6);
  VP9_ADDBLK_CLIP_AND_STORE_8_BYTES_4((dest + (12 * dest_stride)),
                                      dest_stride, reg7, reg9, reg1, reg15);
}

void vp9_idct16x16_256_add_msa(const int16_t *input, uint8_t *dest,
                               int32_t dest_stride) {
  int32_t i;
  DECLARE_ALIGNED(32, int16_t, out_arr[16 * 16]);
  int16_t *out = out_arr;

  /* transform rows */
  for (i = 0; i < 2; ++i) {
    /* process 16 * 8 block */
    vp9_idct16_1d_rows_msa((input + (i << 7)), (out + (i << 7)));
  }

  /* transform columns */
  for (i = 0; i < 2; ++i) {
    /* process 8 * 16 block */
    vp9_idct16_1d_columns_addblk_msa((out + (i << 3)), (dest + (i << 3)),
                                     dest_stride);
  }
}

void vp9_idct16x16_10_add_msa(const int16_t *input, uint8_t *dest,
                              int32_t dest_stride) {
  uint8_t i;
  DECLARE_ALIGNED(32, int16_t, out_arr[16 * 16]);
  int16_t *out = out_arr;

  /* process 16 * 8 block */
  vp9_idct16_1d_rows_msa(input, out);

  /* short case just considers top 4 rows as valid output */
  out += 4 * 16;
  for (i = 12; i--;) {
    __asm__ __volatile__ (
        "sw     $zero,   0(%[out])     \n\t"
        "sw     $zero,   4(%[out])     \n\t"
        "sw     $zero,   8(%[out])     \n\t"
        "sw     $zero,  12(%[out])     \n\t"
        "sw     $zero,  16(%[out])     \n\t"
        "sw     $zero,  20(%[out])     \n\t"
        "sw     $zero,  24(%[out])     \n\t"
        "sw     $zero,  28(%[out])     \n\t"

        :
        : [out] "r" (out)
    );

    out += 16;
  }

  out = out_arr;

  /* transform columns */
  for (i = 0; i < 2; ++i) {
    /* process 8 * 16 block */
    vp9_idct16_1d_columns_addblk_msa((out + (i << 3)), (dest + (i << 3)),
                                     dest_stride);
  }
}

void vp9_idct16x16_1_add_msa(const int16_t *input, uint8_t *dest,
                             int32_t dest_stride) {
  uint8_t i;
  int32_t const1;
  int16_t out;
  v8i16 const2, res0, res1, res2, res3, res4, res5, res6, res7;
  v16u8 dest0, dest1, dest2, dest3;
  v16u8 tmp0, tmp1, tmp2, tmp3;
  v16i8 zero = { 0 };

  out = dct_const_round_shift(input[0] * cospi_16_64);
  out = dct_const_round_shift(out * cospi_16_64);
  const1 = ROUND_POWER_OF_TWO(out, 6);

  const2 = __msa_fill_h(const1);

  for (i = 0; i < 4; ++i) {
    LOAD_4VECS_UB(dest, dest_stride, dest0, dest1, dest2, dest3);

    res0 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest0);
    res1 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest1);
    res2 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest2);
    res3 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest3);
    res4 = (v8i16)__msa_ilvl_b(zero, (v16i8)dest0);
    res5 = (v8i16)__msa_ilvl_b(zero, (v16i8)dest1);
    res6 = (v8i16)__msa_ilvl_b(zero, (v16i8)dest2);
    res7 = (v8i16)__msa_ilvl_b(zero, (v16i8)dest3);

    res0 += const2;
    res1 += const2;
    res2 += const2;
    res3 += const2;
    res4 += const2;
    res5 += const2;
    res6 += const2;
    res7 += const2;

    res0 = CLIP_UNSIGNED_CHAR_H(res0);
    res1 = CLIP_UNSIGNED_CHAR_H(res1);
    res2 = CLIP_UNSIGNED_CHAR_H(res2);
    res3 = CLIP_UNSIGNED_CHAR_H(res3);
    res4 = CLIP_UNSIGNED_CHAR_H(res4);
    res5 = CLIP_UNSIGNED_CHAR_H(res5);
    res6 = CLIP_UNSIGNED_CHAR_H(res6);
    res7 = CLIP_UNSIGNED_CHAR_H(res7);

    tmp0 = (v16u8)__msa_pckev_b((v16i8)res4, (v16i8)res0);
    tmp1 = (v16u8)__msa_pckev_b((v16i8)res5, (v16i8)res1);
    tmp2 = (v16u8)__msa_pckev_b((v16i8)res6, (v16i8)res2);
    tmp3 = (v16u8)__msa_pckev_b((v16i8)res7, (v16i8)res3);

    STORE_4VECS_UB(dest, dest_stride, tmp0, tmp1, tmp2, tmp3);
    dest += (4 * dest_stride);
  }
}

static void vp9_iadst16_1d_rows_msa(const int16_t *input, int16_t *output) {
  v8i16 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;
  v8i16 l0, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15;

  /* load input data */
  LOAD_16VECS_SH(input, 8,
                 l0, l8, l1, l9, l2, l10, l3, l11,
                 l4, l12, l5, l13, l6, l14, l7, l15);

  TRANSPOSE8x8_H_SH(l0, l1, l2, l3, l4, l5, l6, l7,
                    l0, l1, l2, l3, l4, l5, l6, l7);

  TRANSPOSE8x8_H_SH(l8, l9, l10, l11, l12, l13, l14, l15,
                    l8, l9, l10, l11, l12, l13, l14, l15);

  /* ADST in horizontal */
  VP9_IADST8x16_1D(l0, l1, l2, l3, l4, l5, l6, l7,
                   l8, l9, l10, l11, l12, l13, l14, l15,
                   r0, r1, r2, r3, r4, r5, r6, r7,
                   r8, r9, r10, r11, r12, r13, r14, r15);

  l1 = -r8;
  l3 = -r4;
  l13 = -r13;
  l15 = -r1;

  TRANSPOSE8x8_H_SH(r0, l1, r12, l3, r6, r14, r10, r2,
                    l0, l1, l2, l3, l4, l5, l6, l7);

  STORE_8VECS_SH(output, 16, l0, l1, l2, l3, l4, l5, l6, l7);

  TRANSPOSE8x8_H_SH(r3, r11, r15, r7, r5, l13, r9, l15,
                    l8, l9, l10, l11, l12, l13, l14, l15);

  STORE_8VECS_SH((output + 8), 16, l8, l9, l10, l11, l12, l13, l14, l15);
}

static void vp9_iadst16_1d_columns_addblk_msa(int16_t *input, uint8_t *dest,
                                              int32_t dest_stride) {
  v8i16 v0, v2, v4, v6, k0, k1, k2, k3;
  v8i16 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;
  v8i16 out0, out1, out2, out3, out4, out5, out6, out7;
  v8i16 out8, out9, out10, out11, out12, out13, out14, out15;
  v8i16 g0, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15;
  v8i16 h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11;
  v8i16 res0, res1, res2, res3, res4, res5, res6, res7;
  v8i16 res8, res9, res10, res11, res12, res13, res14, res15;
  v16u8 dest0, dest1, dest2, dest3, dest4, dest5, dest6, dest7;
  v16u8 dest8, dest9, dest10, dest11, dest12, dest13, dest14, dest15;
  v16i8 zero = { 0 };

  r0 = LOAD_SH(input + 0 * 16);
  r3 = LOAD_SH(input + 3 * 16);
  r4 = LOAD_SH(input + 4 * 16);
  r7 = LOAD_SH(input + 7 * 16);
  r8 = LOAD_SH(input + 8 * 16);
  r11 = LOAD_SH(input + 11 * 16);
  r12 = LOAD_SH(input + 12 * 16);
  r15 = LOAD_SH(input + 15 * 16);

  /* stage 1 */
  k0 = SET_COSPI_PAIR(cospi_1_64, cospi_31_64);
  k1 = SET_COSPI_PAIR(cospi_31_64, -cospi_1_64);
  k2 = SET_COSPI_PAIR(cospi_17_64, cospi_15_64);
  k3 = SET_COSPI_PAIR(cospi_15_64, -cospi_17_64);
  VP9_MADD_BF(r15, r0, r7, r8, k0, k1, k2, k3, g0, g1, g2, g3);

  k0 = SET_COSPI_PAIR(cospi_9_64, cospi_23_64);
  k1 = SET_COSPI_PAIR(cospi_23_64, -cospi_9_64);
  k2 = SET_COSPI_PAIR(cospi_25_64, cospi_7_64);
  k3 = SET_COSPI_PAIR(cospi_7_64, -cospi_25_64);
  VP9_MADD_BF(r11, r4, r3, r12, k0, k1, k2, k3, g8, g9, g10, g11);

  BUTTERFLY_4(g0, g2, g10, g8, h8, h9, v2, v0);

  k0 = SET_COSPI_PAIR(cospi_4_64, cospi_28_64);
  k1 = SET_COSPI_PAIR(cospi_28_64, -cospi_4_64);
  k2 = SET_COSPI_PAIR(-cospi_28_64, cospi_4_64);
  VP9_MADD_BF(g1, g3, g9, g11, k0, k1, k2, k0, h0, h1, h2, h3);

  r1 = LOAD_SH(input + 1 * 16);
  r2 = LOAD_SH(input + 2 * 16);
  r5 = LOAD_SH(input + 5 * 16);
  r6 = LOAD_SH(input + 6 * 16);
  r9 = LOAD_SH(input + 9 * 16);
  r10 = LOAD_SH(input + 10 * 16);
  r13 = LOAD_SH(input + 13 * 16);
  r14 = LOAD_SH(input + 14 * 16);

  k0 = SET_COSPI_PAIR(cospi_5_64, cospi_27_64);
  k1 = SET_COSPI_PAIR(cospi_27_64, -cospi_5_64);
  k2 = SET_COSPI_PAIR(cospi_21_64, cospi_11_64);
  k3 = SET_COSPI_PAIR(cospi_11_64, -cospi_21_64);
  VP9_MADD_BF(r13, r2, r5, r10, k0, k1, k2, k3, g4, g5, g6, g7);

  k0 = SET_COSPI_PAIR(cospi_13_64, cospi_19_64);
  k1 = SET_COSPI_PAIR(cospi_19_64, -cospi_13_64);
  k2 = SET_COSPI_PAIR(cospi_29_64, cospi_3_64);
  k3 = SET_COSPI_PAIR(cospi_3_64, -cospi_29_64);
  VP9_MADD_BF(r9, r6, r1, r14, k0, k1, k2, k3, g12, g13, g14, g15);

  BUTTERFLY_4(g4, g6, g14, g12, h10, h11, v6, v4);

  BUTTERFLY_4(h8, h9, h11, h10, out0, out1, h11, h10);
  out1 = -out1;
  out0 = __msa_srari_h(out0, 6);
  out1 = __msa_srari_h(out1, 6);
  dest0 = LOAD_UB(dest + 0 * dest_stride);
  dest1 = LOAD_UB(dest + 15 * dest_stride);
  res0 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest0);
  res1 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest1);
  res0 += out0;
  res1 += out1;
  res0 = CLIP_UNSIGNED_CHAR_H(res0);
  res1 = CLIP_UNSIGNED_CHAR_H(res1);
  res0 = (v8i16)__msa_pckev_b((v16i8)res0, (v16i8)res0);
  res1 = (v8i16)__msa_pckev_b((v16i8)res1, (v16i8)res1);
  STORE_DWORD(dest, __msa_copy_u_d((v2i64)res0, 0));
  STORE_DWORD(dest + 15 * dest_stride, __msa_copy_u_d((v2i64)res1, 0));

  k0 = SET_COSPI_PAIR(cospi_12_64, cospi_20_64);
  k1 = SET_COSPI_PAIR(-cospi_20_64, cospi_12_64);
  k2 = SET_COSPI_PAIR(cospi_20_64, -cospi_12_64);
  VP9_MADD_BF(g7, g5, g15, g13, k0, k1, k2, k0, h4, h5, h6, h7);

  BUTTERFLY_4(h0, h2, h6, h4, out8, out9, out11, out10);
  out8 = -out8;

  out8 = __msa_srari_h(out8, 6);
  out9 = __msa_srari_h(out9, 6);
  dest8 = LOAD_UB(dest + 1 * dest_stride);
  dest9 = LOAD_UB(dest + 14 * dest_stride);
  res8 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest8);
  res9 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest9);
  res8 += out8;
  res9 += out9;
  res8 = CLIP_UNSIGNED_CHAR_H(res8);
  res9 = CLIP_UNSIGNED_CHAR_H(res9);
  res8 = (v8i16)__msa_pckev_b((v16i8)res8, (v16i8)res8);
  res9 = (v8i16)__msa_pckev_b((v16i8)res9, (v16i8)res9);
  STORE_DWORD(dest + dest_stride, __msa_copy_u_d((v2i64)res8, 0));
  STORE_DWORD(dest + 14 * dest_stride, __msa_copy_u_d((v2i64)res9, 0));

  k0 = SET_COSPI_PAIR(cospi_8_64, cospi_24_64);
  k1 = SET_COSPI_PAIR(cospi_24_64, -cospi_8_64);
  k2 = SET_COSPI_PAIR(-cospi_24_64, cospi_8_64);
  VP9_MADD_BF(v0, v2, v4, v6, k0, k1, k2, k0, out4, out6, out5, out7);
  out4 = -out4;
  out4 = __msa_srari_h(out4, 6);
  out5 = __msa_srari_h(out5, 6);
  dest4 = LOAD_UB(dest + 3 * dest_stride);
  dest5 = LOAD_UB(dest + 12 * dest_stride);
  res4 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest4);
  res5 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest5);
  res4 += out4;
  res5 += out5;
  res4 = CLIP_UNSIGNED_CHAR_H(res4);
  res5 = CLIP_UNSIGNED_CHAR_H(res5);
  res4 = (v8i16)__msa_pckev_b((v16i8)res4, (v16i8)res4);
  res5 = (v8i16)__msa_pckev_b((v16i8)res5, (v16i8)res5);
  STORE_DWORD(dest + 3 * dest_stride, __msa_copy_u_d((v2i64)res4, 0));
  STORE_DWORD(dest + 12 * dest_stride, __msa_copy_u_d((v2i64)res5, 0));

  VP9_MADD_BF(h1, h3, h5, h7, k0, k1, k2, k0, out12, out14, out13, out15);
  out13 = -out13;
  out12 = __msa_srari_h(out12, 6);
  out13 = __msa_srari_h(out13, 6);
  dest12 = LOAD_UB(dest + 2 * dest_stride);
  dest13 = LOAD_UB(dest + 13 * dest_stride);
  res12 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest12);
  res13 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest13);
  res12 += out12;
  res13 += out13;
  res12 = CLIP_UNSIGNED_CHAR_H(res12);
  res13 = CLIP_UNSIGNED_CHAR_H(res13);
  res12 = (v8i16)__msa_pckev_b((v16i8)res12, (v16i8)res12);
  res13 = (v8i16)__msa_pckev_b((v16i8)res13, (v16i8)res13);
  STORE_DWORD(dest + 2 * dest_stride, __msa_copy_u_d((v2i64)res12, 0));
  STORE_DWORD(dest + 13 * dest_stride, __msa_copy_u_d((v2i64)res13, 0));

  k0 = SET_COSPI_PAIR(cospi_16_64, cospi_16_64);
  k3 = SET_COSPI_PAIR(-cospi_16_64, cospi_16_64);
  VP9_MADD_SHORT(out6, out7, k0, k3, out6, out7);
  out6 = __msa_srari_h(out6, 6);
  out7 = __msa_srari_h(out7, 6);
  dest6 = LOAD_UB(dest + 4 * dest_stride);
  dest7 = LOAD_UB(dest + 11 * dest_stride);
  res6 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest6);
  res7 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest7);
  res6 += out6;
  res7 += out7;
  res6 = CLIP_UNSIGNED_CHAR_H(res6);
  res7 = CLIP_UNSIGNED_CHAR_H(res7);
  res6 = (v8i16)__msa_pckev_b((v16i8)res6, (v16i8)res6);
  res7 = (v8i16)__msa_pckev_b((v16i8)res7, (v16i8)res7);
  STORE_DWORD(dest + 4 * dest_stride, __msa_copy_u_d((v2i64)res6, 0));
  STORE_DWORD(dest + 11 * dest_stride, __msa_copy_u_d((v2i64)res7, 0));

  VP9_MADD_SHORT(out10, out11, k0, k3, out10, out11);
  out10 = __msa_srari_h(out10, 6);
  out11 = __msa_srari_h(out11, 6);
  dest10 = LOAD_UB(dest + 6 * dest_stride);
  dest11 = LOAD_UB(dest + 9 * dest_stride);
  res10 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest10);
  res11 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest11);
  res10 += out10;
  res11 += out11;
  res10 = CLIP_UNSIGNED_CHAR_H(res10);
  res11 = CLIP_UNSIGNED_CHAR_H(res11);
  res10 = (v8i16)__msa_pckev_b((v16i8)res10, (v16i8)res10);
  res11 = (v8i16)__msa_pckev_b((v16i8)res11, (v16i8)res11);
  STORE_DWORD(dest + 6 * dest_stride, __msa_copy_u_d((v2i64)res10, 0));
  STORE_DWORD(dest + 9 * dest_stride, __msa_copy_u_d((v2i64)res11, 0));

  k1 = SET_COSPI_PAIR(-cospi_16_64, -cospi_16_64);
  k2 = SET_COSPI_PAIR(cospi_16_64, -cospi_16_64);
  VP9_MADD_SHORT(h10, h11, k1, k2, out2, out3);
  out2 = __msa_srari_h(out2, 6);
  out3 = __msa_srari_h(out3, 6);
  dest2 = LOAD_UB(dest + 7 * dest_stride);
  dest3 = LOAD_UB(dest + 8 * dest_stride);
  res2 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest2);
  res3 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest3);
  res2 += out2;
  res3 += out3;
  res2 = CLIP_UNSIGNED_CHAR_H(res2);
  res3 = CLIP_UNSIGNED_CHAR_H(res3);
  res2 = (v8i16)__msa_pckev_b((v16i8)res2, (v16i8)res2);
  res3 = (v8i16)__msa_pckev_b((v16i8)res3, (v16i8)res3);
  STORE_DWORD(dest + 7 * dest_stride, __msa_copy_u_d((v2i64)res2, 0));
  STORE_DWORD(dest + 8 * dest_stride, __msa_copy_u_d((v2i64)res3, 0));

  VP9_MADD_SHORT(out14, out15, k1, k2, out14, out15);
  out14 = __msa_srari_h(out14, 6);
  out15 = __msa_srari_h(out15, 6);
  dest14 = LOAD_UB(dest + 5 * dest_stride);
  dest15 = LOAD_UB(dest + 10 * dest_stride);
  res14 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest14);
  res15 = (v8i16)__msa_ilvr_b(zero, (v16i8)dest15);
  res14 += out14;
  res15 += out15;
  res14 = CLIP_UNSIGNED_CHAR_H(res14);
  res15 = CLIP_UNSIGNED_CHAR_H(res15);
  res14 = (v8i16)__msa_pckev_b((v16i8)res14, (v16i8)res14);
  res15 = (v8i16)__msa_pckev_b((v16i8)res15, (v16i8)res15);
  STORE_DWORD(dest + 5 * dest_stride, __msa_copy_u_d((v2i64)res14, 0));
  STORE_DWORD(dest + 10 * dest_stride, __msa_copy_u_d((v2i64)res15, 0));
}

void vp9_iht16x16_256_add_msa(const int16_t *input, uint8_t *dest,
                              int32_t dest_stride, int32_t tx_type) {
  int32_t i;
  DECLARE_ALIGNED(32, int16_t, out[16 * 16]);
  int16_t *out_ptr = &out[0];

  switch (tx_type) {
    case DCT_DCT:
      /* transform rows */
      for (i = 0; i < 2; ++i) {
        /* process 16 * 8 block */
        vp9_idct16_1d_rows_msa((input + (i << 7)), (out_ptr + (i << 7)));
      }

      /* transform columns */
      for (i = 0; i < 2; ++i) {
        /* process 8 * 16 block */
        vp9_idct16_1d_columns_addblk_msa((out_ptr + (i << 3)),
                                         (dest + (i << 3)), dest_stride);
      }
      break;
    case ADST_DCT:
      /* transform rows */
      for (i = 0; i < 2; ++i) {
        /* process 16 * 8 block */
        vp9_idct16_1d_rows_msa((input + (i << 7)), (out_ptr + (i << 7)));
      }

      /* transform columns */
      for (i = 0; i < 2; ++i) {
        vp9_iadst16_1d_columns_addblk_msa((out_ptr + (i << 3)),
                                          (dest + (i << 3)), dest_stride);
      }
      break;
    case DCT_ADST:
      /* transform rows */
      for (i = 0; i < 2; ++i) {
        /* process 16 * 8 block */
        vp9_iadst16_1d_rows_msa((input + (i << 7)), (out_ptr + (i << 7)));
      }

      /* transform columns */
      for (i = 0; i < 2; ++i) {
        /* process 8 * 16 block */
        vp9_idct16_1d_columns_addblk_msa((out_ptr + (i << 3)),
                                         (dest + (i << 3)), dest_stride);
      }
      break;
    case ADST_ADST:
      /* transform rows */
      for (i = 0; i < 2; ++i) {
        /* process 16 * 8 block */
        vp9_iadst16_1d_rows_msa((input + (i << 7)), (out_ptr + (i << 7)));
      }

      /* transform columns */
      for (i = 0; i < 2; ++i) {
        vp9_iadst16_1d_columns_addblk_msa((out_ptr + (i << 3)),
                                          (dest + (i << 3)), dest_stride);
      }
      break;
    default:
      assert(0);
      break;
  }
}
