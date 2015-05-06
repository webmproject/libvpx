/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_ports/mem.h"
#include "vp9/common/vp9_idct.h"
#include "vp9/common/mips/msa/vp9_macros_msa.h"

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

#define VP9_ADDBLK_CLIP_AND_STORE_OFF_4H_VECS(dest, dest_stride,     \
                                              in0, in1, in2, in3) {  \
  uint64_t out0_m, out1_m, out2_m, out3_m;                           \
  v8i16 res0_m, res1_m, res2_m, res3_m;                              \
  v16u8 dest0_m, dest1_m, dest2_m, dest3_m;                          \
  v16i8 tmp0_m, tmp1_m;                                              \
  v16i8 zero_m = { 0 };                                              \
  uint8_t *dst_m = (uint8_t *)(dest);                                \
                                                                     \
  dest0_m = LOAD_UB(dst_m);                                          \
  dest1_m = LOAD_UB(dst_m + 4 * dest_stride);                        \
  dest2_m = LOAD_UB(dst_m + 8 * dest_stride);                        \
  dest3_m = LOAD_UB(dst_m + 12 * dest_stride);                       \
                                                                     \
  res0_m = (v8i16)__msa_ilvr_b(zero_m, (v16i8)dest0_m);              \
  res1_m = (v8i16)__msa_ilvr_b(zero_m, (v16i8)dest1_m);              \
  res2_m = (v8i16)__msa_ilvr_b(zero_m, (v16i8)dest2_m);              \
  res3_m = (v8i16)__msa_ilvr_b(zero_m, (v16i8)dest3_m);              \
                                                                     \
  res0_m += (v8i16)(in0);                                            \
  res1_m += (v8i16)(in1);                                            \
  res2_m += (v8i16)(in2);                                            \
  res3_m += (v8i16)(in3);                                            \
                                                                     \
  res0_m = CLIP_UNSIGNED_CHAR_H(res0_m);                             \
  res1_m = CLIP_UNSIGNED_CHAR_H(res1_m);                             \
  res2_m = CLIP_UNSIGNED_CHAR_H(res2_m);                             \
  res3_m = CLIP_UNSIGNED_CHAR_H(res3_m);                             \
                                                                     \
  tmp0_m = __msa_pckev_b((v16i8)res1_m, (v16i8)res0_m);              \
  tmp1_m = __msa_pckev_b((v16i8)res3_m, (v16i8)res2_m);              \
                                                                     \
  out0_m = __msa_copy_u_d((v2i64)tmp0_m, 0);                         \
  out1_m = __msa_copy_u_d((v2i64)tmp0_m, 1);                         \
  out2_m = __msa_copy_u_d((v2i64)tmp1_m, 0);                         \
  out3_m = __msa_copy_u_d((v2i64)tmp1_m, 1);                         \
                                                                     \
  STORE_DWORD(dst_m, out0_m);                                        \
  dst_m += (4 * dest_stride);                                        \
  STORE_DWORD(dst_m, out1_m);                                        \
  dst_m += (4 * dest_stride);                                        \
  STORE_DWORD(dst_m, out2_m);                                        \
  dst_m += (4 * dest_stride);                                        \
  STORE_DWORD(dst_m, out3_m);                                        \
}

static void vp9_idct32x8_row_transpose_store(const int16_t *input,
                                             int16_t *tmp_buf) {
  v8i16 m0, m1, m2, m3, m4, m5, m6, m7;
  v8i16 n0, n1, n2, n3, n4, n5, n6, n7;

  /* 1st & 2nd 8x8 */
  LOAD_8VECS_SH(input, 32, m0, n0, m1, n1, m2, n2, m3, n3);
  LOAD_8VECS_SH((input + 8), 32, m4, n4, m5, n5, m6, n6, m7, n7);
  TRANSPOSE8x8_H_SH(m0, n0, m1, n1, m2, n2, m3, n3,
                    m0, n0, m1, n1, m2, n2, m3, n3);
  TRANSPOSE8x8_H_SH(m4, n4, m5, n5, m6, n6, m7, n7,
                    m4, n4, m5, n5, m6, n6, m7, n7);
  STORE_4VECS_SH((tmp_buf), 8, m0, n0, m1, n1);
  STORE_4VECS_SH((tmp_buf + 4 * 8), 8, m2, n2, m3, n3);
  STORE_4VECS_SH((tmp_buf + 8 * 8), 8, m4, n4, m5, n5);
  STORE_4VECS_SH((tmp_buf + 12 * 8), 8, m6, n6, m7, n7);

  /* 3rd & 4th 8x8 */
  LOAD_8VECS_SH((input + 16), 32, m0, n0, m1, n1, m2, n2, m3, n3);
  LOAD_8VECS_SH((input + 24), 32, m4, n4, m5, n5, m6, n6, m7, n7);
  TRANSPOSE8x8_H_SH(m0, n0, m1, n1, m2, n2, m3, n3,
                    m0, n0, m1, n1, m2, n2, m3, n3);
  TRANSPOSE8x8_H_SH(m4, n4, m5, n5, m6, n6, m7, n7,
                    m4, n4, m5, n5, m6, n6, m7, n7);
  STORE_4VECS_SH((tmp_buf + 16 * 8), 8, m0, n0, m1, n1);
  STORE_4VECS_SH((tmp_buf + 20 * 8), 8, m2, n2, m3, n3);
  STORE_4VECS_SH((tmp_buf + 24 * 8), 8, m4, n4, m5, n5);
  STORE_4VECS_SH((tmp_buf + 28 * 8), 8, m6, n6, m7, n7);
}

static void vp9_idct32x8_row_even_process_store(int16_t *tmp_buf,
                                                int16_t *tmp_eve_buf) {
  v8i16 vec0, vec1, vec2, vec3, loc0, loc1, loc2, loc3;
  v8i16 reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7;
  v8i16 stp0, stp1, stp2, stp3, stp4, stp5, stp6, stp7;

  /* Even stage 1 */
  LOAD_8VECS_SH(tmp_buf, 32, reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7);

  DOTP_CONST_PAIR(reg1, reg7, cospi_28_64, cospi_4_64, reg1, reg7);
  DOTP_CONST_PAIR(reg5, reg3, cospi_12_64, cospi_20_64, reg5, reg3);

  vec0 = reg1 - reg5;
  vec1 = reg1 + reg5;
  vec2 = reg7 - reg3;
  vec3 = reg7 + reg3;

  DOTP_CONST_PAIR(vec2, vec0, cospi_16_64, cospi_16_64, loc2, loc3);

  loc1 = vec3;
  loc0 = vec1;

  DOTP_CONST_PAIR(reg0, reg4, cospi_16_64, cospi_16_64, reg0, reg4);
  DOTP_CONST_PAIR(reg2, reg6, cospi_24_64, cospi_8_64, reg2, reg6);

  vec0 = reg4 - reg6;
  vec1 = reg4 + reg6;
  vec2 = reg0 - reg2;
  vec3 = reg0 + reg2;

  stp4 = vec0 - loc0;
  stp3 = vec0 + loc0;
  stp7 = vec1 - loc1;
  stp0 = vec1 + loc1;
  stp5 = vec2 - loc2;
  stp2 = vec2 + loc2;
  stp6 = vec3 - loc3;
  stp1 = vec3 + loc3;

  /* Even stage 2 */
  LOAD_8VECS_SH((tmp_buf + 16), 32,
                reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7);

  DOTP_CONST_PAIR(reg0, reg7, cospi_30_64, cospi_2_64, reg0, reg7);
  DOTP_CONST_PAIR(reg4, reg3, cospi_14_64, cospi_18_64, reg4, reg3);
  DOTP_CONST_PAIR(reg2, reg5, cospi_22_64, cospi_10_64, reg2, reg5);
  DOTP_CONST_PAIR(reg6, reg1, cospi_6_64, cospi_26_64, reg6, reg1);

  vec0 = reg0 + reg4;
  reg0 = reg0 - reg4;
  reg4 = reg6 + reg2;
  reg6 = reg6 - reg2;
  reg2 = reg1 + reg5;
  reg1 = reg1 - reg5;
  reg5 = reg7 + reg3;
  reg7 = reg7 - reg3;
  reg3 = vec0;

  vec1 = reg2;
  reg2 = reg3 + reg4;
  reg3 = reg3 - reg4;
  reg4 = reg5 - vec1;
  reg5 = reg5 + vec1;

  DOTP_CONST_PAIR(reg7, reg0, cospi_24_64, cospi_8_64, reg0, reg7);
  DOTP_CONST_PAIR((-reg6), reg1, cospi_24_64, cospi_8_64, reg6, reg1);

  vec0 = reg0 - reg6;
  reg0 = reg0 + reg6;
  vec1 = reg7 - reg1;
  reg7 = reg7 + reg1;

  DOTP_CONST_PAIR(vec1, vec0, cospi_16_64, cospi_16_64, reg6, reg1);
  DOTP_CONST_PAIR(reg4, reg3, cospi_16_64, cospi_16_64, reg3, reg4);

  /* Even stage 3 : Dependency on Even stage 1 & Even stage 2 */
  loc0 = stp0 - reg5;
  loc1 = stp0 + reg5;
  loc2 = stp1 - reg7;
  loc3 = stp1 + reg7;
  STORE_SH(loc0, (tmp_eve_buf + 15 * 8));
  STORE_SH(loc1, (tmp_eve_buf));
  STORE_SH(loc2, (tmp_eve_buf + 14 * 8));
  STORE_SH(loc3, (tmp_eve_buf + 8));

  loc0 = stp2 - reg1;
  loc1 = stp2 + reg1;
  loc2 = stp3 - reg4;
  loc3 = stp3 + reg4;
  STORE_SH(loc0, (tmp_eve_buf + 13 * 8));
  STORE_SH(loc1, (tmp_eve_buf + 2 * 8));
  STORE_SH(loc2, (tmp_eve_buf + 12 * 8));
  STORE_SH(loc3, (tmp_eve_buf + 3 * 8));

  /* Store 8 */
  loc0 = stp4 - reg3;
  loc1 = stp4 + reg3;
  loc2 = stp5 - reg6;
  loc3 = stp5 + reg6;
  STORE_SH(loc0, (tmp_eve_buf + 11 * 8));
  STORE_SH(loc1, (tmp_eve_buf + 4 * 8));
  STORE_SH(loc2, (tmp_eve_buf + 10 * 8));
  STORE_SH(loc3, (tmp_eve_buf + 5 * 8));

  loc0 = stp6 - reg0;
  loc1 = stp6 + reg0;
  loc2 = stp7 - reg2;
  loc3 = stp7 + reg2;
  STORE_SH(loc0, (tmp_eve_buf + 9 * 8));
  STORE_SH(loc1, (tmp_eve_buf + 6 * 8));
  STORE_SH(loc2, (tmp_eve_buf + 8 * 8));
  STORE_SH(loc3, (tmp_eve_buf + 7 * 8));
}

static void vp9_idct32x8_row_odd_process_store(int16_t *tmp_buf,
                                               int16_t *tmp_odd_buf) {
  v8i16 vec0, vec1, vec2, vec3, loc0, loc1, loc2, loc3;
  v8i16 reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7;

  /* Odd stage 1 */
  reg0 = LOAD_SH(tmp_buf + 8);
  reg1 = LOAD_SH(tmp_buf + 7 * 8);
  reg2 = LOAD_SH(tmp_buf + 9 * 8);
  reg3 = LOAD_SH(tmp_buf + 15 * 8);
  reg4 = LOAD_SH(tmp_buf + 17 * 8);
  reg5 = LOAD_SH(tmp_buf + 23 * 8);
  reg6 = LOAD_SH(tmp_buf + 25 * 8);
  reg7 = LOAD_SH(tmp_buf + 31 * 8);

  DOTP_CONST_PAIR(reg0, reg7, cospi_31_64, cospi_1_64, reg0, reg7);
  DOTP_CONST_PAIR(reg4, reg3, cospi_15_64, cospi_17_64, reg3, reg4);
  DOTP_CONST_PAIR(reg2, reg5, cospi_23_64, cospi_9_64, reg2, reg5);
  DOTP_CONST_PAIR(reg6, reg1, cospi_7_64, cospi_25_64, reg1, reg6);

  vec0 = reg0 + reg3;
  reg0 = reg0 - reg3;
  reg3 = reg7 + reg4;
  reg7 = reg7 - reg4;
  reg4 = reg1 + reg2;
  reg1 = reg1 - reg2;
  reg2 = reg6 + reg5;
  reg6 = reg6 - reg5;
  reg5 = vec0;

  /* 4 Stores */
  vec0 = reg5 + reg4;
  vec1 = reg3 + reg2;
  STORE_SH(vec0, (tmp_odd_buf + 4 * 8));
  STORE_SH(vec1, (tmp_odd_buf + 5 * 8));

  vec0 = reg5 - reg4;
  vec1 = reg3 - reg2;
  DOTP_CONST_PAIR(vec1, vec0, cospi_24_64, cospi_8_64, vec0, vec1);
  STORE_SH(vec0, (tmp_odd_buf));
  STORE_SH(vec1, (tmp_odd_buf + 8));

  /* 4 Stores */
  DOTP_CONST_PAIR(reg7, reg0, cospi_28_64, cospi_4_64, reg0, reg7);
  DOTP_CONST_PAIR(reg6, reg1, -cospi_4_64, cospi_28_64, reg1, reg6);

  vec0 = reg0 + reg1;
  vec2 = reg7 - reg6;
  vec1 = reg7 + reg6;
  vec3 = reg0 - reg1;
  STORE_SH(vec0, (tmp_odd_buf + 6 * 8));
  STORE_SH(vec1, (tmp_odd_buf + 7 * 8));

  DOTP_CONST_PAIR(vec2, vec3, cospi_24_64, cospi_8_64, vec2, vec3);
  STORE_SH(vec2, (tmp_odd_buf + 2 * 8));
  STORE_SH(vec3, (tmp_odd_buf + 3 * 8));

  /* Odd stage 2 */

  /* 8 loads */
  reg0 = LOAD_SH(tmp_buf + 3 * 8);
  reg1 = LOAD_SH(tmp_buf + 5 * 8);
  reg2 = LOAD_SH(tmp_buf + 11 * 8);
  reg3 = LOAD_SH(tmp_buf + 13 * 8);
  reg4 = LOAD_SH(tmp_buf + 19 * 8);
  reg5 = LOAD_SH(tmp_buf + 21 * 8);
  reg6 = LOAD_SH(tmp_buf + 27 * 8);
  reg7 = LOAD_SH(tmp_buf + 29 * 8);

  DOTP_CONST_PAIR(reg1, reg6, cospi_27_64, cospi_5_64, reg1, reg6);
  DOTP_CONST_PAIR(reg5, reg2, cospi_11_64, cospi_21_64, reg2, reg5);
  DOTP_CONST_PAIR(reg3, reg4, cospi_19_64, cospi_13_64, reg3, reg4);
  DOTP_CONST_PAIR(reg7, reg0, cospi_3_64, cospi_29_64, reg0, reg7);

  /* 4 Stores */
  vec0 = reg1 - reg2;
  vec1 = reg6 - reg5;
  vec2 = reg0 - reg3;
  vec3 = reg7 - reg4;
  DOTP_CONST_PAIR(vec1, vec0, cospi_12_64, cospi_20_64, loc0, loc1);
  DOTP_CONST_PAIR(vec3, vec2, -cospi_20_64, cospi_12_64, loc2, loc3);

  vec2 = loc2 - loc0;
  vec3 = loc3 - loc1;
  vec0 = loc2 + loc0;
  vec1 = loc3 + loc1;
  STORE_SH(vec0, (tmp_odd_buf + 12 * 8));
  STORE_SH(vec1, (tmp_odd_buf + 15 * 8));

  DOTP_CONST_PAIR(vec3, vec2, -cospi_8_64, cospi_24_64, vec0, vec1);

  STORE_SH(vec0, (tmp_odd_buf + 10 * 8));
  STORE_SH(vec1, (tmp_odd_buf + 11 * 8));

  /* 4 Stores */
  vec0 = reg0 + reg3;
  vec1 = reg1 + reg2;
  vec2 = reg6 + reg5;
  vec3 = reg7 + reg4;
  reg0 = vec0 + vec1;
  reg1 = vec3 + vec2;
  reg2 = vec0 - vec1;
  reg3 = vec3 - vec2;
  STORE_SH(reg0, (tmp_odd_buf + 13 * 8));
  STORE_SH(reg1, (tmp_odd_buf + 14 * 8));

  DOTP_CONST_PAIR(reg3, reg2, -cospi_8_64, cospi_24_64, reg0, reg1);

  STORE_SH(reg0, (tmp_odd_buf + 8 * 8));
  STORE_SH(reg1, (tmp_odd_buf + 9 * 8));

  /* Odd stage 3 : Dependency on Odd stage 1 & Odd stage 2 */

  /* Load 8 & Store 8 */
  reg0 = LOAD_SH(tmp_odd_buf);
  reg1 = LOAD_SH(tmp_odd_buf + 1 * 8);
  reg2 = LOAD_SH(tmp_odd_buf + 2 * 8);
  reg3 = LOAD_SH(tmp_odd_buf + 3 * 8);
  reg4 = LOAD_SH(tmp_odd_buf + 8 * 8);
  reg5 = LOAD_SH(tmp_odd_buf + 9 * 8);
  reg6 = LOAD_SH(tmp_odd_buf + 10 * 8);
  reg7 = LOAD_SH(tmp_odd_buf + 11 * 8);

  loc0 = reg0 + reg4;
  loc1 = reg1 + reg5;
  loc2 = reg2 + reg6;
  loc3 = reg3 + reg7;
  STORE_SH(loc0, (tmp_odd_buf));
  STORE_SH(loc1, (tmp_odd_buf + 1 * 8));
  STORE_SH(loc2, (tmp_odd_buf + 2 * 8));
  STORE_SH(loc3, (tmp_odd_buf + 3 * 8));

  vec0 = reg0 - reg4;
  vec1 = reg1 - reg5;
  DOTP_CONST_PAIR(vec1, vec0, cospi_16_64, cospi_16_64, loc0, loc1);

  vec0 = reg2 - reg6;
  vec1 = reg3 - reg7;
  DOTP_CONST_PAIR(vec1, vec0, cospi_16_64, cospi_16_64, loc2, loc3);

  STORE_SH(loc0, (tmp_odd_buf + 8 * 8));
  STORE_SH(loc1, (tmp_odd_buf + 9 * 8));
  STORE_SH(loc2, (tmp_odd_buf + 10 * 8));
  STORE_SH(loc3, (tmp_odd_buf + 11 * 8));

  /* Load 8 & Store 8 */
  reg1 = LOAD_SH(tmp_odd_buf + 4 * 8);
  reg2 = LOAD_SH(tmp_odd_buf + 5 * 8);
  reg0 = LOAD_SH(tmp_odd_buf + 6 * 8);
  reg3 = LOAD_SH(tmp_odd_buf + 7 * 8);
  reg4 = LOAD_SH(tmp_odd_buf + 12 * 8);
  reg5 = LOAD_SH(tmp_odd_buf + 13 * 8);
  reg6 = LOAD_SH(tmp_odd_buf + 14 * 8);
  reg7 = LOAD_SH(tmp_odd_buf + 15 * 8);

  loc0 = reg0 + reg4;
  loc1 = reg1 + reg5;
  loc2 = reg2 + reg6;
  loc3 = reg3 + reg7;
  STORE_SH(loc0, (tmp_odd_buf + 4 * 8));
  STORE_SH(loc1, (tmp_odd_buf + 5 * 8));
  STORE_SH(loc2, (tmp_odd_buf + 6 * 8));
  STORE_SH(loc3, (tmp_odd_buf + 7 * 8));

  vec0 = reg0 - reg4;
  vec1 = reg3 - reg7;
  DOTP_CONST_PAIR(vec1, vec0, cospi_16_64, cospi_16_64, loc0, loc1);

  vec0 = reg1 - reg5;
  vec1 = reg2 - reg6;
  DOTP_CONST_PAIR(vec1, vec0, cospi_16_64, cospi_16_64, loc2, loc3);

  STORE_SH(loc0, (tmp_odd_buf + 12 * 8));
  STORE_SH(loc1, (tmp_odd_buf + 13 * 8));
  STORE_SH(loc2, (tmp_odd_buf + 14 * 8));
  STORE_SH(loc3, (tmp_odd_buf + 15 * 8));
}

static void vp9_idct_butterfly_transpose_store(int16_t *tmp_buf,
                                               int16_t *tmp_eve_buf,
                                               int16_t *tmp_odd_buf,
                                               int16_t *dest) {
  v8i16 vec0, vec1, vec2, vec3, loc0, loc1, loc2, loc3;
  v8i16 m0, m1, m2, m3, m4, m5, m6, m7;
  v8i16 n0, n1, n2, n3, n4, n5, n6, n7;

  /* FINAL BUTTERFLY : Dependency on Even & Odd */
  /* Total: 32 loads, 32 stores */
  vec0 = LOAD_SH(tmp_odd_buf);
  vec1 = LOAD_SH(tmp_odd_buf + 9 * 8);
  vec2 = LOAD_SH(tmp_odd_buf + 14 * 8);
  vec3 = LOAD_SH(tmp_odd_buf + 6 * 8);
  loc0 = LOAD_SH(tmp_eve_buf);
  loc1 = LOAD_SH(tmp_eve_buf + 8 * 8);
  loc2 = LOAD_SH(tmp_eve_buf + 4 * 8);
  loc3 = LOAD_SH(tmp_eve_buf + 12 * 8);

  m0 = (loc0 + vec3);
  STORE_SH((loc0 - vec3), (tmp_buf + 31 * 8));
  STORE_SH((loc1 - vec2), (tmp_buf + 23 * 8));
  m4 = (loc1 + vec2);
  STORE_SH((loc2 - vec1), (tmp_buf + 27 * 8));
  m2 = (loc2 + vec1);
  STORE_SH((loc3 - vec0), (tmp_buf + 19 * 8));
  m6 = (loc3 + vec0);

  /* Load 8 & Store 8 */
  vec0 = LOAD_SH(tmp_odd_buf + 4 * 8);
  vec1 = LOAD_SH(tmp_odd_buf + 13 * 8);
  vec2 = LOAD_SH(tmp_odd_buf + 10 * 8);
  vec3 = LOAD_SH(tmp_odd_buf + 3 * 8);
  loc0 = LOAD_SH(tmp_eve_buf + 2 * 8);
  loc1 = LOAD_SH(tmp_eve_buf + 10 * 8);
  loc2 = LOAD_SH(tmp_eve_buf + 6 * 8);
  loc3 = LOAD_SH(tmp_eve_buf + 14 * 8);

  m1 = (loc0 + vec3);
  STORE_SH((loc0 - vec3), (tmp_buf + 29 * 8));
  STORE_SH((loc1 - vec2), (tmp_buf + 21 * 8));
  m5 = (loc1 + vec2);
  STORE_SH((loc2 - vec1), (tmp_buf + 25 * 8));
  m3 = (loc2 + vec1);
  STORE_SH((loc3 - vec0), (tmp_buf + 17 * 8));
  m7 = (loc3 + vec0);

  /* Load 8 & Store 8 */
  vec0 = LOAD_SH(tmp_odd_buf + 2 * 8);
  vec1 = LOAD_SH(tmp_odd_buf + 11 * 8);
  vec2 = LOAD_SH(tmp_odd_buf + 12 * 8);
  vec3 = LOAD_SH(tmp_odd_buf + 7 * 8);
  loc0 = LOAD_SH(tmp_eve_buf + 1 * 8);
  loc1 = LOAD_SH(tmp_eve_buf + 9 * 8);
  loc2 = LOAD_SH(tmp_eve_buf + 5 * 8);
  loc3 = LOAD_SH(tmp_eve_buf + 13 * 8);

  n0 = (loc0 + vec3);
  STORE_SH((loc0 - vec3), (tmp_buf + 30 * 8));
  STORE_SH((loc1 - vec2), (tmp_buf + 22 * 8));
  n4 = (loc1 + vec2);
  STORE_SH((loc2 - vec1), (tmp_buf + 26 * 8));
  n2 = (loc2 + vec1);
  STORE_SH((loc3 - vec0), (tmp_buf + 18 * 8));
  n6 = (loc3 + vec0);

  /* Load 8 & Store 8 */
  vec0 = LOAD_SH(tmp_odd_buf + 5 * 8);
  vec1 = LOAD_SH(tmp_odd_buf + 15 * 8);
  vec2 = LOAD_SH(tmp_odd_buf + 8 * 8);
  vec3 = LOAD_SH(tmp_odd_buf + 1 * 8);
  loc0 = LOAD_SH(tmp_eve_buf + 3 * 8);
  loc1 = LOAD_SH(tmp_eve_buf + 11 * 8);
  loc2 = LOAD_SH(tmp_eve_buf + 7 * 8);
  loc3 = LOAD_SH(tmp_eve_buf + 15 * 8);

  n1 = (loc0 + vec3);
  STORE_SH((loc0 - vec3), (tmp_buf + 28 * 8));
  STORE_SH((loc1 - vec2), (tmp_buf + 20 * 8));
  n5 = (loc1 + vec2);
  STORE_SH((loc2 - vec1), (tmp_buf + 24 * 8));
  n3 = (loc2 + vec1);
  STORE_SH((loc3 - vec0), (tmp_buf + 16 * 8));
  n7 = (loc3 + vec0);

  /* Transpose : 16 vectors */
  /* 1st & 2nd 8x8 */
  TRANSPOSE8x8_H_SH(m0, n0, m1, n1, m2, n2, m3, n3,
                    m0, n0, m1, n1, m2, n2, m3, n3);
  STORE_4VECS_SH((dest + 0), 32, m0, n0, m1, n1);
  STORE_4VECS_SH((dest + 4 * 32), 32, m2, n2, m3, n3);

  TRANSPOSE8x8_H_SH(m4, n4, m5, n5, m6, n6, m7, n7,
                    m4, n4, m5, n5, m6, n6, m7, n7);
  STORE_4VECS_SH((dest + 8), 32, m4, n4, m5, n5);
  STORE_4VECS_SH((dest + 8 + 4 * 32), 32, m6, n6, m7, n7);

  /* 3rd & 4th 8x8 */
  LOAD_8VECS_SH((tmp_buf + 8 * 16), 8, m0, n0, m1, n1, m2, n2, m3, n3);
  LOAD_8VECS_SH((tmp_buf + 12 * 16), 8, m4, n4, m5, n5, m6, n6, m7, n7);
  TRANSPOSE8x8_H_SH(m0, n0, m1, n1, m2, n2, m3, n3,
                    m0, n0, m1, n1, m2, n2, m3, n3);
  STORE_4VECS_SH((dest + 16), 32, m0, n0, m1, n1);
  STORE_4VECS_SH((dest + 16 + 4 * 32), 32, m2, n2, m3, n3);

  TRANSPOSE8x8_H_SH(m4, n4, m5, n5, m6, n6, m7, n7,
                    m4, n4, m5, n5, m6, n6, m7, n7);
  STORE_4VECS_SH((dest + 24), 32, m4, n4, m5, n5);
  STORE_4VECS_SH((dest + 24 + 4 * 32), 32, m6, n6, m7, n7);
}

static void vp9_idct32x8_1d_rows_msa(const int16_t *input, int16_t *output) {
  DECLARE_ALIGNED(32, int16_t, tmp_buf[8 * 32]);
  DECLARE_ALIGNED(32, int16_t, tmp_odd_buf[16 * 8]);
  DECLARE_ALIGNED(32, int16_t, tmp_eve_buf[16 * 8]);

  vp9_idct32x8_row_transpose_store(input, &tmp_buf[0]);

  vp9_idct32x8_row_even_process_store(&tmp_buf[0], &tmp_eve_buf[0]);

  vp9_idct32x8_row_odd_process_store(&tmp_buf[0], &tmp_odd_buf[0]);

  vp9_idct_butterfly_transpose_store(&tmp_buf[0], &tmp_eve_buf[0],
                                     &tmp_odd_buf[0], output);
}

static void vp9_idct8x32_column_even_process_store(int16_t *tmp_buf,
                                                   int16_t *tmp_eve_buf) {
  v8i16 vec0, vec1, vec2, vec3, loc0, loc1, loc2, loc3;
  v8i16 reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7;
  v8i16 stp0, stp1, stp2, stp3, stp4, stp5, stp6, stp7;

  /* Even stage 1 */
  LOAD_8VECS_SH(tmp_buf, (4 * 32),
                reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7);

  DOTP_CONST_PAIR(reg1, reg7, cospi_28_64, cospi_4_64, reg1, reg7);
  DOTP_CONST_PAIR(reg5, reg3, cospi_12_64, cospi_20_64, reg5, reg3);

  vec0 = reg1 - reg5;
  vec1 = reg1 + reg5;
  vec2 = reg7 - reg3;
  vec3 = reg7 + reg3;

  DOTP_CONST_PAIR(vec2, vec0, cospi_16_64, cospi_16_64, loc2, loc3);

  loc1 = vec3;
  loc0 = vec1;

  DOTP_CONST_PAIR(reg0, reg4, cospi_16_64, cospi_16_64, reg0, reg4);
  DOTP_CONST_PAIR(reg2, reg6, cospi_24_64, cospi_8_64, reg2, reg6);

  vec0 = reg4 - reg6;
  vec1 = reg4 + reg6;
  vec2 = reg0 - reg2;
  vec3 = reg0 + reg2;

  stp4 = vec0 - loc0;
  stp3 = vec0 + loc0;
  stp7 = vec1 - loc1;
  stp0 = vec1 + loc1;
  stp5 = vec2 - loc2;
  stp2 = vec2 + loc2;
  stp6 = vec3 - loc3;
  stp1 = vec3 + loc3;

  /* Even stage 2 */
  /* Load 8 */
  LOAD_8VECS_SH((tmp_buf + 2 * 32), (4 * 32),
                reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7);

  DOTP_CONST_PAIR(reg0, reg7, cospi_30_64, cospi_2_64, reg0, reg7);
  DOTP_CONST_PAIR(reg4, reg3, cospi_14_64, cospi_18_64, reg4, reg3);
  DOTP_CONST_PAIR(reg2, reg5, cospi_22_64, cospi_10_64, reg2, reg5);
  DOTP_CONST_PAIR(reg6, reg1, cospi_6_64, cospi_26_64, reg6, reg1);

  vec0 = reg0 + reg4;
  reg0 = reg0 - reg4;
  reg4 = reg6 + reg2;
  reg6 = reg6 - reg2;
  reg2 = reg1 + reg5;
  reg1 = reg1 - reg5;
  reg5 = reg7 + reg3;
  reg7 = reg7 - reg3;
  reg3 = vec0;

  vec1 = reg2;
  reg2 = reg3 + reg4;
  reg3 = reg3 - reg4;
  reg4 = reg5 - vec1;
  reg5 = reg5 + vec1;

  DOTP_CONST_PAIR(reg7, reg0, cospi_24_64, cospi_8_64, reg0, reg7);
  DOTP_CONST_PAIR((-reg6), reg1, cospi_24_64, cospi_8_64, reg6, reg1);

  vec0 = reg0 - reg6;
  reg0 = reg0 + reg6;
  vec1 = reg7 - reg1;
  reg7 = reg7 + reg1;

  DOTP_CONST_PAIR(vec1, vec0, cospi_16_64, cospi_16_64, reg6, reg1);
  DOTP_CONST_PAIR(reg4, reg3, cospi_16_64, cospi_16_64, reg3, reg4);

  /* Even stage 3 : Dependency on Even stage 1 & Even stage 2 */
  /* Store 8 */
  loc0 = stp0 - reg5;
  loc1 = stp0 + reg5;
  loc2 = stp1 - reg7;
  loc3 = stp1 + reg7;
  STORE_SH(loc0, (tmp_eve_buf + 15 * 8));
  STORE_SH(loc1, (tmp_eve_buf));
  STORE_SH(loc2, (tmp_eve_buf + 14 * 8));
  STORE_SH(loc3, (tmp_eve_buf + 1 * 8));

  loc0 = stp2 - reg1;
  loc1 = stp2 + reg1;
  loc2 = stp3 - reg4;
  loc3 = stp3 + reg4;
  STORE_SH(loc0, (tmp_eve_buf + 13 * 8));
  STORE_SH(loc1, (tmp_eve_buf + 2 * 8));
  STORE_SH(loc2, (tmp_eve_buf + 12 * 8));
  STORE_SH(loc3, (tmp_eve_buf + 3 * 8));

  /* Store 8 */
  loc0 = stp4 - reg3;
  loc1 = stp4 + reg3;
  loc2 = stp5 - reg6;
  loc3 = stp5 + reg6;
  STORE_SH(loc0, (tmp_eve_buf + 11 * 8));
  STORE_SH(loc1, (tmp_eve_buf + 4 * 8));
  STORE_SH(loc2, (tmp_eve_buf + 10 * 8));
  STORE_SH(loc3, (tmp_eve_buf + 5 * 8));

  loc0 = stp6 - reg0;
  loc1 = stp6 + reg0;
  loc2 = stp7 - reg2;
  loc3 = stp7 + reg2;
  STORE_SH(loc0, (tmp_eve_buf + 9 * 8));
  STORE_SH(loc1, (tmp_eve_buf + 6 * 8));
  STORE_SH(loc2, (tmp_eve_buf + 8 * 8));
  STORE_SH(loc3, (tmp_eve_buf + 7 * 8));
}

static void vp9_idct8x32_column_odd_process_store(int16_t *tmp_buf,
                                                  int16_t *tmp_odd_buf) {
  v8i16 vec0, vec1, vec2, vec3, loc0, loc1, loc2, loc3;
  v8i16 reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7;

  /* Odd stage 1 */
  reg0 = LOAD_SH(tmp_buf + 32);
  reg1 = LOAD_SH(tmp_buf + 7 * 32);
  reg2 = LOAD_SH(tmp_buf + 9 * 32);
  reg3 = LOAD_SH(tmp_buf + 15 * 32);
  reg4 = LOAD_SH(tmp_buf + 17 * 32);
  reg5 = LOAD_SH(tmp_buf + 23 * 32);
  reg6 = LOAD_SH(tmp_buf + 25 * 32);
  reg7 = LOAD_SH(tmp_buf + 31 * 32);

  DOTP_CONST_PAIR(reg0, reg7, cospi_31_64, cospi_1_64, reg0, reg7);
  DOTP_CONST_PAIR(reg4, reg3, cospi_15_64, cospi_17_64, reg3, reg4);
  DOTP_CONST_PAIR(reg2, reg5, cospi_23_64, cospi_9_64, reg2, reg5);
  DOTP_CONST_PAIR(reg6, reg1, cospi_7_64, cospi_25_64, reg1, reg6);

  vec0 = reg0 + reg3;
  reg0 = reg0 - reg3;
  reg3 = reg7 + reg4;
  reg7 = reg7 - reg4;
  reg4 = reg1 + reg2;
  reg1 = reg1 - reg2;
  reg2 = reg6 + reg5;
  reg6 = reg6 - reg5;
  reg5 = vec0;

  /* 4 Stores */
  vec0 = reg5 + reg4;
  vec1 = reg3 + reg2;
  STORE_SH(vec0, (tmp_odd_buf + 4 * 8));
  STORE_SH(vec1, (tmp_odd_buf + 5 * 8));

  vec0 = reg5 - reg4;
  vec1 = reg3 - reg2;
  DOTP_CONST_PAIR(vec1, vec0, cospi_24_64, cospi_8_64, vec0, vec1);
  STORE_SH(vec0, (tmp_odd_buf));
  STORE_SH(vec1, (tmp_odd_buf + 1 * 8));

  /* 4 Stores */
  DOTP_CONST_PAIR(reg7, reg0, cospi_28_64, cospi_4_64, reg0, reg7);
  DOTP_CONST_PAIR(reg6, reg1, -cospi_4_64, cospi_28_64, reg1, reg6);

  vec0 = reg0 + reg1;
  vec2 = reg7 - reg6;
  vec1 = reg7 + reg6;
  vec3 = reg0 - reg1;
  STORE_SH(vec0, (tmp_odd_buf + 6 * 8));
  STORE_SH(vec1, (tmp_odd_buf + 7 * 8));

  DOTP_CONST_PAIR(vec2, vec3, cospi_24_64, cospi_8_64, vec2, vec3);
  STORE_SH(vec2, (tmp_odd_buf + 2 * 8));
  STORE_SH(vec3, (tmp_odd_buf + 3 * 8));

  /* Odd stage 2 */
  /* 8 loads */
  reg0 = LOAD_SH(tmp_buf + 3 * 32);
  reg1 = LOAD_SH(tmp_buf + 5 * 32);
  reg2 = LOAD_SH(tmp_buf + 11 * 32);
  reg3 = LOAD_SH(tmp_buf + 13 * 32);
  reg4 = LOAD_SH(tmp_buf + 19 * 32);
  reg5 = LOAD_SH(tmp_buf + 21 * 32);
  reg6 = LOAD_SH(tmp_buf + 27 * 32);
  reg7 = LOAD_SH(tmp_buf + 29 * 32);

  DOTP_CONST_PAIR(reg1, reg6, cospi_27_64, cospi_5_64, reg1, reg6);
  DOTP_CONST_PAIR(reg5, reg2, cospi_11_64, cospi_21_64, reg2, reg5);
  DOTP_CONST_PAIR(reg3, reg4, cospi_19_64, cospi_13_64, reg3, reg4);
  DOTP_CONST_PAIR(reg7, reg0, cospi_3_64, cospi_29_64, reg0, reg7);

  /* 4 Stores */
  vec0 = reg1 - reg2;
  vec1 = reg6 - reg5;
  vec2 = reg0 - reg3;
  vec3 = reg7 - reg4;
  DOTP_CONST_PAIR(vec1, vec0, cospi_12_64, cospi_20_64, loc0, loc1);
  DOTP_CONST_PAIR(vec3, vec2, -cospi_20_64, cospi_12_64, loc2, loc3);

  vec2 = loc2 - loc0;
  vec3 = loc3 - loc1;
  vec0 = loc2 + loc0;
  vec1 = loc3 + loc1;
  STORE_SH(vec0, (tmp_odd_buf + 12 * 8));
  STORE_SH(vec1, (tmp_odd_buf + 15 * 8));

  DOTP_CONST_PAIR(vec3, vec2, -cospi_8_64, cospi_24_64, vec0, vec1);

  STORE_SH(vec0, (tmp_odd_buf + 10 * 8));
  STORE_SH(vec1, (tmp_odd_buf + 11 * 8));

  /* 4 Stores */
  vec0 = reg0 + reg3;
  vec1 = reg1 + reg2;
  vec2 = reg6 + reg5;
  vec3 = reg7 + reg4;
  reg0 = vec0 + vec1;
  reg1 = vec3 + vec2;
  reg2 = vec0 - vec1;
  reg3 = vec3 - vec2;
  STORE_SH(reg0, (tmp_odd_buf + 13 * 8));
  STORE_SH(reg1, (tmp_odd_buf + 14 * 8));

  DOTP_CONST_PAIR(reg3, reg2, -cospi_8_64, cospi_24_64, reg0, reg1);

  STORE_SH(reg0, (tmp_odd_buf + 8 * 8));
  STORE_SH(reg1, (tmp_odd_buf + 9 * 8));

  /* Odd stage 3 : Dependency on Odd stage 1 & Odd stage 2 */
  /* Load 8 & Store 8 */
  reg0 = LOAD_SH(tmp_odd_buf);
  reg1 = LOAD_SH(tmp_odd_buf + 1 * 8);
  reg2 = LOAD_SH(tmp_odd_buf + 2 * 8);
  reg3 = LOAD_SH(tmp_odd_buf + 3 * 8);
  reg4 = LOAD_SH(tmp_odd_buf + 8 * 8);
  reg5 = LOAD_SH(tmp_odd_buf + 9 * 8);
  reg6 = LOAD_SH(tmp_odd_buf + 10 * 8);
  reg7 = LOAD_SH(tmp_odd_buf + 11 * 8);

  loc0 = reg0 + reg4;
  loc1 = reg1 + reg5;
  loc2 = reg2 + reg6;
  loc3 = reg3 + reg7;
  STORE_SH(loc0, (tmp_odd_buf));
  STORE_SH(loc1, (tmp_odd_buf + 1 * 8));
  STORE_SH(loc2, (tmp_odd_buf + 2 * 8));
  STORE_SH(loc3, (tmp_odd_buf + 3 * 8));

  vec0 = reg0 - reg4;
  vec1 = reg1 - reg5;
  DOTP_CONST_PAIR(vec1, vec0, cospi_16_64, cospi_16_64, loc0, loc1);

  vec0 = reg2 - reg6;
  vec1 = reg3 - reg7;
  DOTP_CONST_PAIR(vec1, vec0, cospi_16_64, cospi_16_64, loc2, loc3);

  STORE_SH(loc0, (tmp_odd_buf + 8 * 8));
  STORE_SH(loc1, (tmp_odd_buf + 9 * 8));
  STORE_SH(loc2, (tmp_odd_buf + 10 * 8));
  STORE_SH(loc3, (tmp_odd_buf + 11 * 8));

  /* Load 8 & Store 8 */
  reg1 = LOAD_SH(tmp_odd_buf + 4 * 8);
  reg2 = LOAD_SH(tmp_odd_buf + 5 * 8);
  reg0 = LOAD_SH(tmp_odd_buf + 6 * 8);
  reg3 = LOAD_SH(tmp_odd_buf + 7 * 8);
  reg4 = LOAD_SH(tmp_odd_buf + 12 * 8);
  reg5 = LOAD_SH(tmp_odd_buf + 13 * 8);
  reg6 = LOAD_SH(tmp_odd_buf + 14 * 8);
  reg7 = LOAD_SH(tmp_odd_buf + 15 * 8);

  loc0 = reg0 + reg4;
  loc1 = reg1 + reg5;
  loc2 = reg2 + reg6;
  loc3 = reg3 + reg7;
  STORE_SH(loc0, (tmp_odd_buf + 4 * 8));
  STORE_SH(loc1, (tmp_odd_buf + 5 * 8));
  STORE_SH(loc2, (tmp_odd_buf + 6 * 8));
  STORE_SH(loc3, (tmp_odd_buf + 7 * 8));

  vec0 = reg0 - reg4;
  vec1 = reg3 - reg7;
  DOTP_CONST_PAIR(vec1, vec0, cospi_16_64, cospi_16_64, loc0, loc1);

  vec0 = reg1 - reg5;
  vec1 = reg2 - reg6;
  DOTP_CONST_PAIR(vec1, vec0, cospi_16_64, cospi_16_64, loc2, loc3);

  STORE_SH(loc0, (tmp_odd_buf + 12 * 8));
  STORE_SH(loc1, (tmp_odd_buf + 13 * 8));
  STORE_SH(loc2, (tmp_odd_buf + 14 * 8));
  STORE_SH(loc3, (tmp_odd_buf + 15 * 8));
}

static void vp9_idct8x32_column_butterfly_addblk(int16_t *tmp_eve_buf,
                                                 int16_t *tmp_odd_buf,
                                                 uint8_t *dest,
                                                 int32_t dest_stride) {
  v8i16 vec0, vec1, vec2, vec3, loc0, loc1, loc2, loc3;
  v8i16 m0, m1, m2, m3, m4, m5, m6, m7;
  v8i16 n0, n1, n2, n3, n4, n5, n6, n7;

  /* FINAL BUTTERFLY : Dependency on Even & Odd */
  vec0 = LOAD_SH(tmp_odd_buf);
  vec1 = LOAD_SH(tmp_odd_buf + 9 * 8);
  vec2 = LOAD_SH(tmp_odd_buf + 14 * 8);
  vec3 = LOAD_SH(tmp_odd_buf + 6 * 8);
  loc0 = LOAD_SH(tmp_eve_buf);
  loc1 = LOAD_SH(tmp_eve_buf + 8 * 8);
  loc2 = LOAD_SH(tmp_eve_buf + 4 * 8);
  loc3 = LOAD_SH(tmp_eve_buf + 12 * 8);

  m0 = (loc0 + vec3);
  m4 = (loc1 + vec2);
  m2 = (loc2 + vec1);
  m6 = (loc3 + vec0);
  SRARI_H_4VECS_SH(m0, m2, m4, m6, m0, m2, m4, m6, 6);
  VP9_ADDBLK_CLIP_AND_STORE_OFF_4H_VECS(dest, dest_stride, m0, m2, m4, m6);

  m6 = (loc0 - vec3);
  m2 = (loc1 - vec2);
  m4 = (loc2 - vec1);
  m0 = (loc3 - vec0);
  SRARI_H_4VECS_SH(m0, m2, m4, m6, m0, m2, m4, m6, 6);
  VP9_ADDBLK_CLIP_AND_STORE_OFF_4H_VECS((dest + 19 * dest_stride),
                                        dest_stride, m0, m2, m4, m6);

  /* Load 8 & Store 8 */
  vec0 = LOAD_SH(tmp_odd_buf + 4 * 8);
  vec1 = LOAD_SH(tmp_odd_buf + 13 * 8);
  vec2 = LOAD_SH(tmp_odd_buf + 10 * 8);
  vec3 = LOAD_SH(tmp_odd_buf + 3 * 8);
  loc0 = LOAD_SH(tmp_eve_buf + 2 * 8);
  loc1 = LOAD_SH(tmp_eve_buf + 10 * 8);
  loc2 = LOAD_SH(tmp_eve_buf + 6 * 8);
  loc3 = LOAD_SH(tmp_eve_buf + 14 * 8);

  m1 = (loc0 + vec3);
  m5 = (loc1 + vec2);
  m3 = (loc2 + vec1);
  m7 = (loc3 + vec0);
  SRARI_H_4VECS_SH(m1, m3, m5, m7, m1, m3, m5, m7, 6);
  VP9_ADDBLK_CLIP_AND_STORE_OFF_4H_VECS((dest + 2 * dest_stride),
                                        dest_stride, m1, m3, m5, m7);

  m7 = (loc0 - vec3);
  m3 = (loc1 - vec2);
  m5 = (loc2 - vec1);
  m1 = (loc3 - vec0);
  SRARI_H_4VECS_SH(m1, m3, m5, m7, m1, m3, m5, m7, 6);
  VP9_ADDBLK_CLIP_AND_STORE_OFF_4H_VECS((dest + 17 * dest_stride),
                                        dest_stride, m1, m3, m5, m7);

  /* Load 8 & Store 8 */
  vec0 = LOAD_SH(tmp_odd_buf + 2 * 8);
  vec1 = LOAD_SH(tmp_odd_buf + 11 * 8);
  vec2 = LOAD_SH(tmp_odd_buf + 12 * 8);
  vec3 = LOAD_SH(tmp_odd_buf + 7 * 8);
  loc0 = LOAD_SH(tmp_eve_buf + 1 * 8);
  loc1 = LOAD_SH(tmp_eve_buf + 9 * 8);
  loc2 = LOAD_SH(tmp_eve_buf + 5 * 8);
  loc3 = LOAD_SH(tmp_eve_buf + 13 * 8);

  n0 = (loc0 + vec3);
  n4 = (loc1 + vec2);
  n2 = (loc2 + vec1);
  n6 = (loc3 + vec0);
  SRARI_H_4VECS_SH(n0, n2, n4, n6, n0, n2, n4, n6, 6);
  VP9_ADDBLK_CLIP_AND_STORE_OFF_4H_VECS((dest + 1 * dest_stride),
                                        dest_stride, n0, n2, n4, n6);

  n6 = (loc0 - vec3);
  n2 = (loc1 - vec2);
  n4 = (loc2 - vec1);
  n0 = (loc3 - vec0);
  SRARI_H_4VECS_SH(n0, n2, n4, n6, n0, n2, n4, n6, 6);
  VP9_ADDBLK_CLIP_AND_STORE_OFF_4H_VECS((dest + 18 * dest_stride),
                                        dest_stride, n0, n2, n4, n6);

  /* Load 8 & Store 8 */
  vec0 = LOAD_SH(tmp_odd_buf + 5 * 8);
  vec1 = LOAD_SH(tmp_odd_buf + 15 * 8);
  vec2 = LOAD_SH(tmp_odd_buf + 8 * 8);
  vec3 = LOAD_SH(tmp_odd_buf + 1 * 8);
  loc0 = LOAD_SH(tmp_eve_buf + 3 * 8);
  loc1 = LOAD_SH(tmp_eve_buf + 11 * 8);
  loc2 = LOAD_SH(tmp_eve_buf + 7 * 8);
  loc3 = LOAD_SH(tmp_eve_buf + 15 * 8);

  n1 = (loc0 + vec3);
  n5 = (loc1 + vec2);
  n3 = (loc2 + vec1);
  n7 = (loc3 + vec0);
  SRARI_H_4VECS_SH(n1, n3, n5, n7, n1, n3, n5, n7, 6);
  VP9_ADDBLK_CLIP_AND_STORE_OFF_4H_VECS((dest + 3 * dest_stride),
                                        dest_stride, n1, n3, n5, n7);

  n7 = (loc0 - vec3);
  n3 = (loc1 - vec2);
  n5 = (loc2 - vec1);
  n1 = (loc3 - vec0);
  SRARI_H_4VECS_SH(n1, n3, n5, n7, n1, n3, n5, n7, 6);
  VP9_ADDBLK_CLIP_AND_STORE_OFF_4H_VECS((dest + 16 * dest_stride),
                                        dest_stride, n1, n3, n5, n7);
}

static void vp9_idct8x32_1d_columns_addblk_msa(int16_t *input, uint8_t *dest,
                                               int32_t dest_stride) {
  DECLARE_ALIGNED(32, int16_t, tmp_odd_buf[16 * 8]);
  DECLARE_ALIGNED(32, int16_t, tmp_eve_buf[16 * 8]);

  vp9_idct8x32_column_even_process_store(input, &tmp_eve_buf[0]);

  vp9_idct8x32_column_odd_process_store(input, &tmp_odd_buf[0]);

  vp9_idct8x32_column_butterfly_addblk(&tmp_eve_buf[0], &tmp_odd_buf[0],
                                       dest, dest_stride);
}

void vp9_idct32x32_1024_add_msa(const int16_t *input, uint8_t *dest,
                                int32_t dest_stride) {
  int32_t i;
  DECLARE_ALIGNED(32, int16_t, out_arr[32 * 32]);
  int16_t *out_ptr = out_arr;

  /* transform rows */
  for (i = 0; i < 4; ++i) {
    /* process 32 * 8 block */
    vp9_idct32x8_1d_rows_msa((input + (i << 8)), (out_ptr + (i << 8)));
  }

  /* transform columns */
  for (i = 0; i < 4; ++i) {
    /* process 8 * 32 block */
    vp9_idct8x32_1d_columns_addblk_msa((out_ptr + (i << 3)), (dest + (i << 3)),
                                       dest_stride);
  }
}

void vp9_idct32x32_34_add_msa(const int16_t *input, uint8_t *dest,
                              int32_t dest_stride) {
  int32_t i;
  DECLARE_ALIGNED(32, int16_t, out_arr[32 * 32]);
  int16_t *out_ptr = out_arr;

  for (i = 32; i--;) {
    __asm__ __volatile__ (
        "sw     $zero,      0(%[out_ptr])     \n\t"
        "sw     $zero,      4(%[out_ptr])     \n\t"
        "sw     $zero,      8(%[out_ptr])     \n\t"
        "sw     $zero,     12(%[out_ptr])     \n\t"
        "sw     $zero,     16(%[out_ptr])     \n\t"
        "sw     $zero,     20(%[out_ptr])     \n\t"
        "sw     $zero,     24(%[out_ptr])     \n\t"
        "sw     $zero,     28(%[out_ptr])     \n\t"
        "sw     $zero,     32(%[out_ptr])     \n\t"
        "sw     $zero,     36(%[out_ptr])     \n\t"
        "sw     $zero,     40(%[out_ptr])     \n\t"
        "sw     $zero,     44(%[out_ptr])     \n\t"
        "sw     $zero,     48(%[out_ptr])     \n\t"
        "sw     $zero,     52(%[out_ptr])     \n\t"
        "sw     $zero,     56(%[out_ptr])     \n\t"
        "sw     $zero,     60(%[out_ptr])     \n\t"

        :
        : [out_ptr] "r" (out_ptr)
    );

    out_ptr += 32;
  }

  out_ptr = out_arr;

  /* rows: only upper-left 8x8 has non-zero coeff */
  vp9_idct32x8_1d_rows_msa(input, out_ptr);

  /* transform columns */
  for (i = 0; i < 4; ++i) {
    /* process 8 * 32 block */
    vp9_idct8x32_1d_columns_addblk_msa((out_ptr + (i << 3)), (dest + (i << 3)),
                                       dest_stride);
  }
}

void vp9_idct32x32_1_add_msa(const int16_t *input, uint8_t *dest,
                             int32_t dest_stride) {
  int32_t i, const1;
  v8i16 const2;
  int16_t out;
  v8i16 res0, res1, res2, res3, res4, res5, res6, res7;
  v16u8 dest0, dest1, dest2, dest3;
  v16u8 tmp0, tmp1, tmp2, tmp3;
  v16i8 zero = { 0 };

  out = dct_const_round_shift(input[0] * cospi_16_64);
  out = dct_const_round_shift(out * cospi_16_64);
  const1 = ROUND_POWER_OF_TWO(out, 6);

  const2 = __msa_fill_h(const1);

  for (i = 0; i < 16; ++i) {
    dest0 = LOAD_UB(dest);
    dest1 = LOAD_UB(dest + 16);
    dest2 = LOAD_UB(dest + dest_stride);
    dest3 = LOAD_UB(dest + dest_stride + 16);

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

    STORE_UB(tmp0, dest);
    STORE_UB(tmp1, dest + 16);
    dest += dest_stride;
    STORE_UB(tmp2, dest);
    STORE_UB(tmp3, dest + 16);
    dest += dest_stride;
  }
}
