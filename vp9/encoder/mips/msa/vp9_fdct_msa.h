/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_MIPS_MSA_VP9_FDCT_MSA_H_
#define VP9_ENCODER_MIPS_MSA_VP9_FDCT_MSA_H_

#include "vpx_ports/mem.h"
#include "vp9/common/vp9_idct.h"
#include "vp9/common/mips/msa/vp9_macros_msa.h"

#define VP9_DOTP_CONST_PAIR(reg0, reg1, cnst0, cnst1, out0, out1) {  \
  v8i16 k0_m = __msa_fill_h(cnst0);                                  \
  v4i32 s0_m, s1_m, s2_m, s3_m;                                      \
                                                                     \
  s0_m = (v4i32)__msa_fill_h(cnst1);                                 \
  k0_m = __msa_ilvev_h((v8i16)s0_m, k0_m);                           \
                                                                     \
  ILVRL_H2_SW((-reg1), reg0, s1_m, s0_m);                            \
  ILVRL_H2_SW(reg0, reg1, s3_m, s2_m);                               \
  DOTP_SH2_SW(s1_m, s0_m, k0_m, k0_m, s1_m, s0_m);                   \
  SRARI_W2_SW(s1_m, s0_m, DCT_CONST_BITS);                           \
  out0 = __msa_pckev_h((v8i16)s0_m, (v8i16)s1_m);                    \
                                                                     \
  DOTP_SH2_SW(s3_m, s2_m, k0_m, k0_m, s1_m, s0_m);                   \
  SRARI_W2_SW(s1_m, s0_m, DCT_CONST_BITS);                           \
  out1 = __msa_pckev_h((v8i16)s0_m, (v8i16)s1_m);                    \
}

#define VP9_DOT_ADD_SUB_SRARI_PCK(in0, in1, in2, in3, in4, in5, in6, in7,  \
                                  dst0, dst1, dst2, dst3) {                \
  v4i32 tp0_m, tp1_m, tp2_m, tp3_m, tp4_m;                                 \
  v4i32 tp5_m, tp6_m, tp7_m, tp8_m, tp9_m;                                 \
                                                                           \
  DOTP_SH4_SW(in0, in1, in0, in1, in4, in4, in5, in5,                      \
              tp0_m, tp2_m, tp3_m, tp4_m);                                 \
  DOTP_SH4_SW(in2, in3, in2, in3, in6, in6, in7, in7,                      \
              tp5_m, tp6_m, tp7_m, tp8_m);                                 \
  BUTTERFLY_4(tp0_m, tp3_m, tp7_m, tp5_m, tp1_m, tp9_m, tp7_m, tp5_m);     \
  BUTTERFLY_4(tp2_m, tp4_m, tp8_m, tp6_m, tp3_m, tp0_m, tp4_m, tp2_m);     \
  SRARI_W4_SW(tp1_m, tp9_m, tp7_m, tp5_m, DCT_CONST_BITS);                 \
  SRARI_W4_SW(tp3_m, tp0_m, tp4_m, tp2_m, DCT_CONST_BITS);                 \
  PCKEV_H4_SH(tp1_m, tp3_m, tp9_m, tp0_m, tp7_m, tp4_m, tp5_m, tp2_m,      \
              dst0, dst1, dst2, dst3);                                     \
}

#define VP9_DOT_SHIFT_RIGHT_PCK_H(in0, in1, in2) ({   \
  v8i16 dst_m;                                        \
  v4i32 tp0_m, tp1_m;                                 \
                                                      \
  DOTP_SH2_SW(in0, in1, in2, in2, tp1_m, tp0_m);      \
  SRARI_W2_SW(tp1_m, tp0_m, DCT_CONST_BITS);          \
  dst_m = __msa_pckev_h((v8i16)tp1_m, (v8i16)tp0_m);  \
                                                      \
  dst_m;                                              \
})

#define VP9_ADST8(in0, in1, in2, in3, in4, in5, in6, in7,                   \
                  out0, out1, out2, out3, out4, out5, out6, out7) {         \
  v8i16 cnst0_m, cnst1_m, cnst2_m, cnst3_m, cnst4_m;                        \
  v8i16 vec0_m, vec1_m, vec2_m, vec3_m, s0_m, s1_m;                         \
  v8i16 coeff0_m = { cospi_2_64, cospi_6_64, cospi_10_64, cospi_14_64,      \
                     cospi_18_64, cospi_22_64, cospi_26_64, cospi_30_64 };  \
  v8i16 coeff1_m = { cospi_8_64, -cospi_8_64, cospi_16_64, -cospi_16_64,    \
                     cospi_24_64, -cospi_24_64, 0, 0 };                     \
                                                                            \
  SPLATI_H2_SH(coeff0_m, 0, 7, cnst0_m, cnst1_m);                           \
  cnst2_m = -cnst0_m;                                                       \
  ILVEV_H2_SH(cnst0_m, cnst1_m, cnst1_m, cnst2_m, cnst0_m, cnst1_m);        \
  SPLATI_H2_SH(coeff0_m, 4, 3, cnst2_m, cnst3_m);                           \
  cnst4_m = -cnst2_m;                                                       \
  ILVEV_H2_SH(cnst2_m, cnst3_m, cnst3_m, cnst4_m, cnst2_m, cnst3_m);        \
                                                                            \
  ILVRL_H2_SH(in0, in7, vec1_m, vec0_m);                                    \
  ILVRL_H2_SH(in4, in3, vec3_m, vec2_m);                                    \
  VP9_DOT_ADD_SUB_SRARI_PCK(vec0_m, vec1_m, vec2_m, vec3_m, cnst0_m,        \
                            cnst1_m, cnst2_m, cnst3_m, in7, in0,            \
                            in4, in3);                                      \
                                                                            \
  SPLATI_H2_SH(coeff0_m, 2, 5, cnst0_m, cnst1_m);                           \
  cnst2_m = -cnst0_m;                                                       \
  ILVEV_H2_SH(cnst0_m, cnst1_m, cnst1_m, cnst2_m, cnst0_m, cnst1_m);        \
  SPLATI_H2_SH(coeff0_m, 6, 1, cnst2_m, cnst3_m);                           \
  cnst4_m = -cnst2_m;                                                       \
  ILVEV_H2_SH(cnst2_m, cnst3_m, cnst3_m, cnst4_m, cnst2_m, cnst3_m);        \
                                                                            \
  ILVRL_H2_SH(in2, in5, vec1_m, vec0_m);                                    \
  ILVRL_H2_SH(in6, in1, vec3_m, vec2_m);                                    \
                                                                            \
  VP9_DOT_ADD_SUB_SRARI_PCK(vec0_m, vec1_m, vec2_m, vec3_m, cnst0_m,        \
                            cnst1_m, cnst2_m, cnst3_m, in5, in2,            \
                            in6, in1);                                      \
  BUTTERFLY_4(in7, in0, in2, in5, s1_m, s0_m, in2, in5);                    \
  out7 = -s0_m;                                                             \
  out0 = s1_m;                                                              \
                                                                            \
  SPLATI_H4_SH(coeff1_m, 0, 4, 1, 5, cnst0_m, cnst1_m, cnst2_m, cnst3_m);   \
                                                                            \
  ILVEV_H2_SH(cnst3_m, cnst0_m, cnst1_m, cnst2_m, cnst3_m, cnst2_m);        \
  cnst0_m = __msa_ilvev_h(cnst1_m, cnst0_m);                                \
  cnst1_m = cnst0_m;                                                        \
                                                                            \
  ILVRL_H2_SH(in4, in3, vec1_m, vec0_m);                                    \
  ILVRL_H2_SH(in6, in1, vec3_m, vec2_m);                                    \
  VP9_DOT_ADD_SUB_SRARI_PCK(vec0_m, vec1_m, vec2_m, vec3_m, cnst0_m,        \
                            cnst2_m, cnst3_m, cnst1_m, out1, out6,          \
                            s0_m, s1_m);                                    \
                                                                            \
  SPLATI_H2_SH(coeff1_m, 2, 3, cnst0_m, cnst1_m);                           \
  cnst1_m = __msa_ilvev_h(cnst1_m, cnst0_m);                                \
                                                                            \
  ILVRL_H2_SH(in2, in5, vec1_m, vec0_m);                                    \
  ILVRL_H2_SH(s0_m, s1_m, vec3_m, vec2_m);                                  \
  out3 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0_m, vec1_m, cnst0_m);                \
  out4 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0_m, vec1_m, cnst1_m);                \
  out2 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec2_m, vec3_m, cnst0_m);                \
  out5 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec2_m, vec3_m, cnst1_m);                \
                                                                            \
  out1 = -out1;                                                             \
  out3 = -out3;                                                             \
  out5 = -out5;                                                             \
}

#define VP9_MADD_SHORT(m0, m1, c0, c1, res0, res1) {                \
  v4i32 madd0_m, madd1_m, madd2_m, madd3_m;                         \
  v8i16 madd_s0_m, madd_s1_m;                                       \
                                                                    \
  ILVRL_H2_SH(m1, m0, madd_s0_m, madd_s1_m);                        \
  DOTP_SH4_SW(madd_s0_m, madd_s1_m, madd_s0_m, madd_s1_m,           \
              c0, c0, c1, c1, madd0_m, madd1_m, madd2_m, madd3_m);  \
  SRARI_W4_SW(madd0_m, madd1_m, madd2_m, madd3_m, DCT_CONST_BITS);  \
  PCKEV_H2_SH(madd1_m, madd0_m, madd3_m, madd2_m, res0, res1);      \
}

#define VP9_MADD_BF(inp0, inp1, inp2, inp3, cst0, cst1, cst2, cst3,     \
                    out0, out1, out2, out3) {                           \
  v8i16 madd_s0_m, madd_s1_m, madd_s2_m, madd_s3_m;                     \
  v4i32 tmp0_m, tmp1_m, tmp2_m, tmp3_m, m4_m, m5_m;                     \
                                                                        \
  ILVRL_H2_SH(inp1, inp0, madd_s0_m, madd_s1_m);                        \
  ILVRL_H2_SH(inp3, inp2, madd_s2_m, madd_s3_m);                        \
  DOTP_SH4_SW(madd_s0_m, madd_s1_m, madd_s2_m, madd_s3_m,               \
              cst0, cst0, cst2, cst2, tmp0_m, tmp1_m, tmp2_m, tmp3_m);  \
  BUTTERFLY_4(tmp0_m, tmp1_m, tmp3_m, tmp2_m,                           \
              m4_m, m5_m, tmp3_m, tmp2_m);                              \
  SRARI_W4_SW(m4_m, m5_m, tmp2_m, tmp3_m, DCT_CONST_BITS);              \
  PCKEV_H2_SH(m5_m, m4_m, tmp3_m, tmp2_m, out0, out1);                  \
  DOTP_SH4_SW(madd_s0_m, madd_s1_m, madd_s2_m, madd_s3_m,               \
              cst1, cst1, cst3, cst3, tmp0_m, tmp1_m, tmp2_m, tmp3_m);  \
  BUTTERFLY_4(tmp0_m, tmp1_m, tmp3_m, tmp2_m,                           \
              m4_m, m5_m, tmp3_m, tmp2_m);                              \
  SRARI_W4_SW(m4_m, m5_m, tmp2_m, tmp3_m, DCT_CONST_BITS);              \
  PCKEV_H2_SH(m5_m, m4_m, tmp3_m, tmp2_m, out2, out3);                  \
}

#define VP9_LD_HADD(psrc, stride) ({                                  \
  v8i16 in0_m, in1_m, in2_m, in3_m, in4_m, in5_m, in6_m, in7_m;       \
  v4i32 vec_w_m;                                                      \
                                                                      \
  LD_SH4((psrc), stride, in0_m, in1_m, in2_m, in3_m);                 \
  ADD2(in0_m, in1_m, in2_m, in3_m, in0_m, in2_m);                     \
  LD_SH4(((psrc) + 4 * stride), stride, in4_m, in5_m, in6_m, in7_m);  \
  ADD4(in4_m, in5_m, in6_m, in7_m, in0_m, in2_m, in4_m, in6_m,        \
       in4_m, in6_m, in0_m, in4_m);                                   \
  in0_m += in4_m;                                                     \
                                                                      \
  vec_w_m = __msa_hadd_s_w(in0_m, in0_m);                             \
  HADD_SW_S32(vec_w_m);                                               \
})

#define VP9_FDCT_POSTPROC_2V_NEG_H(vec0, vec1) {  \
  v8i16 tp0_m, tp1_m;                             \
  v8i16 one_m = __msa_ldi_h(1);                   \
                                                  \
  tp0_m = __msa_clti_s_h(vec0, 0);                \
  tp1_m = __msa_clti_s_h(vec1, 0);                \
  vec0 += 1;                                      \
  vec1 += 1;                                      \
  tp0_m = one_m & tp0_m;                          \
  tp1_m = one_m & tp1_m;                          \
  vec0 += tp0_m;                                  \
  vec1 += tp1_m;                                  \
  vec0 >>= 2;                                     \
  vec1 >>= 2;                                     \
}

#define VP9_FDCT4(in0, in1, in2, in3, out0, out1, out2, out3) {     \
  v8i16 cnst0_m, cnst1_m, cnst2_m, cnst3_m;                         \
  v8i16 vec0_m, vec1_m, vec2_m, vec3_m;                             \
  v4i32 vec4_m, vec5_m, vec6_m, vec7_m;                             \
  v8i16 coeff_m = { cospi_16_64, -cospi_16_64, cospi_8_64,          \
                    cospi_24_64, -cospi_8_64, 0, 0, 0 };            \
                                                                    \
  BUTTERFLY_4(in0, in1, in2, in3, vec0_m, vec1_m, vec2_m, vec3_m);  \
  ILVR_H2_SH(vec1_m, vec0_m, vec3_m, vec2_m, vec0_m, vec2_m);       \
  SPLATI_H2_SH(coeff_m, 0, 1, cnst0_m, cnst1_m);                    \
  cnst1_m = __msa_ilvev_h(cnst1_m, cnst0_m);                        \
  vec5_m = __msa_dotp_s_w(vec0_m, cnst1_m);                         \
                                                                    \
  SPLATI_H2_SH(coeff_m, 4, 3, cnst2_m, cnst3_m);                    \
  cnst2_m = __msa_ilvev_h(cnst3_m, cnst2_m);                        \
  vec7_m = __msa_dotp_s_w(vec2_m, cnst2_m);                         \
                                                                    \
  vec4_m = __msa_dotp_s_w(vec0_m, cnst0_m);                         \
  cnst2_m = __msa_splati_h(coeff_m, 2);                             \
  cnst2_m = __msa_ilvev_h(cnst2_m, cnst3_m);                        \
  vec6_m = __msa_dotp_s_w(vec2_m, cnst2_m);                         \
                                                                    \
  SRARI_W4_SW(vec4_m, vec5_m, vec6_m, vec7_m, DCT_CONST_BITS);      \
  PCKEV_H4_SH(vec4_m, vec4_m, vec5_m, vec5_m, vec6_m, vec6_m,       \
              vec7_m, vec7_m, out0, out2, out1, out3);              \
}

#define VP9_FADST4(in0, in1, in2, in3, out0, out1, out2, out3) {  \
  v4i32 s0_m, s1_m, s2_m, s3_m, constant_m;                       \
  v4i32 in0_r_m, in1_r_m, in2_r_m, in3_r_m;                       \
                                                                  \
  UNPCK_R_SH_SW(in0, in0_r_m);                                    \
  UNPCK_R_SH_SW(in1, in1_r_m);                                    \
  UNPCK_R_SH_SW(in2, in2_r_m);                                    \
  UNPCK_R_SH_SW(in3, in3_r_m);                                    \
                                                                  \
  constant_m = __msa_fill_w(sinpi_4_9);                           \
  MUL2(in0_r_m, constant_m, in3_r_m, constant_m, s1_m, s0_m);     \
                                                                  \
  constant_m = __msa_fill_w(sinpi_1_9);                           \
  s0_m += in0_r_m * constant_m;                                   \
  s1_m -= in1_r_m * constant_m;                                   \
                                                                  \
  constant_m = __msa_fill_w(sinpi_2_9);                           \
  s0_m += in1_r_m * constant_m;                                   \
  s1_m += in3_r_m * constant_m;                                   \
                                                                  \
  s2_m = in0_r_m + in1_r_m - in3_r_m;                             \
                                                                  \
  constant_m = __msa_fill_w(sinpi_3_9);                           \
  MUL2(in2_r_m, constant_m, s2_m, constant_m, s3_m, in1_r_m);     \
                                                                  \
  in0_r_m = s0_m + s3_m;                                          \
  s2_m = s1_m - s3_m;                                             \
  s3_m = s1_m - s0_m + s3_m;                                      \
                                                                  \
  SRARI_W4_SW(in0_r_m, in1_r_m, s2_m, s3_m, DCT_CONST_BITS);      \
  PCKEV_H4_SH(in0_r_m, in0_r_m, in1_r_m, in1_r_m, s2_m, s2_m,     \
              s3_m, s3_m, out0, out1, out2, out3);                \
}

#define VP9_FDCT8(in0, in1, in2, in3, in4, in5, in6, in7,            \
                  out0, out1, out2, out3, out4, out5, out6, out7) {  \
  v8i16 s0_m, s1_m, s2_m, s3_m, s4_m, s5_m, s6_m;                    \
  v8i16 s7_m, x0_m, x1_m, x2_m, x3_m;                                \
  v8i16 coeff_m = { cospi_16_64, -cospi_16_64, cospi_8_64,           \
                    cospi_24_64, cospi_4_64, cospi_28_64,            \
                    cospi_12_64, cospi_20_64 };                      \
                                                                     \
  /* FDCT stage1 */                                                  \
  BUTTERFLY_8(in0, in1, in2, in3, in4, in5, in6, in7,                \
              s0_m, s1_m, s2_m, s3_m, s4_m, s5_m, s6_m, s7_m);       \
  BUTTERFLY_4(s0_m, s1_m, s2_m, s3_m, x0_m, x1_m, x2_m, x3_m);       \
  ILVL_H2_SH(x1_m, x0_m, x3_m, x2_m, s0_m, s2_m);                    \
  ILVR_H2_SH(x1_m, x0_m, x3_m, x2_m, s1_m, s3_m);                    \
  SPLATI_H2_SH(coeff_m, 0, 1, x0_m, x1_m);                           \
  x1_m = __msa_ilvev_h(x1_m, x0_m);                                  \
  out4 = VP9_DOT_SHIFT_RIGHT_PCK_H(s0_m, s1_m, x1_m);                \
                                                                     \
  SPLATI_H2_SH(coeff_m, 2, 3, x2_m, x3_m);                           \
  x2_m = -x2_m;                                                      \
  x2_m = __msa_ilvev_h(x3_m, x2_m);                                  \
  out6 = VP9_DOT_SHIFT_RIGHT_PCK_H(s2_m, s3_m, x2_m);                \
                                                                     \
  out0 = VP9_DOT_SHIFT_RIGHT_PCK_H(s0_m, s1_m, x0_m);                \
  x2_m = __msa_splati_h(coeff_m, 2);                                 \
  x2_m = __msa_ilvev_h(x2_m, x3_m);                                  \
  out2 = VP9_DOT_SHIFT_RIGHT_PCK_H(s2_m, s3_m, x2_m);                \
                                                                     \
  /* stage2 */                                                       \
  ILVRL_H2_SH(s5_m, s6_m, s1_m, s0_m);                               \
                                                                     \
  s6_m = VP9_DOT_SHIFT_RIGHT_PCK_H(s0_m, s1_m, x0_m);                \
  s5_m = VP9_DOT_SHIFT_RIGHT_PCK_H(s0_m, s1_m, x1_m);                \
                                                                     \
  /* stage3 */                                                       \
  BUTTERFLY_4(s4_m, s7_m, s6_m, s5_m, x0_m, x3_m, x2_m, x1_m);       \
                                                                     \
  /* stage4 */                                                       \
  ILVL_H2_SH(x3_m, x0_m, x2_m, x1_m, s4_m, s6_m);                    \
  ILVR_H2_SH(x3_m, x0_m, x2_m, x1_m, s5_m, s7_m);                    \
                                                                     \
  SPLATI_H2_SH(coeff_m, 4, 5, x0_m, x1_m);                           \
  x1_m = __msa_ilvev_h(x0_m, x1_m);                                  \
  out1 = VP9_DOT_SHIFT_RIGHT_PCK_H(s4_m, s5_m, x1_m);                \
                                                                     \
  SPLATI_H2_SH(coeff_m, 6, 7, x2_m, x3_m);                           \
  x2_m = __msa_ilvev_h(x3_m, x2_m);                                  \
  out5 = VP9_DOT_SHIFT_RIGHT_PCK_H(s6_m, s7_m, x2_m);                \
                                                                     \
  x1_m = __msa_splati_h(coeff_m, 5);                                 \
  x0_m = -x0_m;                                                      \
  x0_m = __msa_ilvev_h(x1_m, x0_m);                                  \
  out7 = VP9_DOT_SHIFT_RIGHT_PCK_H(s4_m, s5_m, x0_m);                \
                                                                     \
  x2_m = __msa_splati_h(coeff_m, 6);                                 \
  x3_m = -x3_m;                                                      \
  x2_m = __msa_ilvev_h(x2_m, x3_m);                                  \
  out3 = VP9_DOT_SHIFT_RIGHT_PCK_H(s6_m, s7_m, x2_m);                \
}

#define VP9_SRLI_AVE_S_4V_H(in0, in1, in2, in3, in4, in5, in6, in7) {    \
  v8i16 vec0_m, vec1_m, vec2_m, vec3_m, vec4_m, vec5_m, vec6_m, vec7_m;  \
                                                                         \
  SRLI_H4_SH(in0, in1, in2, in3, vec0_m, vec1_m, vec2_m, vec3_m, 15);    \
  SRLI_H4_SH(in4, in5, in6, in7, vec4_m, vec5_m, vec6_m, vec7_m, 15);    \
  AVE_SH4_SH(vec0_m, in0, vec1_m, in1, vec2_m, in2, vec3_m, in3,         \
             in0, in1, in2, in3);                                        \
  AVE_SH4_SH(vec4_m, in4, vec5_m, in5, vec6_m, in6, vec7_m, in7,         \
             in4, in5, in6, in7);                                        \
}

#define VP9_FDCT8x16_EVEN(in0, in1, in2, in3, in4, in5, in6, in7,            \
                          out0, out1, out2, out3, out4, out5, out6, out7) {  \
  v8i16 s0_m, s1_m, s2_m, s3_m, s4_m, s5_m, s6_m, s7_m;                      \
  v8i16 x0_m, x1_m, x2_m, x3_m;                                              \
  v8i16 coeff_m = { cospi_16_64, -cospi_16_64, cospi_8_64, cospi_24_64,      \
                    cospi_4_64, cospi_28_64, cospi_12_64, cospi_20_64 };     \
                                                                             \
  /* FDCT stage1 */                                                          \
  BUTTERFLY_8(in0, in1, in2, in3, in4, in5, in6, in7,                        \
              s0_m, s1_m, s2_m, s3_m, s4_m, s5_m, s6_m, s7_m);               \
  BUTTERFLY_4(s0_m, s1_m, s2_m, s3_m, x0_m, x1_m, x2_m, x3_m);               \
  ILVL_H2_SH(x1_m, x0_m, x3_m, x2_m, s0_m, s2_m);                            \
  ILVR_H2_SH(x1_m, x0_m, x3_m, x2_m, s1_m, s3_m);                            \
  SPLATI_H2_SH(coeff_m, 0, 1, x0_m, x1_m);                                   \
  x1_m = __msa_ilvev_h(x1_m, x0_m);                                          \
  out4 = VP9_DOT_SHIFT_RIGHT_PCK_H(s0_m, s1_m, x1_m);                        \
                                                                             \
  SPLATI_H2_SH(coeff_m, 2, 3, x2_m, x3_m);                                   \
  x2_m = -x2_m;                                                              \
  x2_m = __msa_ilvev_h(x3_m, x2_m);                                          \
  out6 = VP9_DOT_SHIFT_RIGHT_PCK_H(s2_m, s3_m, x2_m);                        \
                                                                             \
  out0 = VP9_DOT_SHIFT_RIGHT_PCK_H(s0_m, s1_m, x0_m);                        \
  x2_m = __msa_splati_h(coeff_m, 2);                                         \
  x2_m = __msa_ilvev_h(x2_m, x3_m);                                          \
  out2 = VP9_DOT_SHIFT_RIGHT_PCK_H(s2_m, s3_m, x2_m);                        \
                                                                             \
  /* stage2 */                                                               \
  ILVRL_H2_SH(s5_m, s6_m, s1_m, s0_m);                                       \
                                                                             \
  s6_m = VP9_DOT_SHIFT_RIGHT_PCK_H(s0_m, s1_m, x0_m);                        \
  s5_m = VP9_DOT_SHIFT_RIGHT_PCK_H(s0_m, s1_m, x1_m);                        \
                                                                             \
  /* stage3 */                                                               \
  BUTTERFLY_4(s4_m, s7_m, s6_m, s5_m, x0_m, x3_m, x2_m, x1_m);               \
                                                                             \
  /* stage4 */                                                               \
  ILVL_H2_SH(x3_m, x0_m, x2_m, x1_m, s4_m, s6_m);                            \
  ILVR_H2_SH(x3_m, x0_m, x2_m, x1_m, s5_m, s7_m);                            \
                                                                             \
  SPLATI_H2_SH(coeff_m, 4, 5, x0_m, x1_m);                                   \
  x1_m = __msa_ilvev_h(x0_m, x1_m);                                          \
  out1 = VP9_DOT_SHIFT_RIGHT_PCK_H(s4_m, s5_m, x1_m);                        \
                                                                             \
  SPLATI_H2_SH(coeff_m, 6, 7, x2_m, x3_m);                                   \
  x2_m = __msa_ilvev_h(x3_m, x2_m);                                          \
  out5 = VP9_DOT_SHIFT_RIGHT_PCK_H(s6_m, s7_m, x2_m);                        \
                                                                             \
  x1_m = __msa_splati_h(coeff_m, 5);                                         \
  x0_m = -x0_m;                                                              \
  x0_m = __msa_ilvev_h(x1_m, x0_m);                                          \
  out7 = VP9_DOT_SHIFT_RIGHT_PCK_H(s4_m, s5_m, x0_m);                        \
                                                                             \
  x2_m = __msa_splati_h(coeff_m, 6);                                         \
  x3_m = -x3_m;                                                              \
  x2_m = __msa_ilvev_h(x2_m, x3_m);                                          \
  out3 = VP9_DOT_SHIFT_RIGHT_PCK_H(s6_m, s7_m, x2_m);                        \
}

#define VP9_FDCT8x16_ODD(input0, input1, input2, input3,           \
                         input4, input5, input6, input7,           \
                         out1, out3, out5, out7,                   \
                         out9, out11, out13, out15) {              \
  v8i16 stp21_m, stp22_m, stp23_m, stp24_m, stp25_m, stp26_m;      \
  v8i16 stp30_m, stp31_m, stp32_m, stp33_m, stp34_m, stp35_m;      \
  v8i16 stp36_m, stp37_m, vec0_m, vec1_m;                          \
  v8i16 vec2_m, vec3_m, vec4_m, vec5_m, vec6_m;                    \
  v8i16 cnst0_m, cnst1_m, cnst4_m, cnst5_m;                        \
  v8i16 coeff_m = { cospi_16_64, -cospi_16_64, cospi_8_64,         \
                    cospi_24_64, -cospi_8_64, -cospi_24_64,        \
                    cospi_12_64, cospi_20_64 };                    \
  v8i16 coeff1_m = { cospi_2_64, cospi_30_64, cospi_14_64,         \
                     cospi_18_64, cospi_10_64, cospi_22_64,        \
                     cospi_6_64, cospi_26_64 };                    \
  v8i16 coeff2_m = { -cospi_2_64, -cospi_10_64, -cospi_18_64,      \
                     -cospi_26_64, 0, 0, 0, 0 };                   \
                                                                   \
  /* stp 1 */                                                      \
  ILVL_H2_SH(input2, input5, input3, input4, vec2_m, vec4_m);      \
  ILVR_H2_SH(input2, input5, input3, input4, vec3_m, vec5_m);      \
                                                                   \
  cnst4_m = __msa_splati_h(coeff_m, 0);                            \
  stp25_m = VP9_DOT_SHIFT_RIGHT_PCK_H(vec2_m, vec3_m, cnst4_m);    \
                                                                   \
  cnst5_m = __msa_splati_h(coeff_m, 1);                            \
  cnst5_m = __msa_ilvev_h(cnst5_m, cnst4_m);                       \
  stp22_m = VP9_DOT_SHIFT_RIGHT_PCK_H(vec2_m, vec3_m, cnst5_m);    \
  stp24_m = VP9_DOT_SHIFT_RIGHT_PCK_H(vec4_m, vec5_m, cnst4_m);    \
  stp23_m = VP9_DOT_SHIFT_RIGHT_PCK_H(vec4_m, vec5_m, cnst5_m);    \
                                                                   \
  /* stp2 */                                                       \
  BUTTERFLY_4(input0, input1, stp22_m, stp23_m,                    \
              stp30_m, stp31_m, stp32_m, stp33_m);                 \
  BUTTERFLY_4(input7, input6, stp25_m, stp24_m,                    \
              stp37_m, stp36_m, stp35_m, stp34_m);                 \
                                                                   \
  ILVL_H2_SH(stp36_m, stp31_m, stp35_m, stp32_m, vec2_m, vec4_m);  \
  ILVR_H2_SH(stp36_m, stp31_m, stp35_m, stp32_m, vec3_m, vec5_m);  \
                                                                   \
  SPLATI_H2_SH(coeff_m, 2, 3, cnst0_m, cnst1_m);                   \
  cnst0_m = __msa_ilvev_h(cnst0_m, cnst1_m);                       \
  stp26_m = VP9_DOT_SHIFT_RIGHT_PCK_H(vec2_m, vec3_m, cnst0_m);    \
                                                                   \
  cnst0_m = __msa_splati_h(coeff_m, 4);                            \
  cnst1_m = __msa_ilvev_h(cnst1_m, cnst0_m);                       \
  stp21_m = VP9_DOT_SHIFT_RIGHT_PCK_H(vec2_m, vec3_m, cnst1_m);    \
                                                                   \
  SPLATI_H2_SH(coeff_m, 5, 2, cnst0_m, cnst1_m);                   \
  cnst1_m = __msa_ilvev_h(cnst0_m, cnst1_m);                       \
  stp25_m = VP9_DOT_SHIFT_RIGHT_PCK_H(vec4_m, vec5_m, cnst1_m);    \
                                                                   \
  cnst0_m = __msa_splati_h(coeff_m, 3);                            \
  cnst1_m = __msa_ilvev_h(cnst1_m, cnst0_m);                       \
  stp22_m = VP9_DOT_SHIFT_RIGHT_PCK_H(vec4_m, vec5_m, cnst1_m);    \
                                                                   \
  /* stp4 */                                                       \
  BUTTERFLY_4(stp30_m, stp37_m, stp26_m, stp21_m,                  \
              vec6_m, vec2_m, vec4_m, vec5_m);                     \
  BUTTERFLY_4(stp33_m, stp34_m, stp25_m, stp22_m,                  \
              stp21_m, stp23_m, stp24_m, stp31_m);                 \
                                                                   \
  ILVRL_H2_SH(vec2_m, vec6_m, vec1_m, vec0_m);                     \
  SPLATI_H2_SH(coeff1_m, 0, 1, cnst0_m, cnst1_m);                  \
  cnst0_m = __msa_ilvev_h(cnst0_m, cnst1_m);                       \
                                                                   \
  out1 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0_m, vec1_m, cnst0_m);       \
                                                                   \
  cnst0_m = __msa_splati_h(coeff2_m, 0);                           \
  cnst0_m = __msa_ilvev_h(cnst1_m, cnst0_m);                       \
  out15 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0_m, vec1_m, cnst0_m);      \
                                                                   \
  ILVRL_H2_SH(vec4_m, vec5_m, vec1_m, vec0_m);                     \
  SPLATI_H2_SH(coeff1_m, 2, 3, cnst0_m, cnst1_m);                  \
  cnst1_m = __msa_ilvev_h(cnst1_m, cnst0_m);                       \
                                                                   \
  out9 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0_m, vec1_m, cnst1_m);       \
                                                                   \
  cnst1_m = __msa_splati_h(coeff2_m, 2);                           \
  cnst0_m = __msa_ilvev_h(cnst0_m, cnst1_m);                       \
  out7 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0_m, vec1_m, cnst0_m);       \
                                                                   \
  ILVRL_H2_SH(stp23_m, stp21_m, vec1_m, vec0_m);                   \
  SPLATI_H2_SH(coeff1_m, 4, 5, cnst0_m, cnst1_m);                  \
  cnst0_m = __msa_ilvev_h(cnst0_m, cnst1_m);                       \
  out5 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0_m, vec1_m, cnst0_m);       \
                                                                   \
  cnst0_m = __msa_splati_h(coeff2_m, 1);                           \
  cnst0_m = __msa_ilvev_h(cnst1_m, cnst0_m);                       \
  out11 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0_m, vec1_m, cnst0_m);      \
                                                                   \
  ILVRL_H2_SH(stp24_m, stp31_m, vec1_m, vec0_m);                   \
  SPLATI_H2_SH(coeff1_m, 6, 7, cnst0_m, cnst1_m);                  \
  cnst1_m = __msa_ilvev_h(cnst1_m, cnst0_m);                       \
                                                                   \
  out13 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0_m, vec1_m, cnst1_m);      \
                                                                   \
  cnst1_m = __msa_splati_h(coeff2_m, 3);                           \
  cnst0_m = __msa_ilvev_h(cnst0_m, cnst1_m);                       \
  out3 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0_m, vec1_m, cnst0_m);       \
}

#define VP9_FDCT32_POSTPROC_NEG_W(vec) {  \
  v4i32 temp_m;                           \
  v4i32 one_m = __msa_ldi_w(1);           \
                                          \
  temp_m = __msa_clti_s_w(vec, 0);        \
  vec += 1;                               \
  temp_m = one_m & temp_m;                \
  vec += temp_m;                          \
  vec >>= 2;                              \
}

#define VP9_FDCT32_POSTPROC_2V_POS_H(vec0, vec1) {  \
  v8i16 tp0_m, tp1_m;                               \
  v8i16 one = __msa_ldi_h(1);                       \
                                                    \
  tp0_m = __msa_clei_s_h(vec0, 0);                  \
  tp1_m = __msa_clei_s_h(vec1, 0);                  \
  tp0_m = (v8i16)__msa_xori_b((v16u8)tp0_m, 255);   \
  tp1_m = (v8i16)__msa_xori_b((v16u8)tp1_m, 255);   \
  vec0 += 1;                                        \
  vec1 += 1;                                        \
  tp0_m = one & tp0_m;                              \
  tp1_m = one & tp1_m;                              \
  vec0 += tp0_m;                                    \
  vec1 += tp1_m;                                    \
  vec0 >>= 2;                                       \
  vec1 >>= 2;                                       \
}

#define VP9_DOTP_CONST_PAIR_W(reg0_left, reg1_left, reg0_right,  \
                              reg1_right, const0, const1,        \
                              out0, out1, out2, out3) {          \
  v4i32 s0_m, s1_m, s2_m, s3_m, s4_m, s5_m, s6_m, s7_m;          \
  v2i64 tp0_m, tp1_m, tp2_m, tp3_m;                              \
  v4i32 k0_m = __msa_fill_w((int32_t) const0);                   \
                                                                 \
  s0_m = __msa_fill_w((int32_t) const1);                         \
  k0_m = __msa_ilvev_w(s0_m, k0_m);                              \
                                                                 \
  ILVRL_W2_SW(-reg1_left, reg0_left, s1_m, s0_m);                \
  ILVRL_W2_SW(reg0_left, reg1_left, s3_m, s2_m);                 \
  ILVRL_W2_SW(-reg1_right, reg0_right, s5_m, s4_m);              \
  ILVRL_W2_SW(reg0_right, reg1_right, s7_m, s6_m);               \
                                                                 \
  DOTP_SW2_SD(s0_m, s1_m, k0_m, k0_m, tp0_m, tp1_m);             \
  DOTP_SW2_SD(s4_m, s5_m, k0_m, k0_m, tp2_m, tp3_m);             \
  tp0_m = __msa_srari_d(tp0_m, DCT_CONST_BITS);                  \
  tp1_m = __msa_srari_d(tp1_m, DCT_CONST_BITS);                  \
  tp2_m = __msa_srari_d(tp2_m, DCT_CONST_BITS);                  \
  tp3_m = __msa_srari_d(tp3_m, DCT_CONST_BITS);                  \
  out0 = __msa_pckev_w((v4i32)tp0_m, (v4i32)tp1_m);              \
  out1 = __msa_pckev_w((v4i32)tp2_m, (v4i32)tp3_m);              \
                                                                 \
  DOTP_SW2_SD(s2_m, s3_m, k0_m, k0_m, tp0_m, tp1_m);             \
  DOTP_SW2_SD(s6_m, s7_m, k0_m, k0_m, tp2_m, tp3_m);             \
  tp0_m = __msa_srari_d(tp0_m, DCT_CONST_BITS);                  \
  tp1_m = __msa_srari_d(tp1_m, DCT_CONST_BITS);                  \
  tp2_m = __msa_srari_d(tp2_m, DCT_CONST_BITS);                  \
  tp3_m = __msa_srari_d(tp3_m, DCT_CONST_BITS);                  \
  out2 = __msa_pckev_w((v4i32)tp0_m, (v4i32)tp1_m);              \
  out3 = __msa_pckev_w((v4i32)tp2_m, (v4i32)tp3_m);              \
}
#endif  /* VP9_ENCODER_MIPS_MSA_VP9_FDCT_MSA_H_ */
