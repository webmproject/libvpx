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

#include "vp9/common/vp9_idct.h"
#include "vpx_dsp/mips/fwd_txfm_msa.h"
#include "vpx_dsp/mips/txfm_macros_msa.h"
#include "vpx_ports/mem.h"

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
  DOT_ADD_SUB_SRARI_PCK(vec0_m, vec1_m, vec2_m, vec3_m, cnst0_m,            \
                        cnst1_m, cnst2_m, cnst3_m, in7, in0,                \
                        in4, in3);                                          \
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
  DOT_ADD_SUB_SRARI_PCK(vec0_m, vec1_m, vec2_m, vec3_m, cnst0_m,            \
                        cnst1_m, cnst2_m, cnst3_m, in5, in2,                \
                        in6, in1);                                          \
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
  DOT_ADD_SUB_SRARI_PCK(vec0_m, vec1_m, vec2_m, vec3_m, cnst0_m,            \
                        cnst2_m, cnst3_m, cnst1_m, out1, out6,              \
                        s0_m, s1_m);                                        \
                                                                            \
  SPLATI_H2_SH(coeff1_m, 2, 3, cnst0_m, cnst1_m);                           \
  cnst1_m = __msa_ilvev_h(cnst1_m, cnst0_m);                                \
                                                                            \
  ILVRL_H2_SH(in2, in5, vec1_m, vec0_m);                                    \
  ILVRL_H2_SH(s0_m, s1_m, vec3_m, vec2_m);                                  \
  out3 = DOT_SHIFT_RIGHT_PCK_H(vec0_m, vec1_m, cnst0_m);                    \
  out4 = DOT_SHIFT_RIGHT_PCK_H(vec0_m, vec1_m, cnst1_m);                    \
  out2 = DOT_SHIFT_RIGHT_PCK_H(vec2_m, vec3_m, cnst0_m);                    \
  out5 = DOT_SHIFT_RIGHT_PCK_H(vec2_m, vec3_m, cnst1_m);                    \
                                                                            \
  out1 = -out1;                                                             \
  out3 = -out3;                                                             \
  out5 = -out5;                                                             \
}

#define LD_HADD(psrc, stride) ({                                      \
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

#define FDCT_POSTPROC_2V_NEG_H(vec0, vec1) {      \
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

#define FDCT32_POSTPROC_NEG_W(vec) {      \
  v4i32 temp_m;                           \
  v4i32 one_m = __msa_ldi_w(1);           \
                                          \
  temp_m = __msa_clti_s_w(vec, 0);        \
  vec += 1;                               \
  temp_m = one_m & temp_m;                \
  vec += temp_m;                          \
  vec >>= 2;                              \
}

#define FDCT32_POSTPROC_2V_POS_H(vec0, vec1) {      \
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

#define DOTP_CONST_PAIR_W(reg0_left, reg1_left, reg0_right,      \
                          reg1_right, const0, const1,            \
                          out0, out1, out2, out3) {              \
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
