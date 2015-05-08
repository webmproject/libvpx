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

#define VP9_SET_CONST_PAIR(mask_h, idx1_h, idx2_h) ({  \
  v8i16 c0_m, c1_m;                                    \
                                                       \
  c0_m = __msa_splati_h((mask_h), (idx1_h));           \
  c1_m = __msa_splati_h((mask_h), (idx2_h));           \
  c0_m = __msa_ilvev_h(c1_m, c0_m);                    \
                                                       \
  c0_m;                                                \
})

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

/* multiply and add macro */
#define VP9_MADD(inp0, inp1, inp2, inp3,                      \
                 cst0, cst1, cst2, cst3,                      \
                 out0, out1, out2, out3) {                    \
  v8i16 madd_s0_m, madd_s1_m, madd_s2_m, madd_s3_m;           \
  v4i32 tmp0_m, tmp1_m, tmp2_m, tmp3_m;                       \
                                                              \
  ILV_H_LRLR_SH(inp0, inp1, inp2, inp3,                       \
                madd_s0_m, madd_s1_m, madd_s2_m, madd_s3_m);  \
                                                              \
  DOTP_S_W_4VECS_SW(madd_s1_m, cst0, madd_s0_m, cst0,         \
                    madd_s1_m, cst1, madd_s0_m, cst1,         \
                    tmp0_m, tmp1_m, tmp2_m, tmp3_m);          \
                                                              \
  SRARI_W_4VECS_SW(tmp0_m, tmp1_m, tmp2_m, tmp3_m,            \
                   tmp0_m, tmp1_m, tmp2_m, tmp3_m,            \
                   DCT_CONST_BITS);                           \
                                                              \
  PCKEV_H_2VECS_SH(tmp1_m, tmp0_m, tmp3_m, tmp2_m,            \
                   out0, out1);                               \
                                                              \
  DOTP_S_W_4VECS_SW(madd_s3_m, cst2, madd_s2_m, cst2,         \
                    madd_s3_m, cst3, madd_s2_m, cst3,         \
                    tmp0_m, tmp1_m, tmp2_m, tmp3_m);          \
                                                              \
  SRARI_W_4VECS_SW(tmp0_m, tmp1_m, tmp2_m, tmp3_m,            \
                   tmp0_m, tmp1_m, tmp2_m, tmp3_m,            \
                   DCT_CONST_BITS);                           \
                                                              \
  PCKEV_H_2VECS_SH(tmp1_m, tmp0_m, tmp3_m, tmp2_m,            \
                   out2, out3);                               \
}

/* idct 8x8 macro */
#define VP9_IDCT8x8_1D_ODD(in1, in3, in5, in7,        \
                           k0, k1, k2, k3, mask,      \
                           out0, out1, out2, out3) {  \
  v8i16 res0_m, res1_m, res2_m, res3_m;               \
  v4i32 tmp0_m, tmp1_m, tmp2_m, tmp3_m;               \
                                                      \
  VP9_MADD(in1, in7, in3, in5, k0, k1, k2, k3,        \
           in1, in7, in3, in5);                       \
                                                      \
  res0_m = in1 - in3;                                 \
  res1_m = in7 - in5;                                 \
                                                      \
  k0 = VP9_SET_CONST_PAIR(mask, 4, 7);                \
  k1 = __msa_splati_h(mask, 4);                       \
                                                      \
  res2_m = __msa_ilvr_h(res0_m, res1_m);              \
  res3_m = __msa_ilvl_h(res0_m, res1_m);              \
                                                      \
  DOTP_S_W_4VECS_SW(res2_m, k0, res3_m, k0,           \
                    res2_m, k1, res3_m, k1,           \
                    tmp0_m, tmp1_m, tmp2_m, tmp3_m);  \
                                                      \
  SRARI_W_4VECS_SW(tmp0_m, tmp1_m, tmp2_m, tmp3_m,    \
                   tmp0_m, tmp1_m, tmp2_m, tmp3_m,    \
                   DCT_CONST_BITS);                   \
  out0 = in1 + in3;                                   \
  PCKEV_H_2VECS_SH(tmp1_m, tmp0_m, tmp3_m, tmp2_m,    \
                   out1, out2);                       \
  out3 = in7 + in5;                                   \
}

#define VP9_IDCT8x8_1D_EVEN(in0, in2, in4, in6,        \
                            k0, k1, k2, k3,            \
                            out0, out1, out2, out3) {  \
  k2 = SET_COSPI_PAIR(cospi_24_64, -cospi_8_64);       \
  k3 = SET_COSPI_PAIR(cospi_8_64, cospi_24_64);        \
                                                       \
  VP9_MADD(in0, in4, in2, in6, k1, k0, k2, k3,         \
           in0, in4, in2, in6);                        \
                                                       \
  out0 = in0 + in6;                                    \
  out1 = in4 + in2;                                    \
  out2 = in4 - in2;                                    \
  out3 = in0 - in6;                                    \
}

#define VP9_IDCT8x8_1D(in0, in1, in2, in3, in4, in5, in6, in7,            \
                       out0, out1, out2, out3, out4, out5, out6, out7) {  \
  v8i16 res0_m, res1_m, res2_m, res3_m, res4_m, res5_m, res6_m, res7_m;   \
  v8i16 k0_m, k1_m, k2_m, k3_m;                                           \
  v8i16 mask_m = { cospi_28_64, cospi_4_64, cospi_20_64, cospi_12_64,     \
    cospi_16_64, -cospi_4_64, -cospi_20_64, -cospi_16_64                  \
  };                                                                      \
                                                                          \
  k0_m = VP9_SET_CONST_PAIR(mask_m, 0, 5);                                \
  k1_m = VP9_SET_CONST_PAIR(mask_m, 1, 0);                                \
  k2_m = VP9_SET_CONST_PAIR(mask_m, 6, 3);                                \
  k3_m = VP9_SET_CONST_PAIR(mask_m, 3, 2);                                \
                                                                          \
  VP9_IDCT8x8_1D_ODD(in1, in3, in5, in7, k0_m, k1_m, k2_m, k3_m, mask_m,  \
                     res4_m, res5_m, res6_m, res7_m);                     \
                                                                          \
  VP9_IDCT8x8_1D_EVEN(in0, in2, in4, in6, k0_m, k1_m, k2_m, k3_m,         \
                      res0_m, res1_m, res2_m, res3_m);                    \
                                                                          \
  BUTTERFLY_8(res0_m, res1_m, res2_m, res3_m,                             \
              res4_m, res5_m, res6_m, res7_m,                             \
              out0, out1, out2, out3,                                     \
              out4, out5, out6, out7);                                    \
}

#define DOT_ADD_SUB_SRARI_PCK(in0, in1, in2, in3, in4, in5, in6, in7,  \
                              dst0, dst1, dst2, dst3) {                \
  v4i32 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;    \
                                                                       \
  tmp0 = __msa_dotp_s_w((in0), (in4));                                 \
  tmp2 = __msa_dotp_s_w((in1), (in4));                                 \
  tmp3 = __msa_dotp_s_w((in0), (in5));                                 \
  tmp4 = __msa_dotp_s_w((in1), (in5));                                 \
  tmp5 = __msa_dotp_s_w((in2), (in6));                                 \
  tmp6 = __msa_dotp_s_w((in3), (in6));                                 \
  tmp7 = __msa_dotp_s_w((in2), (in7));                                 \
  tmp8 = __msa_dotp_s_w((in3), (in7));                                 \
                                                                       \
  BUTTERFLY_4(tmp0, tmp3, tmp7, tmp5, tmp1, tmp9, tmp7, tmp5);         \
  BUTTERFLY_4(tmp2, tmp4, tmp8, tmp6, tmp3, tmp0, tmp4, tmp2);         \
                                                                       \
  SRARI_W_4VECS_SW(tmp1, tmp9, tmp7, tmp5, tmp1, tmp9, tmp7, tmp5,     \
                   DCT_CONST_BITS);                                    \
  SRARI_W_4VECS_SW(tmp3, tmp0, tmp4, tmp2, tmp3, tmp0, tmp4, tmp2,     \
                   DCT_CONST_BITS);                                    \
                                                                       \
  PCKEV_H_4VECS_SH(tmp1, tmp3, tmp9, tmp0, tmp7, tmp4, tmp5, tmp2,     \
                   dst0, dst1, dst2, dst3);                            \
}

#define DOT_SHIFT_RIGHT_PCK_H(in0, in1, in2) ({       \
  v8i16 dst_m;                                        \
  v4i32 tp0_m, tp1_m;                                 \
                                                      \
  tp1_m = __msa_dotp_s_w((in0), (in2));               \
  tp0_m = __msa_dotp_s_w((in1), (in2));               \
  tp1_m = __msa_srari_w(tp1_m, DCT_CONST_BITS);       \
  tp0_m = __msa_srari_w(tp0_m, DCT_CONST_BITS);       \
  dst_m = __msa_pckev_h((v8i16)tp1_m, (v8i16)tp0_m);  \
                                                      \
  dst_m;                                              \
})

#define VP9_ADST8_ROW(in0, in1, in2, in3, in4, in5, in6, in7,            \
                      out0, out1, out2, out3, out4, out5, out6, out7) {  \
  v8i16 const0_m, const1_m, const2_m, const3_m, const4_m;                \
  v8i16 temp0_m, temp1_m, temp2_m, temp3_m, s0_m, s1_m;                  \
  v8i16 coeff0_m = { cospi_2_64, cospi_6_64, cospi_10_64,                \
    cospi_14_64, cospi_18_64, cospi_22_64, cospi_26_64, cospi_30_64      \
  };                                                                     \
  v8i16 coeff1_m = { cospi_8_64, -cospi_8_64, cospi_16_64,               \
    -cospi_16_64, cospi_24_64, -cospi_24_64, 0, 0                        \
  };                                                                     \
                                                                         \
  const0_m = __msa_splati_h(coeff0_m, 0);                                \
  const1_m = __msa_splati_h(coeff0_m, 7);                                \
  const2_m = -const0_m;                                                  \
  const0_m = __msa_ilvev_h(const1_m, const0_m);                          \
  const1_m = __msa_ilvev_h(const2_m, const1_m);                          \
  const2_m = __msa_splati_h(coeff0_m, 4);                                \
  const3_m = __msa_splati_h(coeff0_m, 3);                                \
  const4_m = -const2_m;                                                  \
  const2_m = __msa_ilvev_h(const3_m, const2_m);                          \
  const3_m = __msa_ilvev_h(const4_m, const3_m);                          \
                                                                         \
  ILV_H_LRLR_SH(in7, in0, in3, in4,                                      \
                temp0_m, temp1_m, temp2_m, temp3_m);                     \
                                                                         \
  DOT_ADD_SUB_SRARI_PCK(temp0_m, temp1_m, temp2_m, temp3_m,              \
                        const0_m, const1_m, const2_m, const3_m,          \
                        in7, in0, in4, in3);                             \
                                                                         \
  const0_m = __msa_splati_h(coeff0_m, 2);                                \
  const1_m = __msa_splati_h(coeff0_m, 5);                                \
  const2_m = -const0_m;                                                  \
  const0_m = __msa_ilvev_h(const1_m, const0_m);                          \
  const1_m = __msa_ilvev_h(const2_m, const1_m);                          \
  const2_m = __msa_splati_h(coeff0_m, 6);                                \
  const3_m = __msa_splati_h(coeff0_m, 1);                                \
  const4_m = -const2_m;                                                  \
  const2_m = __msa_ilvev_h(const3_m, const2_m);                          \
  const3_m = __msa_ilvev_h(const4_m, const3_m);                          \
                                                                         \
  ILV_H_LRLR_SH(in5, in2, in1, in6,                                      \
                temp0_m, temp1_m, temp2_m, temp3_m);                     \
                                                                         \
  DOT_ADD_SUB_SRARI_PCK(temp0_m, temp1_m, temp2_m, temp3_m,              \
                        const0_m, const1_m, const2_m, const3_m,          \
                        in5, in2, in6, in1);                             \
                                                                         \
  BUTTERFLY_4(in7, in0, in2, in5, s1_m, s0_m, in2, in5);                 \
  out7 = -s0_m;                                                          \
  out0 = s1_m;                                                           \
                                                                         \
  SPLATI_H_4VECS_SH(coeff1_m, 0, 4, 1, 5,                                \
                    const0_m, const1_m, const2_m, const3_m);             \
                                                                         \
  const3_m = __msa_ilvev_h(const0_m, const3_m);                          \
  const2_m = __msa_ilvev_h(const2_m, const1_m);                          \
  const0_m = __msa_ilvev_h(const1_m, const0_m);                          \
  const1_m = const0_m;                                                   \
                                                                         \
  ILV_H_LRLR_SH(in3, in4, in1, in6,                                      \
                temp0_m, temp1_m, temp2_m, temp3_m);                     \
                                                                         \
  DOT_ADD_SUB_SRARI_PCK(temp0_m, temp1_m, temp2_m, temp3_m,              \
                        const0_m, const2_m, const3_m, const1_m,          \
                        out1, out6, s0_m, s1_m);                         \
                                                                         \
  const0_m = __msa_splati_h(coeff1_m, 2);                                \
  const1_m = __msa_splati_h(coeff1_m, 3);                                \
  const1_m = __msa_ilvev_h(const1_m, const0_m);                          \
                                                                         \
  ILV_H_LRLR_SH(in5, in2, s1_m, s0_m,                                    \
             temp0_m, temp1_m, temp2_m, temp3_m);                        \
                                                                         \
  out3 = DOT_SHIFT_RIGHT_PCK_H(temp0_m, temp1_m, const0_m);              \
  out4 = DOT_SHIFT_RIGHT_PCK_H(temp0_m, temp1_m, const1_m);              \
  out2 = DOT_SHIFT_RIGHT_PCK_H(temp2_m, temp3_m, const0_m);              \
  out5 = DOT_SHIFT_RIGHT_PCK_H(temp2_m, temp3_m, const1_m);              \
                                                                         \
  out1 = -out1;                                                          \
  out3 = -out3;                                                          \
  out5 = -out5;                                                          \
}

#define VP9_ADST8(in0, in1, in2, in3, in4, in5, in6, in7,            \
                  out0, out1, out2, out3, out4, out5, out6, out7) {  \
  v8i16 const0_m, const1_m, const2_m, const3_m, const4_m;            \
  v8i16 temp0_m, temp1_m, temp2_m, temp3_m, s0_m, s1_m;              \
                                                                     \
  const0_m = __msa_fill_h(cospi_2_64);                               \
  const1_m = __msa_fill_h(cospi_30_64);                              \
  const2_m = -const0_m;                                              \
  const0_m = __msa_ilvev_h(const1_m, const0_m);                      \
  const1_m = __msa_ilvev_h(const2_m, const1_m);                      \
  const2_m = __msa_fill_h(cospi_18_64);                              \
  const3_m = __msa_fill_h(cospi_14_64);                              \
  const4_m = -const2_m;                                              \
  const2_m = __msa_ilvev_h(const3_m, const2_m);                      \
  const3_m = __msa_ilvev_h(const4_m, const3_m);                      \
                                                                     \
  ILV_H_LRLR_SH(in7, in0, in3, in4,                                  \
                temp0_m, temp1_m, temp2_m, temp3_m);                 \
                                                                     \
  DOT_ADD_SUB_SRARI_PCK(temp0_m, temp1_m, temp2_m, temp3_m,          \
                        const0_m, const1_m, const2_m, const3_m,      \
                        in7, in0, in4, in3);                         \
                                                                     \
  const0_m = __msa_fill_h(cospi_10_64);                              \
  const1_m = __msa_fill_h(cospi_22_64);                              \
  const2_m = -const0_m;                                              \
  const0_m = __msa_ilvev_h(const1_m, const0_m);                      \
  const1_m = __msa_ilvev_h(const2_m, const1_m);                      \
  const2_m = __msa_fill_h(cospi_26_64);                              \
  const3_m = __msa_fill_h(cospi_6_64);                               \
  const4_m = -const2_m;                                              \
  const2_m = __msa_ilvev_h(const3_m, const2_m);                      \
  const3_m = __msa_ilvev_h(const4_m, const3_m);                      \
                                                                     \
  ILV_H_LRLR_SH(in5, in2, in1, in6,                                  \
                temp0_m, temp1_m, temp2_m, temp3_m);                 \
                                                                     \
  DOT_ADD_SUB_SRARI_PCK(temp0_m, temp1_m, temp2_m, temp3_m,          \
                        const0_m, const1_m, const2_m, const3_m,      \
                        in5, in2, in6, in1);                         \
                                                                     \
  BUTTERFLY_4(in7, in0, in2, in5, s1_m, s0_m, in2, in5);             \
  out7 = -s0_m;                                                      \
  out0 = s1_m;                                                       \
                                                                     \
  const1_m = __msa_fill_h(cospi_24_64);                              \
  const0_m = __msa_fill_h(cospi_8_64);                               \
  const3_m = -const1_m;                                              \
  const2_m = -const0_m;                                              \
                                                                     \
  const3_m = __msa_ilvev_h(const0_m, const3_m);                      \
  const2_m = __msa_ilvev_h(const2_m, const1_m);                      \
  const0_m = __msa_ilvev_h(const1_m, const0_m);                      \
  const1_m = const0_m;                                               \
                                                                     \
  ILV_H_LRLR_SH(in3, in4, in1, in6,                                  \
                temp0_m, temp1_m, temp2_m, temp3_m);                 \
                                                                     \
  DOT_ADD_SUB_SRARI_PCK(temp0_m, temp1_m, temp2_m, temp3_m,          \
                        const0_m, const2_m, const3_m, const1_m,      \
                        out1, out6, s0_m, s1_m);                     \
                                                                     \
  const0_m = __msa_fill_h(cospi_16_64);                              \
  const1_m = -const0_m;                                              \
  const1_m = __msa_ilvev_h(const1_m, const0_m);                      \
                                                                     \
  ILV_H_LRLR_SH(in5, in2, s1_m, s0_m,                                \
                temp0_m, temp1_m, temp2_m, temp3_m);                 \
                                                                     \
  out3 = DOT_SHIFT_RIGHT_PCK_H(temp0_m, temp1_m, const0_m);          \
  out4 = DOT_SHIFT_RIGHT_PCK_H(temp0_m, temp1_m, const1_m);          \
  out2 = DOT_SHIFT_RIGHT_PCK_H(temp2_m, temp3_m, const0_m);          \
  out5 = DOT_SHIFT_RIGHT_PCK_H(temp2_m, temp3_m, const1_m);          \
                                                                     \
  out1 = -out1;                                                      \
  out3 = -out3;                                                      \
  out5 = -out5;                                                      \
}

void vp9_idct8x8_64_add_msa(const int16_t *input, uint8_t *dest,
                            int32_t dest_stride) {
  v8i16 in0, in1, in2, in3, in4, in5, in6, in7;

  /* load vector elements of 8x8 block */
  LOAD_8VECS_SH(input, 8, in0, in1, in2, in3, in4, in5, in6, in7);

  /* rows transform */
  TRANSPOSE8x8_H_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                    in0, in1, in2, in3, in4, in5, in6, in7);

  /* 1D idct8x8 */
  VP9_IDCT8x8_1D(in0, in1, in2, in3, in4, in5, in6, in7,
                 in0, in1, in2, in3, in4, in5, in6, in7);

  /* columns transform */
  TRANSPOSE8x8_H_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                    in0, in1, in2, in3, in4, in5, in6, in7);

  /* 1D idct8x8 */
  VP9_IDCT8x8_1D(in0, in1, in2, in3, in4, in5, in6, in7,
                 in0, in1, in2, in3, in4, in5, in6, in7);

  /* final rounding (add 2^4, divide by 2^5) and shift */
  SRARI_H_4VECS_SH(in0, in1, in2, in3, in0, in1, in2, in3, 5);
  SRARI_H_4VECS_SH(in4, in5, in6, in7, in4, in5, in6, in7, 5);

  /* add block and store 8x8 */
  VP9_ADDBLK_CLIP_AND_STORE_8_BYTES_4(dest, dest_stride, in0, in1, in2, in3);
  dest += (4 * dest_stride);
  VP9_ADDBLK_CLIP_AND_STORE_8_BYTES_4(dest, dest_stride, in4, in5, in6, in7);
}

void vp9_idct8x8_12_add_msa(const int16_t *input, uint8_t *dest,
                            int32_t dest_stride) {
  v8i16 in0, in1, in2, in3, in4, in5, in6, in7;
  v8i16 s0, s1, s2, s3, s4, s5, s6, s7;
  v8i16 k0, k1, k2, k3, m0, m1, m2, m3;
  v4i32 tmp0, tmp1, tmp2, tmp3;
  v8i16 zero = { 0 };

  /* load vector elements of 8x8 block */
  LOAD_8VECS_SH(input, 8, in0, in1, in2, in3, in4, in5, in6, in7);

  TRANSPOSE8X4_H(in0, in1, in2, in3, in0, in1, in2, in3);

  /* stage1 */
  s0 = __msa_ilvl_h(in3, in0);
  s1 = __msa_ilvl_h(in2, in1);

  k0 = SET_COSPI_PAIR(cospi_28_64, -cospi_4_64);
  k1 = SET_COSPI_PAIR(cospi_4_64, cospi_28_64);
  k2 = SET_COSPI_PAIR(-cospi_20_64, cospi_12_64);
  k3 = SET_COSPI_PAIR(cospi_12_64, cospi_20_64);
  DOTP_S_W_4VECS_SW(s0, k0, s0, k1, s1, k2, s1, k3, tmp0, tmp1, tmp2, tmp3);

  SRARI_W_4VECS_SW(tmp0, tmp1, tmp2, tmp3,
                   tmp0, tmp1, tmp2, tmp3, DCT_CONST_BITS);

  PCKEV_H_2VECS_SH(zero, tmp0, zero, tmp1, s0, s1);
  PCKEV_H_2VECS_SH(zero, tmp2, zero, tmp3, s2, s3);

  BUTTERFLY_4(s0, s1, s3, s2, s4, s7, s6, s5);

  /* stage2 */
  s0 = __msa_ilvr_h(in2, in0);
  s1 = __msa_ilvr_h(in3, in1);

  k0 = SET_COSPI_PAIR(cospi_16_64, cospi_16_64);
  k1 = SET_COSPI_PAIR(cospi_16_64, -cospi_16_64);
  k2 = SET_COSPI_PAIR(cospi_24_64, -cospi_8_64);
  k3 = SET_COSPI_PAIR(cospi_8_64, cospi_24_64);
  DOTP_S_W_4VECS_SW(s0, k0, s0, k1, s1, k2, s1, k3, tmp0, tmp1, tmp2, tmp3);

  SRARI_W_4VECS_SW(tmp0, tmp1, tmp2, tmp3,
                   tmp0, tmp1, tmp2, tmp3, DCT_CONST_BITS);

  PCKEV_H_2VECS_SH(zero, tmp0, zero, tmp1, s0, s1);
  PCKEV_H_2VECS_SH(zero, tmp2, zero, tmp3, s2, s3);

  BUTTERFLY_4(s0, s1, s2, s3, m0, m1, m2, m3);

  /* stage3 */
  s0 = __msa_ilvr_h(s6, s5);

  k1 = SET_COSPI_PAIR(-cospi_16_64, cospi_16_64);
  tmp0 = __msa_dotp_s_w(s0, k1);
  tmp1 = __msa_dotp_s_w(s0, k0);

  tmp0 = __msa_srari_w(tmp0, DCT_CONST_BITS);
  tmp1 = __msa_srari_w(tmp1, DCT_CONST_BITS);

  PCKEV_H_2VECS_SH(zero, tmp0, zero, tmp1, s2, s3);

  /* stage4 */
  BUTTERFLY_8(m0, m1, m2, m3, s4, s2, s3, s7,
              in0, in1, in2, in3, in4, in5, in6, in7);

  TRANSPOSE4X8_H(in0, in1, in2, in3, in4, in5, in6, in7,
                 in0, in1, in2, in3, in4, in5, in6, in7);

  VP9_IDCT8x8_1D(in0, in1, in2, in3, in4, in5, in6, in7,
                 in0, in1, in2, in3, in4, in5, in6, in7);

  /* final rounding (add 2^4, divide by 2^5) and shift */
  SRARI_H_4VECS_SH(in0, in1, in2, in3, in0, in1, in2, in3, 5);
  SRARI_H_4VECS_SH(in4, in5, in6, in7, in4, in5, in6, in7, 5);

  /* add block and store 8x8 */
  VP9_ADDBLK_CLIP_AND_STORE_8_BYTES_4(dest, dest_stride, in0, in1, in2, in3);
  dest += (4 * dest_stride);
  VP9_ADDBLK_CLIP_AND_STORE_8_BYTES_4(dest, dest_stride, in4, in5, in6, in7);
}

void vp9_idct8x8_1_add_msa(const int16_t *input, uint8_t *dest,
                           int32_t dest_stride) {
  int16_t out;
  int32_t const1;
  v8i16 const2;

  out = dct_const_round_shift(input[0] * cospi_16_64);
  out = dct_const_round_shift(out * cospi_16_64);
  const1 = ROUND_POWER_OF_TWO(out, 5);
  const2 = __msa_fill_h(const1);

  VP9_ADDBLK_CLIP_AND_STORE_8_BYTES_4(dest, dest_stride,
                                      const2, const2, const2, const2);
  dest += (4 * dest_stride);
  VP9_ADDBLK_CLIP_AND_STORE_8_BYTES_4(dest, dest_stride,
                                      const2, const2, const2, const2);
}

void vp9_iht8x8_64_add_msa(const int16_t *input, uint8_t *dest,
                           int32_t dest_stride, int32_t tx_type) {
  v8i16 in0, in1, in2, in3, in4, in5, in6, in7;

  /* load vector elements of 8x8 block */
  LOAD_8VECS_SH(input, 8, in0, in1, in2, in3, in4, in5, in6, in7);

  TRANSPOSE8x8_H_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                    in0, in1, in2, in3, in4, in5, in6, in7);

  switch (tx_type) {
    case DCT_DCT:
      /* DCT in horizontal */
      VP9_IDCT8x8_1D(in0, in1, in2, in3, in4, in5, in6, in7,
                     in0, in1, in2, in3, in4, in5, in6, in7);

      /* DCT in vertical */
      TRANSPOSE8x8_H_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                        in0, in1, in2, in3, in4, in5, in6, in7);
      VP9_IDCT8x8_1D(in0, in1, in2, in3, in4, in5, in6, in7,
                     in0, in1, in2, in3, in4, in5, in6, in7);
      break;
    case ADST_DCT:
      /* DCT in horizontal */
      VP9_IDCT8x8_1D(in0, in1, in2, in3, in4, in5, in6, in7,
                     in0, in1, in2, in3, in4, in5, in6, in7);

      /* ADST in vertical */
      TRANSPOSE8x8_H_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                        in0, in1, in2, in3, in4, in5, in6, in7);
      VP9_ADST8(in0, in1, in2, in3, in4, in5, in6, in7,
                in0, in1, in2, in3, in4, in5, in6, in7);
      break;
    case DCT_ADST:
      /* ADST in horizontal */
      VP9_ADST8_ROW(in0, in1, in2, in3, in4, in5, in6, in7,
                    in0, in1, in2, in3, in4, in5, in6, in7);

      /* DCT in vertical */
      TRANSPOSE8x8_H_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                        in0, in1, in2, in3, in4, in5, in6, in7);
      VP9_IDCT8x8_1D(in0, in1, in2, in3, in4, in5, in6, in7,
                     in0, in1, in2, in3, in4, in5, in6, in7);
      break;
    case ADST_ADST:
      /* ADST in horizontal */
      VP9_ADST8(in0, in1, in2, in3, in4, in5, in6, in7,
                in0, in1, in2, in3, in4, in5, in6, in7);

      /* ADST in vertical */
      TRANSPOSE8x8_H_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                        in0, in1, in2, in3, in4, in5, in6, in7);
      VP9_ADST8(in0, in1, in2, in3, in4, in5, in6, in7,
                in0, in1, in2, in3, in4, in5, in6, in7);
      break;
    default:
      assert(0);
      break;
  }

  /* final rounding (add 2^4, divide by 2^5) and shift */
  SRARI_H_4VECS_SH(in0, in1, in2, in3, in0, in1, in2, in3, 5);
  SRARI_H_4VECS_SH(in4, in5, in6, in7, in4, in5, in6, in7, 5);

  /* add block and store 8x8 */
  VP9_ADDBLK_CLIP_AND_STORE_8_BYTES_4(dest, dest_stride, in0, in1, in2, in3);
  dest += (4 * dest_stride);
  VP9_ADDBLK_CLIP_AND_STORE_8_BYTES_4(dest, dest_stride, in4, in5, in6, in7);
}
