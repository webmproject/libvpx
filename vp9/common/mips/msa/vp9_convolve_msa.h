/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_MIPS_MSA_VP9_CONVOLVE_MSA_H_
#define VP9_COMMON_MIPS_MSA_VP9_CONVOLVE_MSA_H_

#include "vp9/common/vp9_filter.h"
#include "vp9/common/mips/msa/vp9_macros_msa.h"

extern uint8_t mc_filt_mask_arr[16 * 3];

#define HORIZ_8TAP_FILT(src, mask0, mask1, mask2, mask3,                   \
                        filt_h0, filt_h1, filt_h2, filt_h3) ({             \
  v8i16 vec0, vec1, vec2, vec3, horiz_out;                                 \
                                                                           \
  vec0 = (v8i16)__msa_vshf_b((v16i8)(mask0), (v16i8)(src), (v16i8)(src));  \
  vec0 = __msa_dotp_s_h((v16i8)vec0, (v16i8)(filt_h0));                    \
  vec1 = (v8i16)__msa_vshf_b((v16i8)(mask1), (v16i8)(src), (v16i8)(src));  \
  vec0 = __msa_dpadd_s_h(vec0, (v16i8)(filt_h1), (v16i8)vec1);             \
  vec2 = (v8i16)__msa_vshf_b((v16i8)(mask2), (v16i8)(src), (v16i8)(src));  \
  vec2 = __msa_dotp_s_h((v16i8)vec2, (v16i8)(filt_h2));                    \
  vec3 = (v8i16)__msa_vshf_b((v16i8)(mask3), (v16i8)(src), (v16i8)(src));  \
  vec2 = __msa_dpadd_s_h(vec2, (v16i8)(filt_h3), (v16i8)vec3);             \
  vec0 = __msa_adds_s_h(vec0, vec2);                                       \
  horiz_out = SRARI_SATURATE_SIGNED_H(vec0, FILTER_BITS, 7);               \
                                                                           \
  horiz_out;                                                               \
})

#define HORIZ_8TAP_FILT_2VECS(src0, src1, mask0, mask1, mask2, mask3,        \
                              filt_h0, filt_h1, filt_h2, filt_h3) ({         \
  v8i16 vec0, vec1, vec2, vec3, horiz_out;                                   \
                                                                             \
  vec0 = (v8i16)__msa_vshf_b((v16i8)(mask0), (v16i8)(src1), (v16i8)(src0));  \
  vec0 = __msa_dotp_s_h((v16i8)vec0, (v16i8)(filt_h0));                      \
  vec1 = (v8i16)__msa_vshf_b((v16i8)(mask1), (v16i8)(src1), (v16i8)(src0));  \
  vec0 = __msa_dpadd_s_h(vec0, (v16i8)(filt_h1), (v16i8)vec1);               \
  vec2 = (v8i16)__msa_vshf_b((v16i8)(mask2), (v16i8)(src1), (v16i8)(src0));  \
  vec2 = __msa_dotp_s_h((v16i8)vec2, (v16i8)(filt_h2));                      \
  vec3 = (v8i16)__msa_vshf_b((v16i8)(mask3), (v16i8)(src1), (v16i8)(src0));  \
  vec2 = __msa_dpadd_s_h(vec2, ((v16i8)filt_h3), (v16i8)vec3);               \
  vec0 = __msa_adds_s_h(vec0, vec2);                                         \
  horiz_out = (v8i16)SRARI_SATURATE_SIGNED_H(vec0, FILTER_BITS, 7);          \
                                                                             \
  horiz_out;                                                                 \
})

#define FILT_8TAP_DPADD_S_H(vec0, vec1, vec2, vec3,             \
                            filt0, filt1, filt2, filt3) ({      \
  v8i16 tmp0, tmp1;                                             \
                                                                \
  tmp0 = __msa_dotp_s_h((v16i8)(vec0), (v16i8)(filt0));         \
  tmp0 = __msa_dpadd_s_h(tmp0, (v16i8)(vec1), (v16i8)(filt1));  \
  tmp1 = __msa_dotp_s_h((v16i8)(vec2), (v16i8)(filt2));         \
  tmp1 = __msa_dpadd_s_h(tmp1, (v16i8)(vec3), ((v16i8)filt3));  \
  tmp0 = __msa_adds_s_h(tmp0, tmp1);                            \
                                                                \
  tmp0;                                                         \
})

#define HORIZ_8TAP_4WID_4VECS_FILT(src0, src1, src2, src3,                     \
                                   mask0, mask1, mask2, mask3,                 \
                                   filt0, filt1, filt2, filt3,                 \
                                   out0, out1) {                               \
  v8i16 vec0_m, vec1_m, vec2_m, vec3_m, vec4_m, vec5_m, vec6_m, vec7_m;        \
  v8i16 res0_m, res1_m, res2_m, res3_m;                                        \
                                                                               \
  vec0_m = (v8i16)__msa_vshf_b((v16i8)(mask0), (v16i8)(src1), (v16i8)(src0));  \
  vec1_m = (v8i16)__msa_vshf_b((v16i8)(mask0), (v16i8)(src3), (v16i8)(src2));  \
                                                                               \
  res0_m = __msa_dotp_s_h((v16i8)vec0_m, (v16i8)(filt0));                      \
  res1_m = __msa_dotp_s_h((v16i8)vec1_m, (v16i8)(filt0));                      \
                                                                               \
  vec2_m = (v8i16)__msa_vshf_b((v16i8)(mask1), (v16i8)(src1), (v16i8)(src0));  \
  vec3_m = (v8i16)__msa_vshf_b((v16i8)(mask1), (v16i8)(src3), (v16i8)(src2));  \
                                                                               \
  res0_m = __msa_dpadd_s_h(res0_m, (filt1), (v16i8)vec2_m);                    \
  res1_m = __msa_dpadd_s_h(res1_m, (filt1), (v16i8)vec3_m);                    \
                                                                               \
  vec4_m = (v8i16)__msa_vshf_b((v16i8)(mask2), (v16i8)(src1), (v16i8)(src0));  \
  vec5_m = (v8i16)__msa_vshf_b((v16i8)(mask2), (v16i8)(src3), (v16i8)(src2));  \
                                                                               \
  res2_m = __msa_dotp_s_h((v16i8)(filt2), (v16i8)vec4_m);                      \
  res3_m = __msa_dotp_s_h((v16i8)(filt2), (v16i8)vec5_m);                      \
                                                                               \
  vec6_m = (v8i16)__msa_vshf_b((v16i8)(mask3), (v16i8)(src1), (v16i8)(src0));  \
  vec7_m = (v8i16)__msa_vshf_b((v16i8)(mask3), (v16i8)(src3), (v16i8)(src2));  \
                                                                               \
  res2_m = __msa_dpadd_s_h(res2_m, (v16i8)(filt3), (v16i8)vec6_m);             \
  res3_m = __msa_dpadd_s_h(res3_m, (v16i8)(filt3), (v16i8)vec7_m);             \
                                                                               \
  out0 = __msa_adds_s_h(res0_m, res2_m);                                       \
  out1 = __msa_adds_s_h(res1_m, res3_m);                                       \
}

#define HORIZ_8TAP_8WID_4VECS_FILT(src0, src1, src2, src3,                     \
                                   mask0, mask1, mask2, mask3,                 \
                                   filt0, filt1, filt2, filt3,                 \
                                   out0, out1, out2, out3) {                   \
  v8i16 vec0_m, vec1_m, vec2_m, vec3_m;                                        \
  v8i16 vec4_m, vec5_m, vec6_m, vec7_m;                                        \
  v8i16 res0_m, res1_m, res2_m, res3_m;                                        \
  v8i16 res4_m, res5_m, res6_m, res7_m;                                        \
                                                                               \
  vec0_m = (v8i16)__msa_vshf_b((v16i8)(mask0), (v16i8)(src0), (v16i8)(src0));  \
  vec1_m = (v8i16)__msa_vshf_b((v16i8)(mask0), (v16i8)(src1), (v16i8)(src1));  \
  vec2_m = (v8i16)__msa_vshf_b((v16i8)(mask0), (v16i8)(src2), (v16i8)(src2));  \
  vec3_m = (v8i16)__msa_vshf_b((v16i8)(mask0), (v16i8)(src3), (v16i8)(src3));  \
                                                                               \
  res0_m = __msa_dotp_s_h((v16i8)vec0_m, (v16i8)(filt0));                      \
  res1_m = __msa_dotp_s_h((v16i8)vec1_m, (v16i8)(filt0));                      \
  res2_m = __msa_dotp_s_h((v16i8)vec2_m, (v16i8)(filt0));                      \
  res3_m = __msa_dotp_s_h((v16i8)vec3_m, (v16i8)(filt0));                      \
                                                                               \
  vec0_m = (v8i16)__msa_vshf_b((v16i8)(mask2), (v16i8)(src0), (v16i8)(src0));  \
  vec1_m = (v8i16)__msa_vshf_b((v16i8)(mask2), (v16i8)(src1), (v16i8)(src1));  \
  vec2_m = (v8i16)__msa_vshf_b((v16i8)(mask2), (v16i8)(src2), (v16i8)(src2));  \
  vec3_m = (v8i16)__msa_vshf_b((v16i8)(mask2), (v16i8)(src3), (v16i8)(src3));  \
                                                                               \
  res4_m = __msa_dotp_s_h((v16i8)vec0_m, (v16i8)(filt2));                      \
  res5_m = __msa_dotp_s_h((v16i8)vec1_m, (v16i8)(filt2));                      \
  res6_m = __msa_dotp_s_h((v16i8)vec2_m, (v16i8)(filt2));                      \
  res7_m = __msa_dotp_s_h((v16i8)vec3_m, (v16i8)(filt2));                      \
                                                                               \
  vec4_m = (v8i16)__msa_vshf_b((v16i8)(mask1), (v16i8)(src0), (v16i8)(src0));  \
  vec5_m = (v8i16)__msa_vshf_b((v16i8)(mask1), (v16i8)(src1), (v16i8)(src1));  \
  vec6_m = (v8i16)__msa_vshf_b((v16i8)(mask1), (v16i8)(src2), (v16i8)(src2));  \
  vec7_m = (v8i16)__msa_vshf_b((v16i8)(mask1), (v16i8)(src3), (v16i8)(src3));  \
                                                                               \
  res0_m = __msa_dpadd_s_h(res0_m, (v16i8)(filt1), (v16i8)vec4_m);             \
  res1_m = __msa_dpadd_s_h(res1_m, (v16i8)(filt1), (v16i8)vec5_m);             \
  res2_m = __msa_dpadd_s_h(res2_m, (v16i8)(filt1), (v16i8)vec6_m);             \
  res3_m = __msa_dpadd_s_h(res3_m, (v16i8)(filt1), (v16i8)vec7_m);             \
                                                                               \
  vec4_m = (v8i16)__msa_vshf_b((v16i8)(mask3), (v16i8)(src0), (v16i8)(src0));  \
  vec5_m = (v8i16)__msa_vshf_b((v16i8)(mask3), (v16i8)(src1), (v16i8)(src1));  \
  vec6_m = (v8i16)__msa_vshf_b((v16i8)(mask3), (v16i8)(src2), (v16i8)(src2));  \
  vec7_m = (v8i16)__msa_vshf_b((v16i8)(mask3), (v16i8)(src3), (v16i8)(src3));  \
                                                                               \
  res4_m = __msa_dpadd_s_h(res4_m, (v16i8)(filt3), (v16i8)vec4_m);             \
  res5_m = __msa_dpadd_s_h(res5_m, (v16i8)(filt3), (v16i8)vec5_m);             \
  res6_m = __msa_dpadd_s_h(res6_m, (v16i8)(filt3), (v16i8)vec6_m);             \
  res7_m = __msa_dpadd_s_h(res7_m, (v16i8)(filt3), (v16i8)vec7_m);             \
                                                                               \
  out0 = __msa_adds_s_h(res0_m, res4_m);                                       \
  out1 = __msa_adds_s_h(res1_m, res5_m);                                       \
  out2 = __msa_adds_s_h(res2_m, res6_m);                                       \
  out3 = __msa_adds_s_h(res3_m, res7_m);                                       \
}
#endif  /* VP9_COMMON_MIPS_MSA_VP9_CONVOLVE_MSA_H_ */
