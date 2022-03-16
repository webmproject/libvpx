/*
 *  Copyright (c) 2022 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VPX_DSP_LOONGARCH_FWD_TXFM_LSX_H_
#define VPX_VPX_DSP_LOONGARCH_FWD_TXFM_LSX_H_

#include "vpx_dsp/loongarch/txfm_macros_lsx.h"
#include "vpx_dsp/txfm_common.h"

#define FDCT32_POSTPROC_2V_POS_H(vec0, vec1) \
  {                                          \
    __m128i tp0_m, tp1_m;                    \
    __m128i one = __lsx_vreplgr2vr_h(1);     \
                                             \
    tp0_m = __lsx_vslei_h(vec0, 0);          \
    tp1_m = __lsx_vslei_h(vec1, 0);          \
    tp0_m = __lsx_vxori_b(tp0_m, 255);       \
    tp1_m = __lsx_vxori_b(tp1_m, 255);       \
    vec0 = __lsx_vadd_h(vec0, one);          \
    vec1 = __lsx_vadd_h(vec1, one);          \
    tp0_m = __lsx_vand_v(one, tp0_m);        \
    tp1_m = __lsx_vand_v(one, tp1_m);        \
    vec0 = __lsx_vadd_h(vec0, tp0_m);        \
    vec1 = __lsx_vadd_h(vec1, tp1_m);        \
    vec0 = __lsx_vsrai_h(vec0, 2);           \
    vec1 = __lsx_vsrai_h(vec1, 2);           \
  }

#define FDCT_POSTPROC_2V_NEG_H(vec0, vec1) \
  {                                        \
    __m128i tp0_m, tp1_m;                  \
    __m128i one_m = __lsx_vldi(0x401);     \
                                           \
    tp0_m = __lsx_vslti_h(vec0, 0);        \
    tp1_m = __lsx_vslti_h(vec1, 0);        \
    vec0 = __lsx_vadd_h(vec0, one_m);      \
    vec1 = __lsx_vadd_h(vec1, one_m);      \
    tp0_m = __lsx_vand_v(one_m, tp0_m);    \
    tp1_m = __lsx_vand_v(one_m, tp1_m);    \
    vec0 = __lsx_vadd_h(vec0, tp0_m);      \
    vec1 = __lsx_vadd_h(vec1, tp1_m);      \
    vec0 = __lsx_vsrai_h(vec0, 2);         \
    vec1 = __lsx_vsrai_h(vec1, 2);         \
  }

#define FDCT32_POSTPROC_NEG_W(vec)         \
  {                                        \
    __m128i temp_m;                        \
    __m128i one_m = __lsx_vreplgr2vr_w(1); \
                                           \
    temp_m = __lsx_vslti_w(vec, 0);        \
    vec = __lsx_vadd_w(vec, one_m);        \
    temp_m = __lsx_vand_v(one_m, temp_m);  \
    vec = __lsx_vadd_w(vec, temp_m);       \
    vec = __lsx_vsrai_w(vec, 2);           \
  }

#define DOTP_CONST_PAIR_W(reg0_left, reg1_left, reg0_right, reg1_right,       \
                          const0, const1, out0, out1, out2, out3)             \
  {                                                                           \
    __m128i s0_m, s1_m, s2_m, s3_m, s4_m, s5_m, s6_m, s7_m;                   \
    __m128i tp0_m, tp1_m, tp2_m, tp3_m, _tmp0, _tmp1;                         \
    __m128i k0_m = __lsx_vreplgr2vr_w((int32_t)const0);                       \
                                                                              \
    s0_m = __lsx_vreplgr2vr_w((int32_t)const1);                               \
    k0_m = __lsx_vpackev_w(s0_m, k0_m);                                       \
                                                                              \
    DUP2_ARG1(__lsx_vneg_w, reg1_left, reg1_right, _tmp0, _tmp1);             \
    s1_m = __lsx_vilvl_w(_tmp0, reg0_left);                                   \
    s0_m = __lsx_vilvh_w(_tmp0, reg0_left);                                   \
    s3_m = __lsx_vilvl_w(reg0_left, reg1_left);                               \
    s2_m = __lsx_vilvh_w(reg0_left, reg1_left);                               \
    s5_m = __lsx_vilvl_w(_tmp1, reg0_right);                                  \
    s4_m = __lsx_vilvh_w(_tmp1, reg0_right);                                  \
    s7_m = __lsx_vilvl_w(reg0_right, reg1_right);                             \
    s6_m = __lsx_vilvh_w(reg0_right, reg1_right);                             \
    DUP2_ARG2(__lsx_vdp2_d_w, s0_m, k0_m, s1_m, k0_m, tp0_m, tp1_m);          \
    DUP2_ARG2(__lsx_vdp2_d_w, s4_m, k0_m, s5_m, k0_m, tp2_m, tp3_m);          \
    DUP2_ARG3(__lsx_vssrarni_w_d, tp0_m, tp1_m, DCT_CONST_BITS, tp2_m, tp3_m, \
              DCT_CONST_BITS, out0, out1);                                    \
    DUP2_ARG2(__lsx_vdp2_d_w, s2_m, k0_m, s3_m, k0_m, tp0_m, tp1_m);          \
    DUP2_ARG2(__lsx_vdp2_d_w, s6_m, k0_m, s7_m, k0_m, tp2_m, tp3_m);          \
    DUP2_ARG3(__lsx_vssrarni_w_d, tp0_m, tp1_m, DCT_CONST_BITS, tp2_m, tp3_m, \
              DCT_CONST_BITS, out2, out3);                                    \
  }

#endif  // VPX_VPX_DSP_LOONGARCH_FWD_TXFM_LSX_H_
