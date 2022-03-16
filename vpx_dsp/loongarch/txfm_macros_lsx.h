/*
 *  Copyright (c) 2022 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VPX_DSP_LOONGARCH_TXFM_MACROS_LSX_H_
#define VPX_VPX_DSP_LOONGARCH_TXFM_MACROS_LSX_H_

#include "vpx_util/loongson_intrinsics.h"

#define DOTP_CONST_PAIR(reg0, reg1, cnst0, cnst1, out0, out1) \
  {                                                           \
    __m128i s0_m, s1_m, s2_m, s3_m, s4_m, s5_m;               \
    __m128i k0_m, k1_m, k2_m, k3_m;                           \
    __m128i zero = __lsx_vldi(0);                             \
                                                              \
    k0_m = __lsx_vreplgr2vr_h(cnst0);                         \
    k1_m = __lsx_vreplgr2vr_h(cnst1);                         \
    k2_m = __lsx_vpackev_h(k1_m, k0_m);                       \
    k0_m = __lsx_vpackev_h(zero, k0_m);                       \
    k1_m = __lsx_vpackev_h(k1_m, zero);                       \
                                                              \
    s5_m = __lsx_vilvl_h(reg1, reg0);                         \
    s4_m = __lsx_vilvh_h(reg1, reg0);                         \
    s3_m = __lsx_vilvl_h(reg0, reg1);                         \
    s2_m = __lsx_vilvh_h(reg0, reg1);                         \
                                                              \
    s1_m = __lsx_vdp2_w_h(s5_m, k0_m);                        \
    s0_m = __lsx_vdp2_w_h(s4_m, k0_m);                        \
    k3_m = __lsx_vdp2_w_h(s5_m, k1_m);                        \
    s1_m = __lsx_vsub_w(s1_m, k3_m);                          \
    k3_m = __lsx_vdp2_w_h(s4_m, k1_m);                        \
    s0_m = __lsx_vsub_w(s0_m, k3_m);                          \
                                                              \
    out0 = __lsx_vssrarni_h_w(s0_m, s1_m, DCT_CONST_BITS);    \
                                                              \
    s1_m = __lsx_vdp2_w_h(s3_m, k2_m);                        \
    s0_m = __lsx_vdp2_w_h(s2_m, k2_m);                        \
    out1 = __lsx_vssrarni_h_w(s0_m, s1_m, DCT_CONST_BITS);    \
  }

#endif  // VPX_VPX_DSP_LOONGARCH_TXFM_MACROS_LSX_H_
