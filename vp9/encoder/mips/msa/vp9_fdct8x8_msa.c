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

#include "./vp9_rtcd.h"
#include "vp9/encoder/mips/msa/vp9_fdct_msa.h"

void vp9_fdct8x8_msa(const int16_t *input, int16_t *output,
                     int32_t src_stride) {
  v8i16 in0, in1, in2, in3, in4, in5, in6, in7;

  LD_SH8(input, src_stride, in0, in1, in2, in3, in4, in5, in6, in7);
  SLLI_4V(in0, in1, in2, in3, 2);
  SLLI_4V(in4, in5, in6, in7, 2);
  VP9_FDCT8(in0, in1, in2, in3, in4, in5, in6, in7,
            in0, in1, in2, in3, in4, in5, in6, in7);
  TRANSPOSE8x8_SH_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                     in0, in1, in2, in3, in4, in5, in6, in7);
  VP9_FDCT8(in0, in1, in2, in3, in4, in5, in6, in7,
            in0, in1, in2, in3, in4, in5, in6, in7);
  TRANSPOSE8x8_SH_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                     in0, in1, in2, in3, in4, in5, in6, in7);
  VP9_SRLI_AVE_S_4V_H(in0, in1, in2, in3, in4, in5, in6, in7);
  ST_SH8(in0, in1, in2, in3, in4, in5, in6, in7, output, 8);
}

void vp9_fdct8x8_1_msa(const int16_t *input, int16_t *out, int32_t stride) {
  out[0] = VP9_LD_HADD(input, stride);
  out[1] = 0;
}

void vp9_fht8x8_msa(const int16_t *input, int16_t *output, int32_t stride,
                    int32_t tx_type) {
  v8i16 in0, in1, in2, in3, in4, in5, in6, in7;

  LD_SH8(input, stride, in0, in1, in2, in3, in4, in5, in6, in7);
  SLLI_4V(in0, in1, in2, in3, 2);
  SLLI_4V(in4, in5, in6, in7, 2);

  switch (tx_type) {
    case DCT_DCT:
      VP9_FDCT8(in0, in1, in2, in3, in4, in5, in6, in7,
                in0, in1, in2, in3, in4, in5, in6, in7);
      TRANSPOSE8x8_SH_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                         in0, in1, in2, in3, in4, in5, in6, in7);
      VP9_FDCT8(in0, in1, in2, in3, in4, in5, in6, in7,
                in0, in1, in2, in3, in4, in5, in6, in7);
      break;
    case ADST_DCT:
      VP9_ADST8(in0, in1, in2, in3, in4, in5, in6, in7,
                in0, in1, in2, in3, in4, in5, in6, in7);
      TRANSPOSE8x8_SH_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                         in0, in1, in2, in3, in4, in5, in6, in7);
      VP9_FDCT8(in0, in1, in2, in3, in4, in5, in6, in7,
                in0, in1, in2, in3, in4, in5, in6, in7);
      break;
    case DCT_ADST:
      VP9_FDCT8(in0, in1, in2, in3, in4, in5, in6, in7,
                in0, in1, in2, in3, in4, in5, in6, in7);
      TRANSPOSE8x8_SH_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                         in0, in1, in2, in3, in4, in5, in6, in7);
      VP9_ADST8(in0, in1, in2, in3, in4, in5, in6, in7,
                in0, in1, in2, in3, in4, in5, in6, in7);
      break;
    case ADST_ADST:
      VP9_ADST8(in0, in1, in2, in3, in4, in5, in6, in7,
                in0, in1, in2, in3, in4, in5, in6, in7);
      TRANSPOSE8x8_SH_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                         in0, in1, in2, in3, in4, in5, in6, in7);
      VP9_ADST8(in0, in1, in2, in3, in4, in5, in6, in7,
                in0, in1, in2, in3, in4, in5, in6, in7);
      break;
    default:
      assert(0);
      break;
  }

  TRANSPOSE8x8_SH_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                     in0, in1, in2, in3, in4, in5, in6, in7);
  VP9_SRLI_AVE_S_4V_H(in0, in1, in2, in3, in4, in5, in6, in7);
  ST_SH8(in0, in1, in2, in3, in4, in5, in6, in7, output, 8);
}
