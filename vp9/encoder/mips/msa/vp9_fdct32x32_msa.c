/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/encoder/mips/msa/vp9_fdct_msa.h"

void vp9_fdct32x32_1_msa(const int16_t *input, int16_t *out, int32_t stride) {
  out[1] = 0;

  out[0] = LD_HADD(input, stride);
  out[0] += LD_HADD(input + 8, stride);
  out[0] += LD_HADD(input + 16, stride);
  out[0] += LD_HADD(input + 24, stride);
  out[0] += LD_HADD(input + 32 * 8, stride);
  out[0] += LD_HADD(input + 32 * 8 + 8, stride);
  out[0] += LD_HADD(input + 32 * 8 + 16, stride);
  out[0] += LD_HADD(input + 32 * 8 + 24, stride);
  out[0] += LD_HADD(input + 32 * 16, stride);
  out[0] += LD_HADD(input + 32 * 16 + 8, stride);
  out[0] += LD_HADD(input + 32 * 16 + 16, stride);
  out[0] += LD_HADD(input + 32 * 16 + 24, stride);
  out[0] += LD_HADD(input + 32 * 24, stride);
  out[0] += LD_HADD(input + 32 * 24 + 8, stride);
  out[0] += LD_HADD(input + 32 * 24 + 16, stride);
  out[0] += LD_HADD(input + 32 * 24 + 24, stride);
  out[0] >>= 3;
}
