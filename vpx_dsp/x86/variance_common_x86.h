/*
 *  Copyright (c) 2026 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VPX_DSP_X86_VARIANCE_COMMON_X86_H_
#define VPX_VPX_DSP_X86_VARIANCE_COMMON_X86_H_

#include "./vpx_config.h"
#include "vpx/vpx_integer.h"

// Computes sse - ((sum^2) >> shift) for 8-bit blocks.
// When w * h <= 256 (shift <= 8), |se| <= 256 * 255 = 65280 so (int64_t)se * se
// fits within uint32_t and can be cast before shifting. Larger blocks can reach
// 40-bit squared sums and must shift in 64-bit arithmetic before truncating.
static INLINE uint32_t variance_from_se_8b(uint32_t sse, int se, int shift) {
  if (shift <= 8) {
    return sse - ((uint32_t)((int64_t)se * se) >> shift);
  }
  return sse - (uint32_t)(((int64_t)se * se) >> shift);
}

#endif  // VPX_VPX_DSP_X86_VARIANCE_COMMON_X86_H_
