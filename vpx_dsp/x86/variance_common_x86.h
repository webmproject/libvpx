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
#include "vpx_ports/mem.h"

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

#if CONFIG_VP9_HIGHBITDEPTH
#define VPX_DISPATCH_SPV_HBD(core_fn, src, ss, xo, yo, ref, rs, sec, sec_s,    \
                             do_sec, w, h, sse, sum)                           \
  do {                                                                         \
    if ((xo) == 0) {                                                           \
      if ((yo) == 0) {                                                         \
        core_fn(src, ss, xo, yo, ref, rs, sec, sec_s, do_sec, w, h, 0, 0, sse, \
                sum);                                                          \
      } else if ((yo) == 4) {                                                  \
        core_fn(src, ss, xo, yo, ref, rs, sec, sec_s, do_sec, w, h, 0, 1, sse, \
                sum);                                                          \
      } else {                                                                 \
        core_fn(src, ss, xo, yo, ref, rs, sec, sec_s, do_sec, w, h, 0, 2, sse, \
                sum);                                                          \
      }                                                                        \
    } else if ((xo) == 4) {                                                    \
      if ((yo) == 0) {                                                         \
        core_fn(src, ss, xo, yo, ref, rs, sec, sec_s, do_sec, w, h, 1, 0, sse, \
                sum);                                                          \
      } else if ((yo) == 4) {                                                  \
        core_fn(src, ss, xo, yo, ref, rs, sec, sec_s, do_sec, w, h, 1, 1, sse, \
                sum);                                                          \
      } else {                                                                 \
        core_fn(src, ss, xo, yo, ref, rs, sec, sec_s, do_sec, w, h, 1, 2, sse, \
                sum);                                                          \
      }                                                                        \
    } else {                                                                   \
      if ((yo) == 0) {                                                         \
        core_fn(src, ss, xo, yo, ref, rs, sec, sec_s, do_sec, w, h, 2, 0, sse, \
                sum);                                                          \
      } else if ((yo) == 4) {                                                  \
        core_fn(src, ss, xo, yo, ref, rs, sec, sec_s, do_sec, w, h, 2, 1, sse, \
                sum);                                                          \
      } else {                                                                 \
        core_fn(src, ss, xo, yo, ref, rs, sec, sec_s, do_sec, w, h, 2, 2, sse, \
                sum);                                                          \
      }                                                                        \
    }                                                                          \
  } while (0)

static INLINE uint32_t finalize_highbd_var_8(uint64_t sse_long,
                                             int64_t sum_long, int shift,
                                             uint32_t *const sse) {
  *sse = (uint32_t)sse_long;
  return variance_from_se_8b(*sse, (int)sum_long, shift);
}

static INLINE uint32_t finalize_highbd_var_10(uint64_t sse_long,
                                              int64_t sum_long, int shift,
                                              uint32_t *const sse) {
  const int sum = (int)ROUND_POWER_OF_TWO(sum_long, 2);
  int64_t var;
  *sse = (uint32_t)ROUND_POWER_OF_TWO(sse_long, 4);
  var = (int64_t)(*sse) - (((int64_t)sum * sum) >> shift);
  return (var >= 0) ? (uint32_t)var : 0;
}

static INLINE uint32_t finalize_highbd_var_12(uint64_t sse_long,
                                              int64_t sum_long, int shift,
                                              uint32_t *const sse) {
  const int sum = (int)ROUND_POWER_OF_TWO(sum_long, 4);
  int64_t var;
  *sse = (uint32_t)ROUND_POWER_OF_TWO(sse_long, 8);
  var = (int64_t)(*sse) - (((int64_t)sum * sum) >> shift);
  return (var >= 0) ? (uint32_t)var : 0;
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

#endif  // VPX_VPX_DSP_X86_VARIANCE_COMMON_X86_H_
