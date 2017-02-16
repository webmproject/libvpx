/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#ifndef VPX_DSP_X86_BITDEPTH_CONVERSION_AVX2_H_
#define VPX_DSP_X86_BITDEPTH_CONVERSION_AVX2_H_

#include <immintrin.h>

#include "./vpx_config.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/vpx_dsp_common.h"

// Load 16 16 bit values. If the source is 32 bits then pack down with
// saturation.
static INLINE __m256i load_tran_low(const tran_low_t *a) {
#if CONFIG_VP9_HIGHBITDEPTH
  const __m256i a_low = _mm256_loadu_si256((const __m256i *)a);
  return _mm256_packs_epi32(a_low, *(const __m256i *)(a + 8));
#else
  return _mm256_loadu_si256((const __m256i *)a);
#endif
}

#endif  // VPX_DSP_X86_BITDEPTH_CONVERSION_AVX2_H_
