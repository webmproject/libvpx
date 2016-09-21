/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_DSP_X86_TXFM_COMMON_AVX2_H
#define AOM_DSP_X86_TXFM_COMMON_AVX2_H

#include <immintrin.h>

#define pair256_set_epi16(a, b)                                            \
  _mm256_set_epi16((int16_t)(b), (int16_t)(a), (int16_t)(b), (int16_t)(a), \
                   (int16_t)(b), (int16_t)(a), (int16_t)(b), (int16_t)(a), \
                   (int16_t)(b), (int16_t)(a), (int16_t)(b), (int16_t)(a), \
                   (int16_t)(b), (int16_t)(a), (int16_t)(b), (int16_t)(a))

#define pair256_set_epi32(a, b)                                                \
  _mm256_set_epi32((int)(b), (int)(a), (int)(b), (int)(a), (int)(b), (int)(a), \
                   (int)(b), (int)(a))

#endif  // AOM_DSP_X86_TXFM_COMMON_AVX2_H
