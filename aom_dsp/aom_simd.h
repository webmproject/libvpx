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

#ifndef _AOM_SIMD_H
#define _AOM_SIMD_H

#ifndef SIMD_INLINE
#ifdef __GNUC__
#define SIMD_INLINE static inline __attribute__((always_inline))
#elif __STDC_VERSION__ >= 199901L
#define SIMD_INLINE static inline
#elif defined(_MSC_VER)
#define SIMD_INLINE static __inline
#else
#define SIMD_INLINE static
#endif
#endif

#include <stdint.h>

#if defined(_WIN32)
#include <intrin.h>
#endif

#include "./aom_config.h"

#if HAVE_NEON
#include "simd/v128_intrinsics_arm.h"
#elif HAVE_SSE2
#include "simd/v128_intrinsics_x86.h"
#else
#include "simd/v128_intrinsics.h"
#endif

#endif /* _AOM_SIMD_H */
