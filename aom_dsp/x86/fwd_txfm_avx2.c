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

#include "./aom_config.h"

#define FDCT32x32_2D_AVX2 aom_fdct32x32_rd_avx2
#define FDCT32x32_HIGH_PRECISION 0
#include "aom_dsp/x86/fwd_dct32x32_impl_avx2.h"
#undef FDCT32x32_2D_AVX2
#undef FDCT32x32_HIGH_PRECISION

// TODO(luoyi): The following macro hides an error. The second parameter type of
// function,
//   void FDCT32x32_2D_AVX2(const int16_t *, int16_t*, int);
// is different from the one in,
//   void aom_fdct32x32_avx2(const int16_t *, tran_low_t*, int);
// In CONFIG_AOM_HIGHBITDEPTH=1 build, the second parameter type should be
// int32_t.
// This function should be removed after av1_fht32x32 scaling/rounding fix.
#define FDCT32x32_2D_AVX2 aom_fdct32x32_avx2
#define FDCT32x32_HIGH_PRECISION 1
#include "aom_dsp/x86/fwd_dct32x32_impl_avx2.h"  // NOLINT
#undef FDCT32x32_2D_AVX2
#undef FDCT32x32_HIGH_PRECISION
