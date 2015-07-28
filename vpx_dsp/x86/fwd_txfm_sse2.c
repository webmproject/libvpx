/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_config.h"

#define DCT_HIGH_BIT_DEPTH 0
#define FDCT4x4_2D vp9_fdct4x4_sse2
#define FDCT8x8_2D vp9_fdct8x8_sse2
#define FDCT16x16_2D vp9_fdct16x16_sse2
#include "vpx_dsp/x86/fwd_txfm_impl_sse2.h"
#undef  FDCT4x4_2D
#undef  FDCT8x8_2D
#undef  FDCT16x16_2D

#define FDCT32x32_2D vp9_fdct32x32_rd_sse2
#define FDCT32x32_HIGH_PRECISION 0
#include "vpx_dsp/x86/fwd_dct32x32_impl_sse2.h"
#undef  FDCT32x32_2D
#undef  FDCT32x32_HIGH_PRECISION

#define FDCT32x32_2D vp9_fdct32x32_sse2
#define FDCT32x32_HIGH_PRECISION 1
#include "vpx_dsp/x86/fwd_dct32x32_impl_sse2.h"  // NOLINT
#undef  FDCT32x32_2D
#undef  FDCT32x32_HIGH_PRECISION
#undef  DCT_HIGH_BIT_DEPTH

#if CONFIG_VP9_HIGHBITDEPTH
#define DCT_HIGH_BIT_DEPTH 1
#define FDCT4x4_2D vp9_highbd_fdct4x4_sse2
#define FDCT8x8_2D vp9_highbd_fdct8x8_sse2
#define FDCT16x16_2D vp9_highbd_fdct16x16_sse2
#include "vpx_dsp/x86/fwd_txfm_impl_sse2.h" // NOLINT
#undef  FDCT4x4_2D
#undef  FDCT8x8_2D
#undef  FDCT16x16_2D

#define FDCT32x32_2D vp9_highbd_fdct32x32_rd_sse2
#define FDCT32x32_HIGH_PRECISION 0
#include "vpx_dsp/x86/fwd_dct32x32_impl_sse2.h" // NOLINT
#undef  FDCT32x32_2D
#undef  FDCT32x32_HIGH_PRECISION

#define FDCT32x32_2D vp9_highbd_fdct32x32_sse2
#define FDCT32x32_HIGH_PRECISION 1
#include "vpx_dsp/x86/fwd_dct32x32_impl_sse2.h" // NOLINT
#undef  FDCT32x32_2D
#undef  FDCT32x32_HIGH_PRECISION
#undef  DCT_HIGH_BIT_DEPTH
#endif  // CONFIG_VP9_HIGHBITDEPTH
