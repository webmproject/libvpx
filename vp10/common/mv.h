/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_COMMON_MV_H_
#define VP10_COMMON_MV_H_

#include "vp10/common/common.h"
#include "vpx_dsp/vpx_filter.h"
#if CONFIG_GLOBAL_MOTION
#include "vp10/common/warped_motion.h"
#endif  // CONFIG_GLOBAL_MOTION

#ifdef __cplusplus
extern "C" {
#endif

typedef struct mv {
  int16_t row;
  int16_t col;
} MV;

typedef union int_mv {
  uint32_t as_int;
  MV as_mv;
} int_mv; /* facilitates faster equality tests and copies */

typedef struct mv32 {
  int32_t row;
  int32_t col;
} MV32;

#if CONFIG_GLOBAL_MOTION
// ALPHA here refers to parameters a and b in rotzoom model:
// | a   b|
// |-b   a|
//
// Anything ending in PREC_BITS is the number of bits of precision
// to maintain when converting from double to integer.
//
// The ABS parameters are used to create an upper and lower bound
// for each parameter. In other words, after a parameter is integerized
// it is clamped between -(1 << ABS_XXX_BITS) and (1 << ABS_XXX_BITS)
//
// XXX_PREC_DIFF and XXX_DECODE_FACTOR are computed once here to
// prevent repetetive computation on the decoder side. These are
// to allow the global motion parameters to be encoded in a lower
// precision than the warped model precision. This means that they
// need to be changed to warped precision when they are decoded.
//
// XX_MIN, XX_MAX are also computed to avoid repeated computation

#define GM_TRANS_PREC_BITS        3
#define GM_TRANS_PREC_DIFF        (WARPEDMODEL_PREC_BITS - GM_TRANS_PREC_BITS)
#define GM_TRANS_DECODE_FACTOR    (1 << GM_TRANS_PREC_DIFF)

#define GM_ALPHA_PREC_BITS        5
#define GM_ALPHA_PREC_DIFF        (WARPEDMODEL_PREC_BITS - GM_ALPHA_PREC_BITS)
#define GM_ALPHA_DECODE_FACTOR    (1 << GM_ALPHA_PREC_DIFF)

#define GM_ABS_ALPHA_BITS         5
#define GM_ABS_TRANS_BITS         5

#define GM_TRANS_MAX              (1 << GM_ABS_TRANS_BITS)
#define GM_ALPHA_MAX              (1 << GM_ABS_ALPHA_BITS)
#define GM_TRANS_MIN              -GM_TRANS_MAX
#define GM_ALPHA_MIN              -GM_ALPHA_MAX

typedef enum {
  GLOBAL_ZERO = 0,
  GLOBAL_TRANSLATION = 1,
  GLOBAL_ROTZOOM = 2,
  GLOBAL_MOTION_TYPES
} GLOBAL_MOTION_TYPE;

typedef struct {
  GLOBAL_MOTION_TYPE gmtype;
  WarpedMotionParams motion_params;
} Global_Motion_Params;

static INLINE GLOBAL_MOTION_TYPE get_gmtype(const Global_Motion_Params *gm) {
  if (gm->motion_params.wmmat[2] == 0 && gm->motion_params.wmmat[3] == 0) {
    return ((gm->motion_params.wmmat[0] | gm->motion_params.wmmat[1]) ?
            GLOBAL_TRANSLATION : GLOBAL_ZERO);
  } else {
    return GLOBAL_ROTZOOM;
  }
}
#endif  // CONFIG_GLOBAL_MOTION

#if CONFIG_REF_MV
typedef struct candidate_mv {
  int_mv this_mv;
  int_mv comp_mv;
  int_mv pred_mv;
  int weight;
} CANDIDATE_MV;
#endif

static INLINE int is_zero_mv(const MV *mv) {
  return *((const uint32_t *)mv) == 0;
}

static INLINE int is_equal_mv(const MV *a, const MV *b) {
  return  *((const uint32_t *)a) == *((const uint32_t *)b);
}

static INLINE void clamp_mv(MV *mv, int min_col, int max_col,
                            int min_row, int max_row) {
  mv->col = clamp(mv->col, min_col, max_col);
  mv->row = clamp(mv->row, min_row, max_row);
}

static INLINE int mv_has_subpel(const MV *mv) {
  return (mv->row & SUBPEL_MASK) || (mv->col & SUBPEL_MASK);
}
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_COMMON_MV_H_
