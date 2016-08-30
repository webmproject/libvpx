/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_COMMON_SCALE_H_
#define AV1_COMMON_SCALE_H_

#include "av1/common/mv.h"
#include "aom_dsp/aom_convolve.h"

#ifdef __cplusplus
extern "C" {
#endif

#define REF_SCALE_SHIFT 14
#define REF_NO_SCALE (1 << REF_SCALE_SHIFT)
#define REF_INVALID_SCALE -1

struct scale_factors {
  int x_scale_fp;  // horizontal fixed point scale factor
  int y_scale_fp;  // vertical fixed point scale factor
  int x_step_q4;
  int y_step_q4;

  int (*scale_value_x)(int val, const struct scale_factors *sf);
  int (*scale_value_y)(int val, const struct scale_factors *sf);

  convolve_fn_t predict[2][2][2];  // horiz, vert, avg
#if CONFIG_AOM_HIGHBITDEPTH
  highbd_convolve_fn_t highbd_predict[2][2][2];  // horiz, vert, avg
#endif                                           // CONFIG_AOM_HIGHBITDEPTH

// Functions for non-interpolating filters (those that filter zero offsets)
#if CONFIG_EXT_INTERP && SUPPORT_NONINTERPOLATING_FILTERS
  convolve_fn_t predict_ni[2][2][2];  // horiz, vert, avg
#if CONFIG_AOM_HIGHBITDEPTH
  highbd_convolve_fn_t highbd_predict_ni[2][2][2];  // horiz, vert, avg
#endif                                              // CONFIG_AOM_HIGHBITDEPTH
#endif  // CONFIG_EXT_INTERP && SUPPORT_NONINTERPOLATING_FILTERS
};

MV32 av1_scale_mv(const MV *mv, int x, int y, const struct scale_factors *sf);

#if CONFIG_AOM_HIGHBITDEPTH
void av1_setup_scale_factors_for_frame(struct scale_factors *sf, int other_w,
                                       int other_h, int this_w, int this_h,
                                       int use_high);
#else
void av1_setup_scale_factors_for_frame(struct scale_factors *sf, int other_w,
                                       int other_h, int this_w, int this_h);
#endif  // CONFIG_AOM_HIGHBITDEPTH

static INLINE int av1_is_valid_scale(const struct scale_factors *sf) {
  return sf->x_scale_fp != REF_INVALID_SCALE &&
         sf->y_scale_fp != REF_INVALID_SCALE;
}

static INLINE int av1_is_scaled(const struct scale_factors *sf) {
  return av1_is_valid_scale(sf) &&
         (sf->x_scale_fp != REF_NO_SCALE || sf->y_scale_fp != REF_NO_SCALE);
}

static INLINE int valid_ref_frame_size(int ref_width, int ref_height,
                                       int this_width, int this_height) {
  return 2 * this_width >= ref_width && 2 * this_height >= ref_height &&
         this_width <= 16 * ref_width && this_height <= 16 * ref_height;
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_SCALE_H_
