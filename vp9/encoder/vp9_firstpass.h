/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_FIRSTPASS_H_
#define VP9_ENCODER_VP9_FIRSTPASS_H_

#include "vp9/encoder/vp9_lookahead.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double frame;
  double intra_error;
  double coded_error;
  double sr_coded_error;
  double pcnt_inter;
  double pcnt_motion;
  double pcnt_second_ref;
  double pcnt_neutral;
  double MVr;
  double mvr_abs;
  double MVc;
  double mvc_abs;
  double MVrv;
  double MVcv;
  double mv_in_out_count;
  double new_mv_count;
  double duration;
  double count;
  int64_t spatial_layer_id;
} FIRSTPASS_STATS;

typedef struct {
  unsigned int section_intra_rating;
  unsigned int next_iiratio;
  FIRSTPASS_STATS total_stats;
  FIRSTPASS_STATS this_frame_stats;
  const FIRSTPASS_STATS *stats_in;
  const FIRSTPASS_STATS *stats_in_start;
  const FIRSTPASS_STATS *stats_in_end;
  FIRSTPASS_STATS total_left_stats;
  int first_pass_done;
  int64_t bits_left;
  double modified_error_min;
  double modified_error_max;
  double modified_error_left;
  double kf_intra_err_min;
  double gf_intra_err_min;

  // Projected total bits available for a key frame group of frames
  int64_t kf_group_bits;

  // Error score of frames still to be coded in kf group
  int64_t kf_group_error_left;
  int sr_update_lag;

  int kf_zeromotion_pct;
  int gf_zeromotion_pct;

  int active_worst_quality;

  int gf_group_index;
  int gf_group_bit_allocation[MAX_LAG_BUFFERS * 2];
} TWO_PASS;

struct VP9_COMP;

void vp9_init_first_pass(struct VP9_COMP *cpi);
void vp9_rc_get_first_pass_params(struct VP9_COMP *cpi);
void vp9_first_pass(struct VP9_COMP *cpi);
void vp9_end_first_pass(struct VP9_COMP *cpi);

void vp9_init_second_pass(struct VP9_COMP *cpi);
void vp9_rc_get_second_pass_params(struct VP9_COMP *cpi);

// Post encode update of the rate control parameters for 2-pass
void vp9_twopass_postencode_update(struct VP9_COMP *cpi);
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_FIRSTPASS_H_
