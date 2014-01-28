/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP9_ENCODER_VP9_RATECTRL_H_
#define VP9_ENCODER_VP9_RATECTRL_H_

#ifdef __cplusplus
extern "C" {
#endif

#define FRAME_OVERHEAD_BITS 200

typedef struct {
  // Rate targetting variables
  int this_frame_target;
  int projected_frame_size;
  int sb64_target_rate;
  int last_q[3];                   // Separate values for Intra/Inter/ARF-GF
  int last_boosted_qindex;         // Last boosted GF/KF/ARF q

  int gfu_boost;
  int last_boost;
  int kf_boost;

  double rate_correction_factor;
  double key_frame_rate_correction_factor;
  double gf_rate_correction_factor;

  unsigned int frames_since_golden;
  unsigned int frames_till_gf_update_due;  // Count down till next GF
  unsigned int max_gf_interval;
  unsigned int baseline_gf_interval;
  unsigned int frames_to_key;
  unsigned int frames_since_key;
  unsigned int this_key_frame_forced;
  unsigned int next_key_frame_forced;
  unsigned int source_alt_ref_pending;
  unsigned int source_alt_ref_active;
  unsigned int is_src_frame_alt_ref;

  int per_frame_bandwidth;        // Current section per frame bandwidth target
  int av_per_frame_bandwidth;     // Average frame size target for clip
  int min_frame_bandwidth;        // Minimum allocation used for any frame
  int max_frame_bandwidth;        // Maximum burst rate allowed for a frame.

  int ni_av_qi;
  int ni_tot_qi;
  int ni_frames;
  int avg_frame_qindex[3];  // 0 - KEY, 1 - INTER, 2 - ARF/GF
  double tot_q;
  double avg_q;

  int buffer_level;
  int bits_off_target;

  int decimation_factor;
  int decimation_count;

  int rolling_target_bits;
  int rolling_actual_bits;

  int long_rolling_target_bits;
  int long_rolling_actual_bits;

  int64_t total_actual_bits;
  int total_target_vs_actual;        // debug stats

  int worst_quality;
  int active_worst_quality;
  int best_quality;
  // int active_best_quality;
} RATE_CONTROL;

struct VP9_COMP;

void vp9_save_coding_context(struct VP9_COMP *cpi);
void vp9_restore_coding_context(struct VP9_COMP *cpi);

void vp9_setup_key_frame(struct VP9_COMP *cpi);
void vp9_setup_inter_frame(struct VP9_COMP *cpi);

double vp9_convert_qindex_to_q(int qindex);

// Updates rate correction factors
void vp9_rc_update_rate_correction_factors(struct VP9_COMP *cpi, int damp_var);

// initialize luts for minq
void vp9_rc_init_minq_luts(void);

// return of 0 means drop frame
// Changes only rc.this_frame_target and rc.sb64_rate_target
int vp9_rc_pick_frame_size_target(struct VP9_COMP *cpi);

void vp9_rc_compute_frame_size_bounds(const struct VP9_COMP *cpi,
                                      int this_frame_target,
                                      int *frame_under_shoot_limit,
                                      int *frame_over_shoot_limit);

// Picks q and q bounds given the target for bits
int vp9_rc_pick_q_and_adjust_q_bounds(const struct VP9_COMP *cpi,
                                      int *bottom_index,
                                      int *top_index);

// Estimates q to achieve a target bits per frame
int vp9_rc_regulate_q(const struct VP9_COMP *cpi, int target_bits_per_frame,
                      int active_best_quality, int active_worst_quality);

// Post encode update of the rate control parameters based
// on bytes used
void vp9_rc_postencode_update(struct VP9_COMP *cpi,
                              uint64_t bytes_used);
// for dropped frames
void vp9_rc_postencode_update_drop_frame(struct VP9_COMP *cpi);

// estimates bits per mb for a given qindex and correction factor
int vp9_rc_bits_per_mb(FRAME_TYPE frame_type, int qindex,
                       double correction_factor);

// Post encode update of the rate control parameters for 2-pass
void vp9_twopass_postencode_update(struct VP9_COMP *cpi,
                                   uint64_t bytes_used);

// Decide if we should drop this frame: For 1-pass CBR.
int vp9_drop_frame(struct VP9_COMP *cpi);

// Update the buffer level.
void vp9_update_buffer_level(struct VP9_COMP *cpi, int encoded_frame_size);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_RATECTRL_H_
