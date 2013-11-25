/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#include "vp9/common/vp9_alloccommon.h"
#include "vp9/common/vp9_common.h"
#include "vp9/encoder/vp9_ratectrl.h"
#include "vp9/common/vp9_entropymode.h"
#include "vpx_mem/vpx_mem.h"
#include "vp9/common/vp9_systemdependent.h"
#include "vp9/encoder/vp9_encodemv.h"
#include "vp9/common/vp9_quant_common.h"
#include "vp9/common/vp9_seg_common.h"

#define MIN_BPB_FACTOR 0.005
#define MAX_BPB_FACTOR 50

// Bits Per MB at different Q (Multiplied by 512)
#define BPER_MB_NORMBITS    9

static const unsigned int prior_key_frame_weight[KEY_FRAME_CONTEXT] =
    { 1, 2, 3, 4, 5 };

// Tables relating active max Q to active min Q
static int kf_low_motion_minq[QINDEX_RANGE];
static int kf_high_motion_minq[QINDEX_RANGE];
static int gf_low_motion_minq[QINDEX_RANGE];
static int gf_high_motion_minq[QINDEX_RANGE];
static int inter_minq[QINDEX_RANGE];
static int afq_low_motion_minq[QINDEX_RANGE];
static int afq_high_motion_minq[QINDEX_RANGE];

// Functions to compute the active minq lookup table entries based on a
// formulaic approach to facilitate easier adjustment of the Q tables.
// The formulae were derived from computing a 3rd order polynomial best
// fit to the original data (after plotting real maxq vs minq (not q index))
static int calculate_minq_index(double maxq,
                                double x3, double x2, double x1, double c) {
  int i;
  const double minqtarget = MIN(((x3 * maxq + x2) * maxq + x1) * maxq + c,
                                maxq);

  // Special case handling to deal with the step from q2.0
  // down to lossless mode represented by q 1.0.
  if (minqtarget <= 2.0)
    return 0;

  for (i = 0; i < QINDEX_RANGE; i++) {
    if (minqtarget <= vp9_convert_qindex_to_q(i))
      return i;
  }

  return QINDEX_RANGE - 1;
}

void vp9_init_minq_luts(void) {
  int i;

  for (i = 0; i < QINDEX_RANGE; i++) {
    const double maxq = vp9_convert_qindex_to_q(i);


    kf_low_motion_minq[i] = calculate_minq_index(maxq,
                                                 0.000001,
                                                 -0.0004,
                                                 0.15,
                                                 0.0);
    kf_high_motion_minq[i] = calculate_minq_index(maxq,
                                                  0.000002,
                                                  -0.0012,
                                                  0.5,
                                                  0.0);

    gf_low_motion_minq[i] = calculate_minq_index(maxq,
                                                 0.0000015,
                                                 -0.0009,
                                                 0.32,
                                                 0.0);
    gf_high_motion_minq[i] = calculate_minq_index(maxq,
                                                  0.0000021,
                                                  -0.00125,
                                                  0.50,
                                                  0.0);
    inter_minq[i] = calculate_minq_index(maxq,
                                         0.00000271,
                                         -0.00113,
                                         0.75,
                                         0.0);
    afq_low_motion_minq[i] = calculate_minq_index(maxq,
                                                  0.0000015,
                                                  -0.0009,
                                                  0.33,
                                                  0.0);
    afq_high_motion_minq[i] = calculate_minq_index(maxq,
                                                   0.0000021,
                                                   -0.00125,
                                                   0.55,
                                                   0.0);
  }
}

// These functions use formulaic calculations to make playing with the
// quantizer tables easier. If necessary they can be replaced by lookup
// tables if and when things settle down in the experimental bitstream
double vp9_convert_qindex_to_q(int qindex) {
  // Convert the index to a real Q value (scaled down to match old Q values)
  return vp9_ac_quant(qindex, 0) / 4.0;
}

int vp9_gfboost_qadjust(int qindex) {
  const double q = vp9_convert_qindex_to_q(qindex);
  return (int)((0.00000828 * q * q * q) +
               (-0.0055 * q * q) +
               (1.32 * q) + 79.3);
}

static int kfboost_qadjust(int qindex) {
  const double q = vp9_convert_qindex_to_q(qindex);
  return (int)((0.00000973 * q * q * q) +
               (-0.00613 * q * q) +
               (1.316 * q) + 121.2);
}

int vp9_bits_per_mb(FRAME_TYPE frame_type, int qindex,
                    double correction_factor) {
  const double q = vp9_convert_qindex_to_q(qindex);
  int enumerator = frame_type == KEY_FRAME ? 3300000 : 2250000;

  // q based adjustment to baseline enumerator
  enumerator += (int)(enumerator * q) >> 12;
  return (int)(0.5 + (enumerator * correction_factor / q));
}

void vp9_save_coding_context(VP9_COMP *cpi) {
  CODING_CONTEXT *const cc = &cpi->coding_context;
  VP9_COMMON *cm = &cpi->common;

  // Stores a snapshot of key state variables which can subsequently be
  // restored with a call to vp9_restore_coding_context. These functions are
  // intended for use in a re-code loop in vp9_compress_frame where the
  // quantizer value is adjusted between loop iterations.
  vp9_copy(cc->nmvjointcost,  cpi->mb.nmvjointcost);
  vp9_copy(cc->nmvcosts,  cpi->mb.nmvcosts);
  vp9_copy(cc->nmvcosts_hp,  cpi->mb.nmvcosts_hp);

  vp9_copy(cc->segment_pred_probs, cm->seg.pred_probs);

  vpx_memcpy(cpi->coding_context.last_frame_seg_map_copy,
             cm->last_frame_seg_map, (cm->mi_rows * cm->mi_cols));

  vp9_copy(cc->last_ref_lf_deltas, cm->lf.last_ref_deltas);
  vp9_copy(cc->last_mode_lf_deltas, cm->lf.last_mode_deltas);

  cc->fc = cm->fc;
}

void vp9_restore_coding_context(VP9_COMP *cpi) {
  CODING_CONTEXT *const cc = &cpi->coding_context;
  VP9_COMMON *cm = &cpi->common;

  // Restore key state variables to the snapshot state stored in the
  // previous call to vp9_save_coding_context.
  vp9_copy(cpi->mb.nmvjointcost, cc->nmvjointcost);
  vp9_copy(cpi->mb.nmvcosts, cc->nmvcosts);
  vp9_copy(cpi->mb.nmvcosts_hp, cc->nmvcosts_hp);

  vp9_copy(cm->seg.pred_probs, cc->segment_pred_probs);

  vpx_memcpy(cm->last_frame_seg_map,
             cpi->coding_context.last_frame_seg_map_copy,
             (cm->mi_rows * cm->mi_cols));

  vp9_copy(cm->lf.last_ref_deltas, cc->last_ref_lf_deltas);
  vp9_copy(cm->lf.last_mode_deltas, cc->last_mode_lf_deltas);

  cm->fc = cc->fc;
}

void vp9_setup_key_frame(VP9_COMP *cpi) {
  VP9_COMMON *cm = &cpi->common;

  vp9_setup_past_independence(cm);

  // interval before next GF
  cpi->rc.frames_till_gf_update_due = cpi->rc.baseline_gf_interval;
  /* All buffers are implicitly updated on key frames. */
  cpi->refresh_golden_frame = 1;
  cpi->refresh_alt_ref_frame = 1;
}

void vp9_setup_inter_frame(VP9_COMP *cpi) {
  VP9_COMMON *cm = &cpi->common;
  if (cm->error_resilient_mode || cm->intra_only)
    vp9_setup_past_independence(cm);

  assert(cm->frame_context_idx < NUM_FRAME_CONTEXTS);
  cm->fc = cm->frame_contexts[cm->frame_context_idx];
}

static int estimate_bits_at_q(int frame_kind, int q, int mbs,
                              double correction_factor) {
  const int bpm = (int)(vp9_bits_per_mb(frame_kind, q, correction_factor));

  // Attempt to retain reasonable accuracy without overflow. The cutoff is
  // chosen such that the maximum product of Bpm and MBs fits 31 bits. The
  // largest Bpm takes 20 bits.
  return (mbs > (1 << 11)) ? (bpm >> BPER_MB_NORMBITS) * mbs
                           : (bpm * mbs) >> BPER_MB_NORMBITS;
}


static void calc_iframe_target_size(VP9_COMP *cpi) {
  // boost defaults to half second
  int target;

  // Clear down mmx registers to allow floating point in what follows
  vp9_clear_system_state();  // __asm emms;

  // New Two pass RC
  target = cpi->rc.per_frame_bandwidth;

  if (cpi->oxcf.rc_max_intra_bitrate_pct) {
    int max_rate = cpi->rc.per_frame_bandwidth
                 * cpi->oxcf.rc_max_intra_bitrate_pct / 100;

    if (target > max_rate)
      target = max_rate;
  }

  cpi->rc.this_frame_target = target;

  // Target rate per SB64 (including partial SB64s.
  cpi->rc.sb64_target_rate = (cpi->rc.this_frame_target * 64 * 64) /
                             (cpi->common.width * cpi->common.height);
}


//  Do the best we can to define the parameters for the next GF based
//  on what information we have available.
//
//  In this experimental code only two pass is supported
//  so we just use the interval determined in the two pass code.
static void calc_gf_params(VP9_COMP *cpi) {
  // Set the gf interval
  cpi->rc.frames_till_gf_update_due = cpi->rc.baseline_gf_interval;
}


static void calc_pframe_target_size(VP9_COMP *cpi) {
  const int min_frame_target = MAX(cpi->rc.min_frame_bandwidth,
                                   cpi->rc.av_per_frame_bandwidth >> 5);
  if (cpi->refresh_alt_ref_frame) {
    // Special alt reference frame case
    // Per frame bit target for the alt ref frame
    cpi->rc.per_frame_bandwidth = cpi->twopass.gf_bits;
    cpi->rc.this_frame_target = cpi->rc.per_frame_bandwidth;
  } else {
    // Normal frames (gf,and inter)
    cpi->rc.this_frame_target = cpi->rc.per_frame_bandwidth;
  }

  // Target rate per SB64 (including partial SB64s.
  cpi->rc.sb64_target_rate = (cpi->rc.this_frame_target * 64 * 64) /
                             (cpi->common.width * cpi->common.height);


  // Check that the total sum of adjustments is not above the maximum allowed.
  // That is, having allowed for the KF and GF penalties, we have not pushed
  // the current inter-frame target too low. If the adjustment we apply here is
  // not capable of recovering all the extra bits we have spent in the KF or GF,
  // then the remainder will have to be recovered over a longer time span via
  // other buffer / rate control mechanisms.
  if (cpi->rc.this_frame_target < min_frame_target)
    cpi->rc.this_frame_target = min_frame_target;

  // Adjust target frame size for Golden Frames:
  if (cpi->rc.frames_till_gf_update_due == 0) {
    cpi->refresh_golden_frame = 1;
    calc_gf_params(cpi);
    // If we are using alternate ref instead of gf then do not apply the boost
    // It will instead be applied to the altref update
    // Jims modified boost
    if (!cpi->source_alt_ref_active) {
      // The spend on the GF is defined in the two pass code
      // for two pass encodes
      cpi->rc.this_frame_target = cpi->rc.per_frame_bandwidth;
    } else {
      // If there is an active ARF at this location use the minimum
      // bits on this frame even if it is a constructed arf.
      // The active maximum quantizer insures that an appropriate
      // number of bits will be spent if needed for constructed ARFs.
      cpi->rc.this_frame_target = 0;
    }
  }
}


void vp9_update_rate_correction_factors(VP9_COMP *cpi, int damp_var) {
  const int q = cpi->common.base_qindex;
  int correction_factor = 100;
  double rate_correction_factor;
  double adjustment_limit;

  int projected_size_based_on_q = 0;

  // Clear down mmx registers to allow floating point in what follows
  vp9_clear_system_state();  // __asm emms;

  if (cpi->common.frame_type == KEY_FRAME) {
    rate_correction_factor = cpi->rc.key_frame_rate_correction_factor;
  } else {
    if (cpi->refresh_alt_ref_frame || cpi->refresh_golden_frame)
      rate_correction_factor = cpi->rc.gf_rate_correction_factor;
    else
      rate_correction_factor = cpi->rc.rate_correction_factor;
  }

  // Work out how big we would have expected the frame to be at this Q given
  // the current correction factor.
  // Stay in double to avoid int overflow when values are large
  projected_size_based_on_q = estimate_bits_at_q(cpi->common.frame_type, q,
                                                 cpi->common.MBs,
                                                 rate_correction_factor);

  // Work out a size correction factor.
  if (projected_size_based_on_q > 0)
    correction_factor =
        (100 * cpi->rc.projected_frame_size) / projected_size_based_on_q;

  // More heavily damped adjustment used if we have been oscillating either side
  // of target.
  switch (damp_var) {
    case 0:
      adjustment_limit = 0.75;
      break;
    case 1:
      adjustment_limit = 0.375;
      break;
    case 2:
    default:
      adjustment_limit = 0.25;
      break;
  }

  // if ( (correction_factor > 102) && (Q < cpi->rc.active_worst_quality) )
  if (correction_factor > 102) {
    // We are not already at the worst allowable quality
    correction_factor =
        (int)(100 + ((correction_factor - 100) * adjustment_limit));
    rate_correction_factor =
        ((rate_correction_factor * correction_factor) / 100);

    // Keep rate_correction_factor within limits
    if (rate_correction_factor > MAX_BPB_FACTOR)
      rate_correction_factor = MAX_BPB_FACTOR;
  } else if (correction_factor < 99) {
    // We are not already at the best allowable quality
    correction_factor =
        (int)(100 - ((100 - correction_factor) * adjustment_limit));
    rate_correction_factor =
        ((rate_correction_factor * correction_factor) / 100);

    // Keep rate_correction_factor within limits
    if (rate_correction_factor < MIN_BPB_FACTOR)
      rate_correction_factor = MIN_BPB_FACTOR;
  }

  if (cpi->common.frame_type == KEY_FRAME) {
    cpi->rc.key_frame_rate_correction_factor = rate_correction_factor;
  } else {
    if (cpi->refresh_alt_ref_frame || cpi->refresh_golden_frame)
      cpi->rc.gf_rate_correction_factor = rate_correction_factor;
    else
      cpi->rc.rate_correction_factor = rate_correction_factor;
  }
}


int vp9_regulate_q(VP9_COMP *cpi, int target_bits_per_frame) {
  int q = cpi->rc.active_worst_quality;

  int i;
  int last_error = INT_MAX;
  int target_bits_per_mb;
  int bits_per_mb_at_this_q;
  double correction_factor;

  // Select the appropriate correction factor based upon type of frame.
  if (cpi->common.frame_type == KEY_FRAME) {
    correction_factor = cpi->rc.key_frame_rate_correction_factor;
  } else {
    if (cpi->refresh_alt_ref_frame || cpi->refresh_golden_frame)
      correction_factor = cpi->rc.gf_rate_correction_factor;
    else
      correction_factor = cpi->rc.rate_correction_factor;
  }

  // Calculate required scaling factor based on target frame size and size of
  // frame produced using previous Q.
  if (target_bits_per_frame >= (INT_MAX >> BPER_MB_NORMBITS))
    target_bits_per_mb =
        (target_bits_per_frame / cpi->common.MBs)
        << BPER_MB_NORMBITS;  // Case where we would overflow int
  else
    target_bits_per_mb =
        (target_bits_per_frame << BPER_MB_NORMBITS) / cpi->common.MBs;

  i = cpi->rc.active_best_quality;

  do {
    bits_per_mb_at_this_q = (int)vp9_bits_per_mb(cpi->common.frame_type, i,
                                                 correction_factor);

    if (bits_per_mb_at_this_q <= target_bits_per_mb) {
      if ((target_bits_per_mb - bits_per_mb_at_this_q) <= last_error)
        q = i;
      else
        q = i - 1;

      break;
    } else {
      last_error = bits_per_mb_at_this_q - target_bits_per_mb;
    }
  } while (++i <= cpi->rc.active_worst_quality);

  return q;
}

static int get_active_quality(int q,
                              int gfu_boost,
                              int low,
                              int high,
                              int *low_motion_minq,
                              int *high_motion_minq) {
  int active_best_quality;
  if (gfu_boost > high) {
    active_best_quality = low_motion_minq[q];
  } else if (gfu_boost < low) {
    active_best_quality = high_motion_minq[q];
  } else {
    const int gap = high - low;
    const int offset = high - gfu_boost;
    const int qdiff = high_motion_minq[q] - low_motion_minq[q];
    const int adjustment = ((offset * qdiff) + (gap >> 1)) / gap;
    active_best_quality = low_motion_minq[q] + adjustment;
  }
  return active_best_quality;
}

int vp9_pick_q_and_adjust_q_bounds(VP9_COMP *cpi,
                                   int * bottom_index, int * top_index) {
  // Set an active best quality and if necessary active worst quality
  int q = cpi->rc.active_worst_quality;
  VP9_COMMON *const cm = &cpi->common;

  if (frame_is_intra_only(cm)) {
#if !CONFIG_MULTIPLE_ARF
    // Handle the special case for key frames forced when we have75 reached
    // the maximum key frame interval. Here force the Q to a range
    // based on the ambient Q to reduce the risk of popping.
    if (cpi->this_key_frame_forced) {
      int delta_qindex;
      int qindex = cpi->rc.last_boosted_qindex;
      double last_boosted_q = vp9_convert_qindex_to_q(qindex);

      delta_qindex = vp9_compute_qdelta(cpi, last_boosted_q,
                                        (last_boosted_q * 0.75));

      cpi->rc.active_best_quality = MAX(qindex + delta_qindex,
                                     cpi->rc.best_quality);
    } else {
      int high = 5000;
      int low = 400;
      double q_adj_factor = 1.0;
      double q_val;

      // Baseline value derived from cpi->active_worst_quality and kf boost
      cpi->rc.active_best_quality = get_active_quality(q, cpi->rc.kf_boost,
                                                    low, high,
                                                    kf_low_motion_minq,
                                                    kf_high_motion_minq);

      // Allow somewhat lower kf minq with small image formats.
      if ((cm->width * cm->height) <= (352 * 288)) {
        q_adj_factor -= 0.25;
      }

      // Make a further adjustment based on the kf zero motion measure.
      q_adj_factor += 0.05 - (0.001 * (double)cpi->kf_zeromotion_pct);

      // Convert the adjustment factor to a qindex delta
      // on active_best_quality.
      q_val = vp9_convert_qindex_to_q(cpi->rc.active_best_quality);
      cpi->rc.active_best_quality +=
          vp9_compute_qdelta(cpi, q_val, (q_val * q_adj_factor));
    }
#else
    double current_q;
    // Force the KF quantizer to be 30% of the active_worst_quality.
    current_q = vp9_convert_qindex_to_q(cpi->rc.active_worst_quality);
    cpi->rc.active_best_quality = cpi->rc.active_worst_quality
        + vp9_compute_qdelta(cpi, current_q, current_q * 0.3);
#endif
  } else if (!cpi->is_src_frame_alt_ref &&
             (cpi->refresh_golden_frame || cpi->refresh_alt_ref_frame)) {
    int high = 2000;
    int low = 400;

    // Use the lower of cpi->rc.active_worst_quality and recent
    // average Q as basis for GF/ARF best Q limit unless last frame was
    // a key frame.
    if (cpi->frames_since_key > 1 &&
        cpi->rc.avg_frame_qindex < cpi->rc.active_worst_quality) {
      q = cpi->rc.avg_frame_qindex;
    }
    // For constrained quality dont allow Q less than the cq level
    if (cpi->oxcf.end_usage == USAGE_CONSTRAINED_QUALITY) {
      if (q < cpi->cq_target_quality)
        q = cpi->cq_target_quality;
      if (cpi->frames_since_key > 1) {
        cpi->rc.active_best_quality = get_active_quality(q, cpi->rc.gfu_boost,
                                                      low, high,
                                                      afq_low_motion_minq,
                                                      afq_high_motion_minq);
      } else {
        cpi->rc.active_best_quality = get_active_quality(q, cpi->rc.gfu_boost,
                                                      low, high,
                                                      gf_low_motion_minq,
                                                      gf_high_motion_minq);
      }
      // Constrained quality use slightly lower active best.
      cpi->rc.active_best_quality = cpi->rc.active_best_quality * 15 / 16;

    } else if (cpi->oxcf.end_usage == USAGE_CONSTANT_QUALITY) {
      if (!cpi->refresh_alt_ref_frame) {
        cpi->rc.active_best_quality = cpi->cq_target_quality;
      } else {
        if (cpi->frames_since_key > 1) {
          cpi->rc.active_best_quality = get_active_quality(q, cpi->rc.gfu_boost,
                                                        low, high,
                                                        afq_low_motion_minq,
                                                        afq_high_motion_minq);
        } else {
          cpi->rc.active_best_quality = get_active_quality(q, cpi->rc.gfu_boost,
                                                        low, high,
                                                        gf_low_motion_minq,
                                                        gf_high_motion_minq);
        }
      }
    } else {
        cpi->rc.active_best_quality = get_active_quality(q, cpi->rc.gfu_boost,
                                                         low, high,
                                                         gf_low_motion_minq,
                                                         gf_high_motion_minq);
    }
  } else {
    if (cpi->oxcf.end_usage == USAGE_CONSTANT_QUALITY) {
      cpi->rc.active_best_quality = cpi->cq_target_quality;
    } else {
      cpi->rc.active_best_quality = inter_minq[q];
      // 1-pass: for now, use the average Q for the active_best, if its lower
      // than active_worst.
      if (cpi->pass == 0 && (cpi->rc.avg_frame_qindex < q))
        cpi->rc.active_best_quality = inter_minq[cpi->rc.avg_frame_qindex];

      // For the constrained quality mode we don't want
      // q to fall below the cq level.
      if ((cpi->oxcf.end_usage == USAGE_CONSTRAINED_QUALITY) &&
          (cpi->rc.active_best_quality < cpi->cq_target_quality)) {
        // If we are strongly undershooting the target rate in the last
        // frames then use the user passed in cq value not the auto
        // cq value.
        if (cpi->rc.rolling_actual_bits < cpi->rc.min_frame_bandwidth)
          cpi->rc.active_best_quality = cpi->oxcf.cq_level;
        else
          cpi->rc.active_best_quality = cpi->cq_target_quality;
      }
    }
  }

  // Clip the active best and worst quality values to limits
  if (cpi->rc.active_worst_quality > cpi->rc.worst_quality)
    cpi->rc.active_worst_quality = cpi->rc.worst_quality;

  if (cpi->rc.active_best_quality < cpi->rc.best_quality)
    cpi->rc.active_best_quality = cpi->rc.best_quality;

  if (cpi->rc.active_best_quality > cpi->rc.worst_quality)
    cpi->rc.active_best_quality = cpi->rc.worst_quality;

  if (cpi->rc.active_worst_quality < cpi->rc.active_best_quality)
    cpi->rc.active_worst_quality = cpi->rc.active_best_quality;

  // Limit Q range for the adaptive loop.
  if (cm->frame_type == KEY_FRAME && !cpi->this_key_frame_forced) {
    *top_index =
      (cpi->rc.active_worst_quality + cpi->rc.active_best_quality * 3) / 4;
    // If this is the first (key) frame in 1-pass, active best is the user
    // best-allowed, and leave the top_index to active_worst.
    if (cpi->pass == 0 && cpi->common.current_video_frame == 0) {
      cpi->rc.active_best_quality = cpi->oxcf.best_allowed_q;
      *top_index = cpi->oxcf.worst_allowed_q;
    }
  } else if (!cpi->is_src_frame_alt_ref &&
             (cpi->oxcf.end_usage != USAGE_STREAM_FROM_SERVER) &&
             (cpi->refresh_golden_frame || cpi->refresh_alt_ref_frame)) {
    *top_index =
      (cpi->rc.active_worst_quality + cpi->rc.active_best_quality) / 2;
  } else {
    *top_index = cpi->rc.active_worst_quality;
  }
  *bottom_index = cpi->rc.active_best_quality;

  if (cpi->oxcf.end_usage == USAGE_CONSTANT_QUALITY) {
    q = cpi->rc.active_best_quality;
  // Special case code to try and match quality with forced key frames
  } else if ((cm->frame_type == KEY_FRAME) && cpi->this_key_frame_forced) {
    q = cpi->rc.last_boosted_qindex;
  } else {
    // Determine initial Q to try.
    if (cpi->pass == 0) {
      // 1-pass: for now, use per-frame-bw for target size of frame, scaled
      // by |x| for key frame.
      int scale = (cm->frame_type == KEY_FRAME) ? 5 : 1;
      q = vp9_regulate_q(cpi, scale * cpi->rc.av_per_frame_bandwidth);
    } else {
      q = vp9_regulate_q(cpi, cpi->rc.this_frame_target);
    }
    if (q > *top_index)
      q = *top_index;
  }

  return q;
}


static int estimate_keyframe_frequency(VP9_COMP *cpi) {
  int i;

  // Average key frame frequency
  int av_key_frame_frequency = 0;

  /* First key frame at start of sequence is a special case. We have no
   * frequency data.
   */
  if (cpi->rc.key_frame_count == 1) {
    /* Assume a default of 1 kf every 2 seconds, or the max kf interval,
     * whichever is smaller.
     */
    int key_freq = cpi->oxcf.key_freq > 0 ? cpi->oxcf.key_freq : 1;
    av_key_frame_frequency = (int)cpi->output_framerate * 2;

    if (cpi->oxcf.auto_key && av_key_frame_frequency > key_freq)
      av_key_frame_frequency = cpi->oxcf.key_freq;

    cpi->rc.prior_key_frame_distance[KEY_FRAME_CONTEXT - 1]
      = av_key_frame_frequency;
  } else {
    unsigned int total_weight = 0;
    int last_kf_interval =
      (cpi->frames_since_key > 0) ? cpi->frames_since_key : 1;

    /* reset keyframe context and calculate weighted average of last
     * KEY_FRAME_CONTEXT keyframes
     */
    for (i = 0; i < KEY_FRAME_CONTEXT; i++) {
      if (i < KEY_FRAME_CONTEXT - 1)
        cpi->rc.prior_key_frame_distance[i]
          = cpi->rc.prior_key_frame_distance[i + 1];
      else
        cpi->rc.prior_key_frame_distance[i] = last_kf_interval;

      av_key_frame_frequency += prior_key_frame_weight[i]
                                * cpi->rc.prior_key_frame_distance[i];
      total_weight += prior_key_frame_weight[i];
    }

    av_key_frame_frequency /= total_weight;
  }
  return av_key_frame_frequency;
}


void vp9_adjust_key_frame_context(VP9_COMP *cpi) {
  // Clear down mmx registers to allow floating point in what follows
  vp9_clear_system_state();

  cpi->frames_since_key = 0;
  cpi->rc.key_frame_count++;
}


void vp9_compute_frame_size_bounds(VP9_COMP *cpi, int *frame_under_shoot_limit,
                                   int *frame_over_shoot_limit) {
  // Set-up bounds on acceptable frame size:
  if (cpi->oxcf.end_usage == USAGE_CONSTANT_QUALITY) {
    *frame_under_shoot_limit = 0;
    *frame_over_shoot_limit  = INT_MAX;
  } else {
    if (cpi->common.frame_type == KEY_FRAME) {
      *frame_over_shoot_limit  = cpi->rc.this_frame_target * 9 / 8;
      *frame_under_shoot_limit = cpi->rc.this_frame_target * 7 / 8;
    } else {
      if (cpi->refresh_alt_ref_frame || cpi->refresh_golden_frame) {
        *frame_over_shoot_limit  = cpi->rc.this_frame_target * 9 / 8;
        *frame_under_shoot_limit = cpi->rc.this_frame_target * 7 / 8;
      } else {
        // Stron overshoot limit for constrained quality
        if (cpi->oxcf.end_usage == USAGE_CONSTRAINED_QUALITY) {
          *frame_over_shoot_limit  = cpi->rc.this_frame_target * 11 / 8;
          *frame_under_shoot_limit = cpi->rc.this_frame_target * 2 / 8;
        } else {
          *frame_over_shoot_limit  = cpi->rc.this_frame_target * 11 / 8;
          *frame_under_shoot_limit = cpi->rc.this_frame_target * 5 / 8;
        }
      }
    }

    // For very small rate targets where the fractional adjustment
    // (eg * 7/8) may be tiny make sure there is at least a minimum
    // range.
    *frame_over_shoot_limit += 200;
    *frame_under_shoot_limit -= 200;
    if (*frame_under_shoot_limit < 0)
      *frame_under_shoot_limit = 0;
  }
}


// return of 0 means drop frame
int vp9_pick_frame_size(VP9_COMP *cpi) {
  VP9_COMMON *cm = &cpi->common;

  if (cm->frame_type == KEY_FRAME)
    calc_iframe_target_size(cpi);
  else
    calc_pframe_target_size(cpi);

  return 1;
}
