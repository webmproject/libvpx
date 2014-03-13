/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>

#include "vp9/encoder/vp9_onyx_int.h"
#include "vp9/encoder/vp9_svc_layercontext.h"

void vp9_init_layer_context(VP9_COMP *const cpi) {
  const VP9_CONFIG *const oxcf = &cpi->oxcf;
  int temporal_layer = 0;
  cpi->svc.spatial_layer_id = 0;
  cpi->svc.temporal_layer_id = 0;
  for (temporal_layer = 0; temporal_layer < cpi->svc.number_temporal_layers;
      ++temporal_layer) {
    LAYER_CONTEXT *const lc = &cpi->svc.layer_context[temporal_layer];
    RATE_CONTROL *const lrc = &lc->rc;
    lrc->avg_frame_qindex[INTER_FRAME] = q_trans[oxcf->worst_allowed_q];
    lrc->last_q[INTER_FRAME] = q_trans[oxcf->worst_allowed_q];
    lrc->ni_av_qi = q_trans[oxcf->worst_allowed_q];
    lrc->total_actual_bits = 0;
    lrc->total_target_vs_actual = 0;
    lrc->ni_tot_qi = 0;
    lrc->tot_q = 0.0;
    lrc->avg_q = 0.0;
    lrc->ni_frames = 0;
    lrc->decimation_count = 0;
    lrc->decimation_factor = 0;
    lrc->rate_correction_factor = 1.0;
    lrc->key_frame_rate_correction_factor = 1.0;
    lc->target_bandwidth = oxcf->ts_target_bitrate[temporal_layer] *
        1000;
    lrc->buffer_level =
        vp9_rescale((int)(oxcf->starting_buffer_level),
                    lc->target_bandwidth, 1000);
    lrc->bits_off_target = lrc->buffer_level;
  }
}

// Update the layer context from a change_config() call.
void vp9_update_layer_context_change_config(VP9_COMP *const cpi,
                                            const int target_bandwidth) {
  const VP9_CONFIG *const oxcf = &cpi->oxcf;
  const RATE_CONTROL *const rc = &cpi->rc;
  int temporal_layer = 0;
  float bitrate_alloc = 1.0;
  for (temporal_layer = 0; temporal_layer < cpi->svc.number_temporal_layers;
      ++temporal_layer) {
    LAYER_CONTEXT *const lc = &cpi->svc.layer_context[temporal_layer];
    RATE_CONTROL *const lrc = &lc->rc;
    lc->target_bandwidth = oxcf->ts_target_bitrate[temporal_layer] * 1000;
    bitrate_alloc = (float)lc->target_bandwidth / (float)target_bandwidth;
    // Update buffer-related quantities.
    lc->starting_buffer_level =
        (int64_t)(oxcf->starting_buffer_level * bitrate_alloc);
    lc->optimal_buffer_level =
        (int64_t)(oxcf->optimal_buffer_level * bitrate_alloc);
    lc->maximum_buffer_size =
        (int64_t)(oxcf->maximum_buffer_size * bitrate_alloc);
    lrc->bits_off_target = MIN(lrc->bits_off_target, lc->maximum_buffer_size);
    lrc->buffer_level = MIN(lrc->buffer_level, lc->maximum_buffer_size);
    // Update framerate-related quantities.
    lc->framerate = oxcf->framerate / oxcf->ts_rate_decimator[temporal_layer];
    lrc->av_per_frame_bandwidth = (int)(lc->target_bandwidth / lc->framerate);
    lrc->max_frame_bandwidth = rc->max_frame_bandwidth;
    // Update qp-related quantities.
    lrc->worst_quality = rc->worst_quality;
    lrc->best_quality = rc->best_quality;
  }
}

void vp9_update_layer_framerate(VP9_COMP *const cpi) {
  int temporal_layer = cpi->svc.temporal_layer_id;
  const VP9_CONFIG *const oxcf = &cpi->oxcf;
  LAYER_CONTEXT *const lc = &cpi->svc.layer_context[temporal_layer];
  RATE_CONTROL *const lrc = &lc->rc;
  lc->framerate = oxcf->framerate / oxcf->ts_rate_decimator[temporal_layer];
  lrc->av_per_frame_bandwidth = (int)(lc->target_bandwidth / lc->framerate);
  lrc->max_frame_bandwidth = cpi->rc.max_frame_bandwidth;
  // Update the average layer frame size (non-cumulative per-frame-bw).
  if (temporal_layer == 0) {
    lc->avg_frame_size = lrc->av_per_frame_bandwidth;
  } else {
    double prev_layer_framerate = oxcf->framerate /
        oxcf->ts_rate_decimator[temporal_layer - 1];
    int prev_layer_target_bandwidth =
        oxcf->ts_target_bitrate[temporal_layer - 1] * 1000;
    lc->avg_frame_size =
        (int)((lc->target_bandwidth - prev_layer_target_bandwidth) /
              (lc->framerate - prev_layer_framerate));
  }
}

void vp9_restore_layer_context(VP9_COMP *const cpi) {
  int temporal_layer = cpi->svc.temporal_layer_id;
  LAYER_CONTEXT *lc = &cpi->svc.layer_context[temporal_layer];
  int frame_since_key = cpi->rc.frames_since_key;
  int frame_to_key = cpi->rc.frames_to_key;
  cpi->rc = lc->rc;
  cpi->oxcf.target_bandwidth = lc->target_bandwidth;
  cpi->oxcf.starting_buffer_level = lc->starting_buffer_level;
  cpi->oxcf.optimal_buffer_level = lc->optimal_buffer_level;
  cpi->oxcf.maximum_buffer_size = lc->maximum_buffer_size;
  cpi->output_framerate = lc->framerate;
  // Reset the frames_since_key and frames_to_key counters to their values
  // before the layer restore. Keep these defined for the stream (not layer).
  cpi->rc.frames_since_key = frame_since_key;
  cpi->rc.frames_to_key = frame_to_key;
}

void vp9_save_layer_context(VP9_COMP *const cpi) {
  int temporal_layer = cpi->svc.temporal_layer_id;
  LAYER_CONTEXT *lc = &cpi->svc.layer_context[temporal_layer];
  lc->rc = cpi->rc;
  lc->target_bandwidth = (int)cpi->oxcf.target_bandwidth;
  lc->starting_buffer_level = cpi->oxcf.starting_buffer_level;
  lc->optimal_buffer_level = cpi->oxcf.optimal_buffer_level;
  lc->maximum_buffer_size = cpi->oxcf.maximum_buffer_size;
  lc->framerate = cpi->output_framerate;
}
