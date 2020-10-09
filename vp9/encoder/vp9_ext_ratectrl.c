/*
 *  Copyright (c) 2020 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/encoder/vp9_ext_ratectrl.h"
#include "vp9/common/vp9_common.h"

void vp9_extrc_init(EXT_RATECTRL *ext_ratectrl) { vp9_zero(*ext_ratectrl); }

void vp9_extrc_create(vpx_rc_funcs_t funcs, vpx_rc_config_t ratectrl_config,
                      EXT_RATECTRL *ext_ratectrl) {
  vpx_rc_firstpass_stats_t *rc_firstpass_stats;
  vp9_extrc_delete(ext_ratectrl);
  ext_ratectrl->funcs = funcs;
  ext_ratectrl->ratectrl_config = ratectrl_config;
  ext_ratectrl->funcs.create_model(ext_ratectrl->funcs.priv,
                                   &ext_ratectrl->ratectrl_config,
                                   ext_ratectrl->model);
  rc_firstpass_stats = &ext_ratectrl->rc_firstpass_stats;
  rc_firstpass_stats->num_frames = ratectrl_config.show_frame_count;
  rc_firstpass_stats->frame_stats =
      vpx_malloc(sizeof(*rc_firstpass_stats->frame_stats) *
                 rc_firstpass_stats->num_frames);
  ext_ratectrl->ready = 1;
}

void vp9_extrc_delete(EXT_RATECTRL *ext_ratectrl) {
  if (ext_ratectrl->ready) {
    ext_ratectrl->funcs.delete_model(ext_ratectrl->model);
  }
  vp9_extrc_init(ext_ratectrl);
}

static void gen_rc_firstpass_stats(const FIRSTPASS_STATS *stats,
                                   double *rc_frame_stats) {
  rc_frame_stats[0] = stats->frame;
  rc_frame_stats[1] = stats->weight;
  rc_frame_stats[2] = stats->intra_error;
  rc_frame_stats[3] = stats->coded_error;
  rc_frame_stats[4] = stats->sr_coded_error;
  rc_frame_stats[5] = stats->frame_noise_energy;
  rc_frame_stats[6] = stats->pcnt_inter;
  rc_frame_stats[7] = stats->pcnt_motion;
  rc_frame_stats[8] = stats->pcnt_second_ref;
  rc_frame_stats[9] = stats->pcnt_neutral;
  rc_frame_stats[10] = stats->pcnt_intra_low;
  rc_frame_stats[11] = stats->pcnt_intra_high;
  rc_frame_stats[12] = stats->intra_skip_pct;
  rc_frame_stats[13] = stats->intra_smooth_pct;
  rc_frame_stats[14] = stats->inactive_zone_rows;
  rc_frame_stats[15] = stats->inactive_zone_cols;
  rc_frame_stats[16] = stats->MVr;
  rc_frame_stats[17] = stats->mvr_abs;
  rc_frame_stats[18] = stats->MVc;
  rc_frame_stats[19] = stats->mvc_abs;
  rc_frame_stats[20] = stats->MVrv;
  rc_frame_stats[21] = stats->MVcv;
  rc_frame_stats[22] = stats->mv_in_out_count;
  rc_frame_stats[23] = stats->duration;
  rc_frame_stats[24] = stats->count;
}

void vp9_extrc_send_firstpass_stats(const FIRST_PASS_INFO *first_pass_info,
                                    EXT_RATECTRL *ext_ratectrl) {
  if (ext_ratectrl->ready) {
    vpx_rc_firstpass_stats_t *rc_firstpass_stats =
        &ext_ratectrl->rc_firstpass_stats;
    int i;
    assert(rc_firstpass_stats->num_frames == first_pass_info->num_frames);
    for (i = 0; i < rc_firstpass_stats->num_frames; ++i) {
      gen_rc_firstpass_stats(&first_pass_info->stats[i],
                             rc_firstpass_stats->frame_stats[i]);
    }
    ext_ratectrl->funcs.send_firstpass_stats(ext_ratectrl->model,
                                             rc_firstpass_stats);
  }
}
