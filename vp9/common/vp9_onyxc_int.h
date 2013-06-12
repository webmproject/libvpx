/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_ONYXC_INT_H_
#define VP9_COMMON_VP9_ONYXC_INT_H_

#include "vpx_config.h"
#include "vpx/internal/vpx_codec_internal.h"
#include "vp9_rtcd.h"
#include "vp9/common/vp9_loopfilter.h"
#include "vp9/common/vp9_entropymv.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_quant_common.h"

#if CONFIG_POSTPROC
#include "vp9/common/vp9_postproc.h"
#endif

/* Create/destroy static data structures. */

// Define the number of candidate reference buffers.
#define NUM_REF_FRAMES 8
#define NUM_REF_FRAMES_LG2 3

#define ALLOWED_REFS_PER_FRAME 3

// 1 scratch frame for the new frame, 3 for scaled references on the encoder
// TODO(jkoleszar): These 3 extra references could probably come from the
// normal reference pool.
#define NUM_YV12_BUFFERS (NUM_REF_FRAMES + 4)

#define NUM_FRAME_CONTEXTS_LG2 2
#define NUM_FRAME_CONTEXTS (1 << NUM_FRAME_CONTEXTS_LG2)

#define MAX_LAG_BUFFERS 25

typedef struct frame_contexts {
  vp9_prob y_mode_prob[BLOCK_SIZE_GROUPS][VP9_INTRA_MODES - 1];
  vp9_prob uv_mode_prob[VP9_INTRA_MODES][VP9_INTRA_MODES - 1];
  vp9_prob partition_prob[NUM_FRAME_TYPES][NUM_PARTITION_CONTEXTS]
                         [PARTITION_TYPES - 1];

  nmv_context nmvc;
  nmv_context pre_nmvc;
  /* interframe intra mode probs */
  vp9_prob pre_y_mode_prob[BLOCK_SIZE_GROUPS][VP9_INTRA_MODES - 1];
  vp9_prob pre_uv_mode_prob[VP9_INTRA_MODES][VP9_INTRA_MODES - 1];
  vp9_prob pre_partition_prob[NUM_PARTITION_CONTEXTS][PARTITION_TYPES - 1];
  /* interframe intra mode probs */
  unsigned int y_mode_counts[BLOCK_SIZE_GROUPS][VP9_INTRA_MODES];
  unsigned int uv_mode_counts[VP9_INTRA_MODES][VP9_INTRA_MODES];
  unsigned int partition_counts[NUM_PARTITION_CONTEXTS][PARTITION_TYPES];

  vp9_coeff_probs_model coef_probs[TX_SIZE_MAX_SB][BLOCK_TYPES];
  vp9_coeff_probs_model pre_coef_probs[TX_SIZE_MAX_SB][BLOCK_TYPES];
  vp9_coeff_count_model coef_counts[TX_SIZE_MAX_SB][BLOCK_TYPES];
  unsigned int eob_branch_counts[TX_SIZE_MAX_SB][BLOCK_TYPES][REF_TYPES]
                                [COEF_BANDS][PREV_COEF_CONTEXTS];

  nmv_context_counts NMVcount;
  vp9_prob switchable_interp_prob[VP9_SWITCHABLE_FILTERS + 1]
                                 [VP9_SWITCHABLE_FILTERS - 1];
  vp9_prob pre_switchable_interp_prob[VP9_SWITCHABLE_FILTERS + 1]
      [VP9_SWITCHABLE_FILTERS - 1];
  unsigned int switchable_interp_count[VP9_SWITCHABLE_FILTERS + 1]
                                      [VP9_SWITCHABLE_FILTERS];

  vp9_prob inter_mode_probs[INTER_MODE_CONTEXTS][VP9_INTER_MODES - 1];
  vp9_prob pre_inter_mode_probs[INTER_MODE_CONTEXTS][VP9_INTER_MODES - 1];
  unsigned int inter_mode_counts[INTER_MODE_CONTEXTS][VP9_INTER_MODES - 1][2];

  vp9_prob intra_inter_prob[INTRA_INTER_CONTEXTS];
  vp9_prob comp_inter_prob[COMP_INTER_CONTEXTS];
  vp9_prob single_ref_prob[REF_CONTEXTS][2];
  vp9_prob comp_ref_prob[REF_CONTEXTS];
  vp9_prob pre_intra_inter_prob[INTRA_INTER_CONTEXTS];
  vp9_prob pre_comp_inter_prob[COMP_INTER_CONTEXTS];
  vp9_prob pre_single_ref_prob[REF_CONTEXTS][2];
  vp9_prob pre_comp_ref_prob[REF_CONTEXTS];
  unsigned int intra_inter_count[INTRA_INTER_CONTEXTS][2];
  unsigned int comp_inter_count[COMP_INTER_CONTEXTS][2];
  unsigned int single_ref_count[REF_CONTEXTS][2][2];
  unsigned int comp_ref_count[REF_CONTEXTS][2];

  vp9_prob tx_probs_32x32p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB - 1];
  vp9_prob tx_probs_16x16p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB - 2];
  vp9_prob tx_probs_8x8p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB - 3];
  vp9_prob pre_tx_probs_32x32p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB - 1];
  vp9_prob pre_tx_probs_16x16p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB - 2];
  vp9_prob pre_tx_probs_8x8p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB - 3];
  unsigned int tx_count_32x32p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB];
  unsigned int tx_count_16x16p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB - 1];
  unsigned int tx_count_8x8p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB - 2];

  vp9_prob mbskip_probs[MBSKIP_CONTEXTS];
  vp9_prob pre_mbskip_probs[MBSKIP_CONTEXTS];
  unsigned int mbskip_count[MBSKIP_CONTEXTS][2];
} FRAME_CONTEXT;

typedef enum {
  SINGLE_PREDICTION_ONLY = 0,
  COMP_PREDICTION_ONLY   = 1,
  HYBRID_PREDICTION      = 2,
  NB_PREDICTION_TYPES    = 3,
} COMPPREDMODE_TYPE;

typedef enum {
  ONLY_4X4            = 0,
  ALLOW_8X8           = 1,
  ALLOW_16X16         = 2,
  ALLOW_32X32         = 3,
  TX_MODE_SELECT      = 4,
  NB_TXFM_MODES       = 5,
} TXFM_MODE;

typedef struct VP9Common {
  struct vpx_internal_error_info  error;

  DECLARE_ALIGNED(16, int16_t, y_dequant[QINDEX_RANGE][2]);
  DECLARE_ALIGNED(16, int16_t, uv_dequant[QINDEX_RANGE][2]);
#if CONFIG_ALPHA
  DECLARE_ALIGNED(16, int16_t, a_dequant[QINDEX_RANGE][2]);
#endif

  int width;
  int height;
  int display_width;
  int display_height;
  int last_width;
  int last_height;

  // TODO(jkoleszar): this implies chroma ss right now, but could vary per
  // plane. Revisit as part of the future change to YV12_BUFFER_CONFIG to
  // support additional planes.
  int subsampling_x;
  int subsampling_y;

  YUV_TYPE clr_type;

  YV12_BUFFER_CONFIG *frame_to_show;

  YV12_BUFFER_CONFIG yv12_fb[NUM_YV12_BUFFERS];
  int fb_idx_ref_cnt[NUM_YV12_BUFFERS]; /* reference counts */
  int ref_frame_map[NUM_REF_FRAMES]; /* maps fb_idx to reference slot */

  // TODO(jkoleszar): could expand active_ref_idx to 4, with 0 as intra, and
  // roll new_fb_idx into it.

  // Each frame can reference ALLOWED_REFS_PER_FRAME buffers
  int active_ref_idx[ALLOWED_REFS_PER_FRAME];
  struct scale_factors active_ref_scale[ALLOWED_REFS_PER_FRAME];
  int new_fb_idx;


  YV12_BUFFER_CONFIG post_proc_buffer;
  YV12_BUFFER_CONFIG temp_scale_frame;


  FRAME_TYPE last_frame_type;  /* Save last frame's frame type for motion search. */
  FRAME_TYPE frame_type;

  int show_frame;
  int last_show_frame;

  // Flag signaling that the frame is encoded using only INTRA modes.
  int intra_only;

  // Flag signaling that the frame context should be reset to default values.
  // 0 or 1 implies don't reset, 2 reset just the context specified in the
  // frame header, 3 reset all contexts.
  int reset_frame_context;

  int frame_flags;
  // MBs, mb_rows/cols is in 16-pixel units; mi_rows/cols is in
  // MODE_INFO (8-pixel) units.
  int MBs;
  int mb_rows, mi_rows;
  int mb_cols, mi_cols;
  int mode_info_stride;

  /* profile settings */
  TXFM_MODE txfm_mode;

  int base_qindex;
  int last_kf_gf_q;  /* Q used on the last GF or KF */

  int y_dc_delta_q;
  int uv_dc_delta_q;
  int uv_ac_delta_q;
#if CONFIG_ALPHA
  int a_dc_delta_q;
  int a_ac_delta_q;
#endif

  unsigned int frames_since_golden;
  unsigned int frames_till_alt_ref_frame;

  /* We allocate a MODE_INFO struct for each macroblock, together with
     an extra row on top and column on the left to simplify prediction. */

  MODE_INFO *mip; /* Base of allocated array */
  MODE_INFO *mi;  /* Corresponds to upper left visible macroblock */
  MODE_INFO *prev_mip; /* MODE_INFO array 'mip' from last decoded frame */
  MODE_INFO *prev_mi;  /* 'mi' from last frame (points into prev_mip) */


  // Persistent mb segment id map used in prediction.
  unsigned char *last_frame_seg_map;

  INTERPOLATIONFILTERTYPE mcomp_filter_type;

  loop_filter_info_n lf_info;

  int filter_level;
  int last_sharpness_level;
  int sharpness_level;

  int refresh_frame_context;    /* Two state 0 = NO, 1 = YES */

  int ref_frame_sign_bias[MAX_REF_FRAMES];    /* Two state 0, 1 */

  /* Y,U,V */
  ENTROPY_CONTEXT *above_context[MAX_MB_PLANE];
  ENTROPY_CONTEXT left_context[MAX_MB_PLANE][16];

  // partition contexts
  PARTITION_CONTEXT *above_seg_context;
  PARTITION_CONTEXT left_seg_context[8];

  /* keyframe block modes are predicted by their above, left neighbors */

  vp9_prob kf_y_mode_prob[VP9_INTRA_MODES]
                         [VP9_INTRA_MODES]
                         [VP9_INTRA_MODES - 1];
  vp9_prob kf_uv_mode_prob[VP9_INTRA_MODES] [VP9_INTRA_MODES - 1];

  // Context probabilities when using predictive coding of segment id
  vp9_prob segment_pred_probs[PREDICTION_PROBS];
  unsigned char temporal_update;

  // Context probabilities for reference frame prediction
  int allow_comp_inter_inter;
  MV_REFERENCE_FRAME comp_fixed_ref;
  MV_REFERENCE_FRAME comp_var_ref[2];
  COMPPREDMODE_TYPE comp_pred_mode;

  FRAME_CONTEXT fc;  /* this frame entropy */
  FRAME_CONTEXT frame_contexts[NUM_FRAME_CONTEXTS];
  unsigned int  frame_context_idx; /* Context to use/update */

  unsigned int current_video_frame;
  int near_boffset[3];
  int version;

  double bitrate;
  double framerate;

#if CONFIG_POSTPROC
  struct postproc_state  postproc_state;
#endif

  int error_resilient_mode;
  int frame_parallel_decoding_mode;

  int tile_columns, log2_tile_columns;
  int cur_tile_mi_col_start, cur_tile_mi_col_end, cur_tile_col_idx;
  int tile_rows, log2_tile_rows;
  int cur_tile_mi_row_start, cur_tile_mi_row_end, cur_tile_row_idx;
} VP9_COMMON;

static int get_free_fb(VP9_COMMON *cm) {
  int i;
  for (i = 0; i < NUM_YV12_BUFFERS; i++)
    if (cm->fb_idx_ref_cnt[i] == 0)
      break;

  assert(i < NUM_YV12_BUFFERS);
  cm->fb_idx_ref_cnt[i] = 1;
  return i;
}

static void ref_cnt_fb(int *buf, int *idx, int new_idx) {
  if (buf[*idx] > 0)
    buf[*idx]--;

  *idx = new_idx;

  buf[new_idx]++;
}

static int mi_cols_aligned_to_sb(VP9_COMMON *cm) {
  return 2 * ((cm->mb_cols + 3) & ~3);
}

static INLINE void set_partition_seg_context(VP9_COMMON *cm,
                                             MACROBLOCKD *xd,
                                             int mi_row, int mi_col) {
  xd->above_seg_context = cm->above_seg_context + mi_col;
  xd->left_seg_context  = cm->left_seg_context + (mi_row & MI_MASK);
}

static int check_bsize_coverage(VP9_COMMON *cm, MACROBLOCKD *xd,
                                int mi_row, int mi_col,
                                BLOCK_SIZE_TYPE bsize) {
  int bsl = mi_width_log2(bsize), bs = 1 << bsl;
  int ms = bs / 2;

  if ((mi_row + ms < cm->mi_rows) && (mi_col + ms < cm->mi_cols))
    return 0;

  // frame width/height are multiples of 8, hence 8x8 block should always
  // pass the above check
  assert(bsize > BLOCK_SIZE_SB8X8);

  // return the node index in the prob tree for binary coding
  // only allow horizontal/split partition types
  if ((mi_col + ms < cm->mi_cols) && (mi_row + ms >= cm->mi_rows))
    return 1;
  // only allow vertical/split partition types
  if ((mi_row + ms < cm->mi_rows) && (mi_col + ms >= cm->mi_cols))
    return 2;

  return -1;
}

static void set_mi_row_col(VP9_COMMON *cm, MACROBLOCKD *xd,
                       int mi_row, int bh,
                       int mi_col, int bw) {
  xd->mb_to_top_edge    = -((mi_row * MI_SIZE) << 3);
  xd->mb_to_bottom_edge = ((cm->mi_rows - bh - mi_row) * MI_SIZE) << 3;
  xd->mb_to_left_edge   = -((mi_col * MI_SIZE) << 3);
  xd->mb_to_right_edge  = ((cm->mi_cols - bw - mi_col) * MI_SIZE) << 3;

  // Are edges available for intra prediction?
  xd->up_available    = (mi_row != 0);
  xd->left_available  = (mi_col > cm->cur_tile_mi_col_start);
  xd->right_available = (mi_col + bw < cm->cur_tile_mi_col_end);
}

static int get_mi_row(const MACROBLOCKD *xd) {
  return ((-xd->mb_to_top_edge) >> (3 + LOG2_MI_SIZE));
}

static int get_mi_col(const MACROBLOCKD *xd) {
  return ((-xd->mb_to_left_edge) >> (3 + LOG2_MI_SIZE));
}

static int get_token_alloc(int mb_rows, int mb_cols) {
  return mb_rows * mb_cols * (48 * 16 + 4);
}

static void set_prev_mi(VP9_COMMON *cm) {
  const int use_prev_in_find_mv_refs = cm->width == cm->last_width &&
                                       cm->height == cm->last_height &&
                                       !cm->error_resilient_mode &&
                                       !cm->intra_only &&
                                       cm->last_show_frame;
  // Special case: set prev_mi to NULL when the previous mode info
  // context cannot be used.
  cm->prev_mi = use_prev_in_find_mv_refs ?
                  cm->prev_mip + cm->mode_info_stride + 1 : NULL;
}
#endif  // VP9_COMMON_VP9_ONYXC_INT_H_
