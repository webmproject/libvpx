/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_ENCODER_H_
#define VP9_ENCODER_VP9_ENCODER_H_

#include <stdio.h>

#include "./vpx_config.h"
#include "vpx_ports/mem.h"
#include "vpx/internal/vpx_codec_internal.h"
#include "vpx/vp8cx.h"

#include "vp9/common/vp9_ppflags.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_onyxc_int.h"

#include "vp9/encoder/vp9_aq_cyclicrefresh.h"
#include "vp9/encoder/vp9_context_tree.h"
#include "vp9/encoder/vp9_encodemb.h"
#include "vp9/encoder/vp9_firstpass.h"
#include "vp9/encoder/vp9_lookahead.h"
#include "vp9/encoder/vp9_mbgraph.h"
#include "vp9/encoder/vp9_mcomp.h"
#include "vp9/encoder/vp9_quantize.h"
#include "vp9/encoder/vp9_ratectrl.h"
#include "vp9/encoder/vp9_rd.h"
#include "vp9/encoder/vp9_speed_features.h"
#include "vp9/encoder/vp9_svc_layercontext.h"
#include "vp9/encoder/vp9_tokenize.h"
#include "vp9/encoder/vp9_variance.h"
#if CONFIG_DENOISING
#include "vp9/encoder/vp9_denoiser.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define DEFAULT_GF_INTERVAL         10

typedef struct {
  int nmvjointcost[MV_JOINTS];
  int nmvcosts[2][MV_VALS];
  int nmvcosts_hp[2][MV_VALS];

  vp9_prob segment_pred_probs[PREDICTION_PROBS];

  unsigned char *last_frame_seg_map_copy;

  // 0 = Intra, Last, GF, ARF
  signed char last_ref_lf_deltas[MAX_REF_LF_DELTAS];
  // 0 = ZERO_MV, MV
  signed char last_mode_lf_deltas[MAX_MODE_LF_DELTAS];

  FRAME_CONTEXT fc;
} CODING_CONTEXT;


typedef enum {
  // encode_breakout is disabled.
  ENCODE_BREAKOUT_DISABLED = 0,
  // encode_breakout is enabled.
  ENCODE_BREAKOUT_ENABLED = 1,
  // encode_breakout is enabled with small max_thresh limit.
  ENCODE_BREAKOUT_LIMITED = 2
} ENCODE_BREAKOUT_TYPE;

typedef enum {
  NORMAL      = 0,
  FOURFIVE    = 1,
  THREEFIVE   = 2,
  ONETWO      = 3
} VPX_SCALING;

typedef enum {
  // Good Quality Fast Encoding. The encoder balances quality with the
  // amount of time it takes to encode the output. (speed setting
  // controls how fast)
  ONE_PASS_GOOD = 1,

  // One Pass - Best Quality. The encoder places priority on the
  // quality of the output over encoding speed. The output is compressed
  // at the highest possible quality. This option takes the longest
  // amount of time to encode. (speed setting ignored)
  ONE_PASS_BEST = 2,

  // Two Pass - First Pass. The encoder generates a file of statistics
  // for use in the second encoding pass. (speed setting controls how fast)
  TWO_PASS_FIRST = 3,

  // Two Pass - Second Pass. The encoder uses the statistics that were
  // generated in the first encoding pass to create the compressed
  // output. (speed setting controls how fast)
  TWO_PASS_SECOND_GOOD = 4,

  // Two Pass - Second Pass Best.  The encoder uses the statistics that
  // were generated in the first encoding pass to create the compressed
  // output using the highest possible quality, and taking a
  // longer amount of time to encode. (speed setting ignored)
  TWO_PASS_SECOND_BEST = 5,

  // Realtime/Live Encoding. This mode is optimized for realtime
  // encoding (for example, capturing a television signal or feed from
  // a live camera). (speed setting controls how fast)
  REALTIME = 6,
} MODE;

typedef enum {
  FRAMEFLAGS_KEY    = 1 << 0,
  FRAMEFLAGS_GOLDEN = 1 << 1,
  FRAMEFLAGS_ALTREF = 1 << 2,
} FRAMETYPE_FLAGS;

typedef enum {
  NO_AQ = 0,
  VARIANCE_AQ = 1,
  COMPLEXITY_AQ = 2,
  CYCLIC_REFRESH_AQ = 3,
  AQ_MODE_COUNT  // This should always be the last member of the enum
} AQ_MODE;


typedef struct VP9EncoderConfig {
  BITSTREAM_PROFILE profile;
  BIT_DEPTH bit_depth;
  int width;  // width of data passed to the compressor
  int height;  // height of data passed to the compressor
  double framerate;  // set to passed in framerate
  int64_t target_bandwidth;  // bandwidth to be used in kilobits per second

  int noise_sensitivity;  // pre processing blur: recommendation 0
  int sharpness;  // sharpening output: recommendation 0:
  int speed;
  unsigned int rc_max_intra_bitrate_pct;

  MODE mode;

  // Key Framing Operations
  int auto_key;  // autodetect cut scenes and set the keyframes
  int key_freq;  // maximum distance to key frame.

  int lag_in_frames;  // how many frames lag before we start encoding

  // ----------------------------------------------------------------
  // DATARATE CONTROL OPTIONS

  // vbr, cbr, constrained quality or constant quality
  enum vpx_rc_mode rc_mode;

  // buffer targeting aggressiveness
  int under_shoot_pct;
  int over_shoot_pct;

  // buffering parameters
  int64_t starting_buffer_level_ms;
  int64_t optimal_buffer_level_ms;
  int64_t maximum_buffer_size_ms;

  // Frame drop threshold.
  int drop_frames_water_mark;

  // controlling quality
  int fixed_q;
  int worst_allowed_q;
  int best_allowed_q;
  int cq_level;
  AQ_MODE aq_mode;  // Adaptive Quantization mode

  // Internal frame size scaling.
  int allow_spatial_resampling;
  int scaled_frame_width;
  int scaled_frame_height;

  // Enable feature to reduce the frame quantization every x frames.
  int frame_periodic_boost;

  // two pass datarate control
  int two_pass_vbrbias;        // two pass datarate control tweaks
  int two_pass_vbrmin_section;
  int two_pass_vbrmax_section;
  // END DATARATE CONTROL OPTIONS
  // ----------------------------------------------------------------

  // Spatial and temporal scalability.
  int ss_number_layers;  // Number of spatial layers.
  int ts_number_layers;  // Number of temporal layers.
  // Bitrate allocation for spatial layers.
  int ss_target_bitrate[VPX_SS_MAX_LAYERS];
  int ss_play_alternate[VPX_SS_MAX_LAYERS];
  // Bitrate allocation (CBR mode) and framerate factor, for temporal layers.
  int ts_target_bitrate[VPX_TS_MAX_LAYERS];
  int ts_rate_decimator[VPX_TS_MAX_LAYERS];

  // these parameters aren't to be used in final build don't use!!!
  int play_alternate;

  int encode_breakout;  // early breakout : for video conf recommend 800

  /* Bitfield defining the error resiliency features to enable.
   * Can provide decodable frames after losses in previous
   * frames and decodable partitions after losses in the same frame.
   */
  unsigned int error_resilient_mode;

  /* Bitfield defining the parallel decoding mode where the
   * decoding in successive frames may be conducted in parallel
   * just by decoding the frame headers.
   */
  unsigned int frame_parallel_decoding_mode;

  int arnr_max_frames;
  int arnr_strength;
  int arnr_type;

  int tile_columns;
  int tile_rows;

  struct vpx_fixed_buf         two_pass_stats_in;
  struct vpx_codec_pkt_list  *output_pkt_list;

#if CONFIG_FP_MB_STATS
  struct vpx_fixed_buf         firstpass_mb_stats_in;
#endif

  vp8e_tuning tuning;
} VP9EncoderConfig;

static INLINE int is_lossless_requested(const VP9EncoderConfig *cfg) {
  return cfg->best_allowed_q == 0 && cfg->worst_allowed_q == 0;
}

static INLINE int is_best_mode(MODE mode) {
  return mode == ONE_PASS_BEST || mode == TWO_PASS_SECOND_BEST;
}

typedef struct VP9_COMP {
  QUANTS quants;
  MACROBLOCK mb;
  VP9_COMMON common;
  VP9EncoderConfig oxcf;
  struct lookahead_ctx    *lookahead;
  struct lookahead_entry  *source;
  struct lookahead_entry  *alt_ref_source;
  struct lookahead_entry  *last_source;

  YV12_BUFFER_CONFIG *Source;
  YV12_BUFFER_CONFIG *Last_Source;  // NULL for first frame and alt_ref frames
  YV12_BUFFER_CONFIG *un_scaled_source;
  YV12_BUFFER_CONFIG scaled_source;
  YV12_BUFFER_CONFIG *unscaled_last_source;
  YV12_BUFFER_CONFIG scaled_last_source;

  int gold_is_last;  // gold same as last frame ( short circuit gold searches)
  int alt_is_last;  // Alt same as last ( short circuit altref search)
  int gold_is_alt;  // don't do both alt and gold search ( just do gold).

  int skippable_frame;

  int scaled_ref_idx[3];
  int lst_fb_idx;
  int gld_fb_idx;
  int alt_fb_idx;

  int refresh_last_frame;
  int refresh_golden_frame;
  int refresh_alt_ref_frame;

  int ext_refresh_frame_flags_pending;
  int ext_refresh_last_frame;
  int ext_refresh_golden_frame;
  int ext_refresh_alt_ref_frame;

  int ext_refresh_frame_context_pending;
  int ext_refresh_frame_context;

  YV12_BUFFER_CONFIG last_frame_uf;

  TOKENEXTRA *tok;
  unsigned int tok_count[4][1 << 6];

  // Ambient reconstruction err target for force key frames
  int ambient_err;

  RD_OPT rd;

  CODING_CONTEXT coding_context;

  int zbin_mode_boost;
  int zbin_mode_boost_enabled;
  int active_arnr_frames;           // <= cpi->oxcf.arnr_max_frames
  int active_arnr_strength;         // <= cpi->oxcf.arnr_max_strength

  int64_t last_time_stamp_seen;
  int64_t last_end_time_stamp_seen;
  int64_t first_time_stamp_ever;

  RATE_CONTROL rc;

  vp9_coeff_count coef_counts[TX_SIZES][PLANE_TYPES];

  struct vpx_codec_pkt_list  *output_pkt_list;

  MBGRAPH_FRAME_STATS mbgraph_stats[MAX_LAG_BUFFERS];
  int mbgraph_n_frames;             // number of frames filled in the above
  int static_mb_pct;                // % forced skip mbs by segmentation

  int pass;

  int ref_frame_flags;

  SPEED_FEATURES sf;

  unsigned int max_mv_magnitude;
  int mv_step_param;

  // Default value is 1. From first pass stats, encode_breakout may be disabled.
  ENCODE_BREAKOUT_TYPE allow_encode_breakout;

  // Get threshold from external input. A suggested threshold is 800 for HD
  // clips, and 300 for < HD clips.
  int encode_breakout;

  unsigned char *segmentation_map;

  // segment threashold for encode breakout
  int  segment_encode_breakout[MAX_SEGMENTS];

  unsigned char *complexity_map;

  CYCLIC_REFRESH *cyclic_refresh;

  fractional_mv_step_fp *find_fractional_mv_step;
  vp9_full_search_fn_t full_search_sad;
  vp9_refining_search_fn_t refining_search_sad;
  vp9_diamond_search_fn_t diamond_search_sad;
  vp9_variance_fn_ptr_t fn_ptr[BLOCK_SIZES];
  uint64_t time_receive_data;
  uint64_t time_compress_data;
  uint64_t time_pick_lpf;
  uint64_t time_encode_sb_row;

#if CONFIG_FP_MB_STATS
  int use_fp_mb_stats;
#endif

  TWO_PASS twopass;

  YV12_BUFFER_CONFIG alt_ref_buffer;
  YV12_BUFFER_CONFIG *frames[MAX_LAG_BUFFERS];

#if CONFIG_INTERNAL_STATS
  unsigned int mode_chosen_counts[MAX_MODES];

  int    count;
  double total_y;
  double total_u;
  double total_v;
  double total;
  uint64_t total_sq_error;
  uint64_t total_samples;

  double totalp_y;
  double totalp_u;
  double totalp_v;
  double totalp;
  uint64_t totalp_sq_error;
  uint64_t totalp_samples;

  int    bytes;
  double summed_quality;
  double summed_weights;
  double summedp_quality;
  double summedp_weights;
  unsigned int tot_recode_hits;


  double total_ssimg_y;
  double total_ssimg_u;
  double total_ssimg_v;
  double total_ssimg_all;

  int b_calculate_ssimg;
#endif
  int b_calculate_psnr;

  int droppable;

  int dummy_packing;    /* flag to indicate if packing is dummy */

  unsigned int tx_stepdown_count[TX_SIZES];

  int initial_width;
  int initial_height;

  int use_svc;

  SVC svc;

  // Store frame variance info in SOURCE_VAR_BASED_PARTITION search type.
  diff *source_diff_var;
  // The threshold used in SOURCE_VAR_BASED_PARTITION search type.
  unsigned int source_var_thresh;
  int frames_till_next_var_check;

  int frame_flags;

  search_site_config ss_cfg;

  int mbmode_cost[INTRA_MODES];
  unsigned inter_mode_cost[INTER_MODE_CONTEXTS][INTER_MODES];
  int intra_uv_mode_cost[FRAME_TYPES][INTRA_MODES];
  int y_mode_costs[INTRA_MODES][INTRA_MODES][INTRA_MODES];
  int switchable_interp_costs[SWITCHABLE_FILTER_CONTEXTS][SWITCHABLE_FILTERS];

  PICK_MODE_CONTEXT *leaf_tree;
  PC_TREE *pc_tree;
  PC_TREE *pc_root;
  int partition_cost[PARTITION_CONTEXTS][PARTITION_TYPES];

  int multi_arf_allowed;
  int multi_arf_enabled;
  int multi_arf_last_grp_enabled;

#if CONFIG_DENOISING
  VP9_DENOISER denoiser;
#endif
} VP9_COMP;

void vp9_initialize_enc();

struct VP9_COMP *vp9_create_compressor(VP9EncoderConfig *oxcf);
void vp9_remove_compressor(VP9_COMP *cpi);

void vp9_change_config(VP9_COMP *cpi, const VP9EncoderConfig *oxcf);

  // receive a frames worth of data. caller can assume that a copy of this
  // frame is made and not just a copy of the pointer..
int vp9_receive_raw_frame(VP9_COMP *cpi, unsigned int frame_flags,
                          YV12_BUFFER_CONFIG *sd, int64_t time_stamp,
                          int64_t end_time_stamp);

int vp9_get_compressed_data(VP9_COMP *cpi, unsigned int *frame_flags,
                            size_t *size, uint8_t *dest,
                            int64_t *time_stamp, int64_t *time_end, int flush);

int vp9_get_preview_raw_frame(VP9_COMP *cpi, YV12_BUFFER_CONFIG *dest,
                              vp9_ppflags_t *flags);

int vp9_use_as_reference(VP9_COMP *cpi, int ref_frame_flags);

void vp9_update_reference(VP9_COMP *cpi, int ref_frame_flags);

int vp9_copy_reference_enc(VP9_COMP *cpi, VP9_REFFRAME ref_frame_flag,
                           YV12_BUFFER_CONFIG *sd);

int vp9_get_reference_enc(VP9_COMP *cpi, int index,
                          YV12_BUFFER_CONFIG **fb);

int vp9_set_reference_enc(VP9_COMP *cpi, VP9_REFFRAME ref_frame_flag,
                          YV12_BUFFER_CONFIG *sd);

int vp9_update_entropy(VP9_COMP *cpi, int update);

int vp9_set_active_map(VP9_COMP *cpi, unsigned char *map, int rows, int cols);

int vp9_set_internal_size(VP9_COMP *cpi,
                          VPX_SCALING horiz_mode, VPX_SCALING vert_mode);

int vp9_set_size_literal(VP9_COMP *cpi, unsigned int width,
                         unsigned int height);

void vp9_set_svc(VP9_COMP *cpi, int use_svc);

int vp9_get_quantizer(struct VP9_COMP *cpi);

static INLINE int get_ref_frame_idx(const VP9_COMP *cpi,
                                    MV_REFERENCE_FRAME ref_frame) {
  if (ref_frame == LAST_FRAME) {
    return cpi->lst_fb_idx;
  } else if (ref_frame == GOLDEN_FRAME) {
    return cpi->gld_fb_idx;
  } else {
    return cpi->alt_fb_idx;
  }
}

static INLINE YV12_BUFFER_CONFIG *get_ref_frame_buffer(
    VP9_COMP *cpi, MV_REFERENCE_FRAME ref_frame) {
  VP9_COMMON * const cm = &cpi->common;
  return &cm->frame_bufs[cm->ref_frame_map[get_ref_frame_idx(cpi, ref_frame)]]
      .buf;
}

// Intra only frames, golden frames (except alt ref overlays) and
// alt ref frames tend to be coded at a higher than ambient quality
static INLINE int frame_is_boosted(const VP9_COMP *cpi) {
  return frame_is_intra_only(&cpi->common) || cpi->refresh_alt_ref_frame ||
         (cpi->refresh_golden_frame && !cpi->rc.is_src_frame_alt_ref) ||
         vp9_is_upper_layer_key_frame(cpi);
}

static INLINE int get_token_alloc(int mb_rows, int mb_cols) {
  // TODO(JBB): double check we can't exceed this token count if we have a
  // 32x32 transform crossing a boundary at a multiple of 16.
  // mb_rows, cols are in units of 16 pixels. We assume 3 planes all at full
  // resolution. We assume up to 1 token per pixel, and then allow
  // a head room of 4.
  return mb_rows * mb_cols * (16 * 16 * 3 + 4);
}

int vp9_get_y_sse(const YV12_BUFFER_CONFIG *a, const YV12_BUFFER_CONFIG *b);

void vp9_alloc_compressor_data(VP9_COMP *cpi);

void vp9_scale_references(VP9_COMP *cpi);

void vp9_update_reference_frames(VP9_COMP *cpi);

int64_t vp9_rescale(int64_t val, int64_t num, int denom);

void vp9_set_high_precision_mv(VP9_COMP *cpi, int allow_high_precision_mv);

YV12_BUFFER_CONFIG *vp9_scale_if_required(VP9_COMMON *cm,
                                          YV12_BUFFER_CONFIG *unscaled,
                                          YV12_BUFFER_CONFIG *scaled);

void vp9_apply_encoding_flags(VP9_COMP *cpi, vpx_enc_frame_flags_t flags);

static INLINE int is_altref_enabled(const VP9_COMP *const cpi) {
  return cpi->oxcf.mode != REALTIME && cpi->oxcf.lag_in_frames > 0 &&
         (cpi->oxcf.play_alternate &&
          (!(cpi->use_svc && cpi->svc.number_temporal_layers == 1) ||
           cpi->oxcf.ss_play_alternate[cpi->svc.spatial_layer_id]));
}

static INLINE void set_ref_ptrs(VP9_COMMON *cm, MACROBLOCKD *xd,
                                MV_REFERENCE_FRAME ref0,
                                MV_REFERENCE_FRAME ref1) {
  xd->block_refs[0] = &cm->frame_refs[ref0 >= LAST_FRAME ? ref0 - LAST_FRAME
                                                         : 0];
  xd->block_refs[1] = &cm->frame_refs[ref1 >= LAST_FRAME ? ref1 - LAST_FRAME
                                                         : 0];
}

static INLINE int get_chessboard_index(const VP9_COMMON *cm) {
  return cm->current_video_frame % 2;
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_ENCODER_H_
