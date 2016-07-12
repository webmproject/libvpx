/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_ENCODER_H_
#define VP10_ENCODER_ENCODER_H_

#include <stdio.h>

#include "./vpx_config.h"
#include "vpx/vp8cx.h"

#include "vp10/common/alloccommon.h"
#include "vp10/common/ppflags.h"
#include "vp10/common/entropymode.h"
#include "vp10/common/thread_common.h"
#include "vp10/common/onyxc_int.h"

#include "vp10/encoder/aq_cyclicrefresh.h"
#if CONFIG_ANS
#include "vp10/encoder/buf_ans.h"
#endif
#include "vp10/encoder/context_tree.h"
#include "vp10/encoder/encodemb.h"
#include "vp10/encoder/firstpass.h"
#include "vp10/encoder/lookahead.h"
#include "vp10/encoder/mbgraph.h"
#include "vp10/encoder/mcomp.h"
#include "vp10/encoder/quantize.h"
#include "vp10/encoder/ratectrl.h"
#include "vp10/encoder/rd.h"
#include "vp10/encoder/speed_features.h"
#include "vp10/encoder/tokenize.h"
#include "vp10/encoder/variance_tree.h"

#if CONFIG_VP9_TEMPORAL_DENOISING
#include "vp10/encoder/denoiser.h"
#endif

#if CONFIG_INTERNAL_STATS
#include "vpx_dsp/ssim.h"
#endif
#include "vpx_dsp/variance.h"
#include "vpx/internal/vpx_codec_internal.h"
#include "vpx_util/vpx_thread.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int nmvjointcost[MV_JOINTS];
  int nmvcosts[2][MV_VALS];
  int nmvcosts_hp[2][MV_VALS];

#if CONFIG_REF_MV
  int nmv_vec_cost[NMV_CONTEXTS][MV_JOINTS];
  int nmv_costs[NMV_CONTEXTS][2][MV_VALS];
  int nmv_costs_hp[NMV_CONTEXTS][2][MV_VALS];
#endif

  unsigned char *last_frame_seg_map_copy;

  // 0 = Intra, Last, GF, ARF
  signed char last_ref_lf_deltas[MAX_REF_FRAMES];
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
  // Good Quality Fast Encoding. The encoder balances quality with the amount of
  // time it takes to encode the output. Speed setting controls how fast.
  GOOD,

  // The encoder places priority on the quality of the output over encoding
  // speed. The output is compressed at the highest possible quality. This
  // option takes the longest amount of time to encode. Speed setting ignored.
  BEST,

  // Realtime/Live Encoding. This mode is optimized for realtime encoding (for
  // example, capturing a television signal or feed from a live camera). Speed
  // setting controls how fast.
  REALTIME
} MODE;

typedef enum {
  FRAMEFLAGS_KEY    = 1 << 0,
  FRAMEFLAGS_GOLDEN = 1 << 1,
#if CONFIG_EXT_REFS
  FRAMEFLAGS_BWDREF = 1 << 2,
  FRAMEFLAGS_ALTREF = 1 << 3,
#else
  FRAMEFLAGS_ALTREF = 1 << 2,
#endif  // CONFIG_EXT_REFS
} FRAMETYPE_FLAGS;

typedef enum {
  NO_AQ = 0,
  VARIANCE_AQ = 1,
  COMPLEXITY_AQ = 2,
  CYCLIC_REFRESH_AQ = 3,
  AQ_MODE_COUNT  // This should always be the last member of the enum
} AQ_MODE;

typedef enum {
  RESIZE_NONE = 0,    // No frame resizing allowed.
  RESIZE_FIXED = 1,   // All frames are coded at the specified dimension.
  RESIZE_DYNAMIC = 2  // Coded size of each frame is determined by the codec.
} RESIZE_TYPE;

typedef struct VP10EncoderConfig {
  BITSTREAM_PROFILE profile;
  vpx_bit_depth_t bit_depth;     // Codec bit-depth.
  int width;  // width of data passed to the compressor
  int height;  // height of data passed to the compressor
  unsigned int input_bit_depth;  // Input bit depth.
  double init_framerate;  // set to passed in framerate
  int64_t target_bandwidth;  // bandwidth to be used in bits per second

  int noise_sensitivity;  // pre processing blur: recommendation 0
  int sharpness;  // sharpening output: recommendation 0:
  int speed;
  // maximum allowed bitrate for any intra frame in % of bitrate target.
  unsigned int rc_max_intra_bitrate_pct;
  // maximum allowed bitrate for any inter frame in % of bitrate target.
  unsigned int rc_max_inter_bitrate_pct;
  // percent of rate boost for golden frame in CBR mode.
  unsigned int gf_cbr_boost_pct;

  MODE mode;
  int pass;

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
  RESIZE_TYPE resize_mode;
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

  int enable_auto_arf;
#if CONFIG_EXT_REFS
  int enable_auto_brf;  // (b)ackward (r)ef (f)rame
#endif  // CONFIG_EXT_REFS

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

  int min_gf_interval;
  int max_gf_interval;

  int tile_columns;
  int tile_rows;

  int max_threads;

  vpx_fixed_buf_t two_pass_stats_in;
  struct vpx_codec_pkt_list *output_pkt_list;

#if CONFIG_FP_MB_STATS
  vpx_fixed_buf_t firstpass_mb_stats_in;
#endif

  vp8e_tuning tuning;
  vp9e_tune_content content;
#if CONFIG_VP9_HIGHBITDEPTH
  int use_highbitdepth;
#endif
  vpx_color_space_t color_space;
  int color_range;
  int render_width;
  int render_height;

#if CONFIG_EXT_PARTITION
  vpx_superblock_size_t superblock_size;
#endif  // CONFIG_EXT_PARTITION
} VP10EncoderConfig;

static INLINE int is_lossless_requested(const VP10EncoderConfig *cfg) {
  return cfg->best_allowed_q == 0 && cfg->worst_allowed_q == 0;
}

// TODO(jingning) All spatially adaptive variables should go to TileDataEnc.
typedef struct TileDataEnc {
  TileInfo tile_info;
  int thresh_freq_fact[BLOCK_SIZES][MAX_MODES];
  int mode_map[BLOCK_SIZES][MAX_MODES];
} TileDataEnc;

typedef struct RD_COUNTS {
  vp10_coeff_count coef_counts[TX_SIZES][PLANE_TYPES];
  int64_t comp_pred_diff[REFERENCE_MODES];
  int m_search_count;
  int ex_search_count;
} RD_COUNTS;

typedef struct ThreadData {
  MACROBLOCK mb;
  RD_COUNTS rd_counts;
  FRAME_COUNTS *counts;

  PICK_MODE_CONTEXT *leaf_tree;
  PC_TREE *pc_tree;
  PC_TREE *pc_root[MAX_MIB_SIZE_LOG2 - MIN_MIB_SIZE_LOG2 + 1];

  VAR_TREE *var_tree;
  VAR_TREE *var_root[MAX_MIB_SIZE_LOG2 - MIN_MIB_SIZE_LOG2 + 1];
} ThreadData;

struct EncWorkerData;

typedef struct ActiveMap {
  int enabled;
  int update;
  unsigned char *map;
} ActiveMap;

typedef enum {
  Y,
  U,
  V,
  ALL
} STAT_TYPE;

typedef struct IMAGE_STAT {
  double stat[ALL+1];
  double worst;
} ImageStat;

typedef struct {
  int ref_count;
  YV12_BUFFER_CONFIG buf;
} EncRefCntBuffer;

#if CONFIG_ENTROPY
typedef struct SUBFRAME_STATS {
  vp10_coeff_probs_model
  coef_probs_buf[COEF_PROBS_BUFS][TX_SIZES][PLANE_TYPES];
  vp10_coeff_count
  coef_counts_buf[COEF_PROBS_BUFS][TX_SIZES][PLANE_TYPES];
  unsigned int
  eob_counts_buf[COEF_PROBS_BUFS]
                [TX_SIZES][PLANE_TYPES][REF_TYPES][COEF_BANDS][COEFF_CONTEXTS];
  vp10_coeff_probs_model enc_starting_coef_probs[TX_SIZES][PLANE_TYPES];
} SUBFRAME_STATS;
#endif  // CONFIG_ENTROPY

typedef struct TileBufferEnc {
  uint8_t *data;
  size_t size;
} TileBufferEnc;

typedef struct VP10_COMP {
  QUANTS quants;
  ThreadData td;
  MB_MODE_INFO_EXT *mbmi_ext_base;
  DECLARE_ALIGNED(16, int16_t, y_dequant[QINDEX_RANGE][8]);   // 8: SIMD width
  DECLARE_ALIGNED(16, int16_t, uv_dequant[QINDEX_RANGE][8]);  // 8: SIMD width
#if CONFIG_NEW_QUANT
  DECLARE_ALIGNED(16, dequant_val_type_nuq,
                  y_dequant_val_nuq[QUANT_PROFILES][QINDEX_RANGE][COEF_BANDS]);
  DECLARE_ALIGNED(16, dequant_val_type_nuq,
                  uv_dequant_val_nuq[QUANT_PROFILES][QINDEX_RANGE][COEF_BANDS]);
#endif  // CONFIG_NEW_QUANT
  VP10_COMMON common;
  VP10EncoderConfig oxcf;
  struct lookahead_ctx    *lookahead;
  struct lookahead_entry  *alt_ref_source;

  YV12_BUFFER_CONFIG *Source;
  YV12_BUFFER_CONFIG *Last_Source;  // NULL for first frame and alt_ref frames
  YV12_BUFFER_CONFIG *un_scaled_source;
  YV12_BUFFER_CONFIG scaled_source;
  YV12_BUFFER_CONFIG *unscaled_last_source;
  YV12_BUFFER_CONFIG scaled_last_source;

  // Up-sampled reference buffers
  EncRefCntBuffer upsampled_ref_bufs[MAX_REF_FRAMES];
  int upsampled_ref_idx[MAX_REF_FRAMES];

  // For a still frame, this flag is set to 1 to skip partition search.
  int partition_search_skippable_frame;

  int scaled_ref_idx[MAX_REF_FRAMES];
#if CONFIG_EXT_REFS
  int lst_fb_idxes[LAST_REF_FRAMES];
#else
  int lst_fb_idx;
#endif  // CONFIG_EXT_REFS
  int gld_fb_idx;
#if CONFIG_EXT_REFS
  int bwd_fb_idx;  // BWD_REF_FRAME
#endif  // CONFIG_EXT_REFS
  int alt_fb_idx;

  int last_show_frame_buf_idx;  // last show frame buffer index

  int refresh_last_frame;
  int refresh_golden_frame;
#if CONFIG_EXT_REFS
  int refresh_bwd_ref_frame;
#endif  // CONFIG_EXT_REFS
  int refresh_alt_ref_frame;

  int ext_refresh_frame_flags_pending;
  int ext_refresh_last_frame;
  int ext_refresh_golden_frame;
  int ext_refresh_alt_ref_frame;

  int ext_refresh_frame_context_pending;
  int ext_refresh_frame_context;

  YV12_BUFFER_CONFIG last_frame_uf;
#if CONFIG_LOOP_RESTORATION
  YV12_BUFFER_CONFIG last_frame_db;
#endif  // CONFIG_LOOP_RESTORATION

  // Ambient reconstruction err target for force key frames
  int64_t ambient_err;

  RD_OPT rd;

  CODING_CONTEXT coding_context;

#if CONFIG_REF_MV
  int *nmv_costs[NMV_CONTEXTS][2];
  int *nmv_costs_hp[NMV_CONTEXTS][2];
#endif

  int *nmvcosts[2];
  int *nmvcosts_hp[2];
  int *nmvsadcosts[2];
  int *nmvsadcosts_hp[2];

  int64_t last_time_stamp_seen;
  int64_t last_end_time_stamp_seen;
  int64_t first_time_stamp_ever;

  RATE_CONTROL rc;
  double framerate;

  int interp_filter_selected[MAX_REF_FRAMES][SWITCHABLE];

  struct vpx_codec_pkt_list  *output_pkt_list;

  MBGRAPH_FRAME_STATS mbgraph_stats[MAX_LAG_BUFFERS];
  int mbgraph_n_frames;             // number of frames filled in the above
  int static_mb_pct;                // % forced skip mbs by segmentation
  int ref_frame_flags;

  SPEED_FEATURES sf;

  unsigned int max_mv_magnitude;
  int mv_step_param;

  int allow_comp_inter_inter;

  // Default value is 1. From first pass stats, encode_breakout may be disabled.
  ENCODE_BREAKOUT_TYPE allow_encode_breakout;

  // Get threshold from external input. A suggested threshold is 800 for HD
  // clips, and 300 for < HD clips.
  int encode_breakout;

  uint8_t *segmentation_map;

  // segment threashold for encode breakout
  int  segment_encode_breakout[MAX_SEGMENTS];

  CYCLIC_REFRESH *cyclic_refresh;
  ActiveMap active_map;

  fractional_mv_step_fp *find_fractional_mv_step;
  vp10_full_search_fn_t full_search_sad;  // It is currently unused.
  vp10_diamond_search_fn_t diamond_search_sad;
  vp10_variance_fn_ptr_t fn_ptr[BLOCK_SIZES];
  uint64_t time_receive_data;
  uint64_t time_compress_data;
  uint64_t time_pick_lpf;
  uint64_t time_encode_sb_row;

#if CONFIG_FP_MB_STATS
  int use_fp_mb_stats;
#endif

  TWO_PASS twopass;

  YV12_BUFFER_CONFIG alt_ref_buffer;

#if CONFIG_INTERNAL_STATS
  unsigned int mode_chosen_counts[MAX_MODES];

  int    count;
  uint64_t total_sq_error;
  uint64_t total_samples;
  ImageStat psnr;

  double total_blockiness;
  double worst_blockiness;

  int    bytes;
  double summed_quality;
  double summed_weights;
  unsigned int tot_recode_hits;
  double worst_ssim;

  ImageStat fastssim;
  ImageStat psnrhvs;

  int b_calculate_blockiness;
  int b_calculate_consistency;

  double total_inconsistency;
  double worst_consistency;
  Ssimv *ssim_vars;
  Metrics metrics;
#endif
  int b_calculate_psnr;

  int droppable;

  int initial_width;
  int initial_height;
  int initial_mbs;  // Number of MBs in the full-size frame; to be used to
                    // normalize the firstpass stats. This will differ from the
                    // number of MBs in the current frame when the frame is
                    // scaled.

  // Store frame variance info in SOURCE_VAR_BASED_PARTITION search type.
  diff *source_diff_var;
  // The threshold used in SOURCE_VAR_BASED_PARTITION search type.
  unsigned int source_var_thresh;
  int frames_till_next_var_check;

  int frame_flags;

  search_site_config ss_cfg;

  int mbmode_cost[BLOCK_SIZE_GROUPS][INTRA_MODES];
#if CONFIG_REF_MV
  int newmv_mode_cost[NEWMV_MODE_CONTEXTS][2];
  int zeromv_mode_cost[ZEROMV_MODE_CONTEXTS][2];
  int refmv_mode_cost[REFMV_MODE_CONTEXTS][2];
  int drl_mode_cost0[DRL_MODE_CONTEXTS][2];
#if CONFIG_EXT_INTER
  int new2mv_mode_cost[2];
#endif  // CONFIG_EXT_INTER
#endif

  unsigned int inter_mode_cost[INTER_MODE_CONTEXTS][INTER_MODES];
#if CONFIG_EXT_INTER
  unsigned int inter_compound_mode_cost[INTER_MODE_CONTEXTS]
                                       [INTER_COMPOUND_MODES];
  unsigned int interintra_mode_cost[BLOCK_SIZE_GROUPS][INTERINTRA_MODES];
#endif  // CONFIG_EXT_INTER
#if CONFIG_OBMC || CONFIG_WARPED_MOTION
  int motvar_cost[BLOCK_SIZES][MOTION_VARIATIONS];
#endif  // CONFIG_OBMC || CONFIG_WARPED_MOTION
  int intra_uv_mode_cost[INTRA_MODES][INTRA_MODES];
  int y_mode_costs[INTRA_MODES][INTRA_MODES][INTRA_MODES];
  int switchable_interp_costs[SWITCHABLE_FILTER_CONTEXTS][SWITCHABLE_FILTERS];
#if CONFIG_EXT_PARTITION_TYPES
  int partition_cost[PARTITION_CONTEXTS][EXT_PARTITION_TYPES];
#else
  int partition_cost[PARTITION_CONTEXTS][PARTITION_TYPES];
#endif
  int palette_y_size_cost[PALETTE_BLOCK_SIZES][PALETTE_SIZES];
  int palette_uv_size_cost[PALETTE_BLOCK_SIZES][PALETTE_SIZES];
  int palette_y_color_cost[PALETTE_MAX_SIZE - 1][PALETTE_COLOR_CONTEXTS]
                                                 [PALETTE_COLORS];
  int palette_uv_color_cost[PALETTE_MAX_SIZE - 1][PALETTE_COLOR_CONTEXTS]
                                                  [PALETTE_COLORS];
  int tx_size_cost[TX_SIZES - 1][TX_SIZE_CONTEXTS][TX_SIZES];
#if CONFIG_EXT_TX
  int inter_tx_type_costs[EXT_TX_SETS_INTER][EXT_TX_SIZES][TX_TYPES];
  int intra_tx_type_costs[EXT_TX_SETS_INTRA][EXT_TX_SIZES][INTRA_MODES]
                                                          [TX_TYPES];
#else
  int intra_tx_type_costs[EXT_TX_SIZES][TX_TYPES][TX_TYPES];
  int inter_tx_type_costs[EXT_TX_SIZES][TX_TYPES];
#endif  // CONFIG_EXT_TX
#if CONFIG_EXT_INTRA
  int intra_filter_cost[INTRA_FILTERS + 1][INTRA_FILTERS];
#endif  // CONFIG_EXT_INTRA

  int multi_arf_allowed;
  int multi_arf_enabled;
  int multi_arf_last_grp_enabled;

  TileDataEnc *tile_data;
  int allocated_tiles;  // Keep track of memory allocated for tiles.

  TOKENEXTRA *tile_tok[MAX_TILE_ROWS][MAX_TILE_COLS];
  unsigned int tok_count[MAX_TILE_ROWS][MAX_TILE_COLS];

  TileBufferEnc tile_buffers[MAX_TILE_ROWS][MAX_TILE_COLS];

#if CONFIG_VP9_TEMPORAL_DENOISING
  VP9_DENOISER denoiser;
#endif

  int resize_pending;
  int resize_state;
  int resize_scale_num;
  int resize_scale_den;
  int resize_avg_qp;
  int resize_buffer_underflow;
  int resize_count;

  // VAR_BASED_PARTITION thresholds
  // 0 - threshold_128x128;
  // 1 - threshold_64x64;
  // 2 - threshold_32x32;
  // 3 - threshold_16x16;
  // 4 - threshold_8x8;
  int64_t vbp_thresholds[5];
  int64_t vbp_threshold_minmax;
  int64_t vbp_threshold_sad;
  BLOCK_SIZE vbp_bsize_min;

  // VARIANCE_AQ segment map refresh
  int vaq_refresh;

  // Multi-threading
  int num_workers;
  VPxWorker *workers;
  struct EncWorkerData *tile_thr_data;
  VP9LfSync lf_row_sync;
#if CONFIG_ENTROPY
  SUBFRAME_STATS subframe_stats;
  // TODO(yaowu): minimize the size of count buffers
  SUBFRAME_STATS wholeframe_stats;
  vp10_coeff_stats branch_ct_buf[COEF_PROBS_BUFS][TX_SIZES][PLANE_TYPES];
#endif  // CONFIG_ENTROPY
#if CONFIG_ANS
  struct BufAnsCoder buf_ans;
#endif
#if CONFIG_EXT_REFS
  int refresh_frame_mask;
  int existing_fb_idx_to_show;
  int is_arf_filter_off;
#endif  // CONFIG_EXT_REFS
} VP10_COMP;

void vp10_initialize_enc(void);

struct VP10_COMP *vp10_create_compressor(VP10EncoderConfig *oxcf,
                                       BufferPool *const pool);
void vp10_remove_compressor(VP10_COMP *cpi);

void vp10_change_config(VP10_COMP *cpi, const VP10EncoderConfig *oxcf);

  // receive a frames worth of data. caller can assume that a copy of this
  // frame is made and not just a copy of the pointer..
int vp10_receive_raw_frame(VP10_COMP *cpi, unsigned int frame_flags,
                          YV12_BUFFER_CONFIG *sd, int64_t time_stamp,
                          int64_t end_time_stamp);

int vp10_get_compressed_data(VP10_COMP *cpi, unsigned int *frame_flags,
                            size_t *size, uint8_t *dest,
                            int64_t *time_stamp, int64_t *time_end, int flush);

int vp10_get_preview_raw_frame(VP10_COMP *cpi, YV12_BUFFER_CONFIG *dest,
                              vp10_ppflags_t *flags);

int vp10_get_last_show_frame(VP10_COMP *cpi, YV12_BUFFER_CONFIG *frame);

int vp10_use_as_reference(VP10_COMP *cpi, int ref_frame_flags);

void vp10_update_reference(VP10_COMP *cpi, int ref_frame_flags);

int vp10_copy_reference_enc(VP10_COMP *cpi, VP9_REFFRAME ref_frame_flag,
                           YV12_BUFFER_CONFIG *sd);

int vp10_set_reference_enc(VP10_COMP *cpi, VP9_REFFRAME ref_frame_flag,
                          YV12_BUFFER_CONFIG *sd);

int vp10_update_entropy(VP10_COMP *cpi, int update);

int vp10_set_active_map(VP10_COMP *cpi, unsigned char *map, int rows, int cols);

int vp10_get_active_map(VP10_COMP *cpi, unsigned char *map, int rows, int cols);

int vp10_set_internal_size(VP10_COMP *cpi,
                          VPX_SCALING horiz_mode, VPX_SCALING vert_mode);

int vp10_set_size_literal(VP10_COMP *cpi, unsigned int width,
                         unsigned int height);

int vp10_get_quantizer(struct VP10_COMP *cpi);

void vp10_full_to_model_counts(vp10_coeff_count_model *model_count,
                               vp10_coeff_count *full_count);

static INLINE int frame_is_kf_gf_arf(const VP10_COMP *cpi) {
  return frame_is_intra_only(&cpi->common) ||
         cpi->refresh_alt_ref_frame ||
         (cpi->refresh_golden_frame && !cpi->rc.is_src_frame_alt_ref);
}

static INLINE int get_ref_frame_map_idx(const VP10_COMP *cpi,
                                        MV_REFERENCE_FRAME ref_frame) {
#if CONFIG_EXT_REFS
  if (ref_frame >= LAST_FRAME && ref_frame <= LAST3_FRAME)
    return cpi->lst_fb_idxes[ref_frame - 1];
#else
  if (ref_frame == LAST_FRAME)
    return cpi->lst_fb_idx;
#endif  // CONFIG_EXT_REFS
  else if (ref_frame == GOLDEN_FRAME)
    return cpi->gld_fb_idx;
#if CONFIG_EXT_REFS
  else if (ref_frame == BWDREF_FRAME)
    return cpi->bwd_fb_idx;
#endif  // CONFIG_EXT_REFS
  else
    return cpi->alt_fb_idx;
}

static INLINE int get_ref_frame_buf_idx(const VP10_COMP *const cpi,
                                        MV_REFERENCE_FRAME ref_frame) {
  const VP10_COMMON *const cm = &cpi->common;
  const int map_idx = get_ref_frame_map_idx(cpi, ref_frame);
  return (map_idx != INVALID_IDX) ? cm->ref_frame_map[map_idx] : INVALID_IDX;
}

static INLINE YV12_BUFFER_CONFIG *get_ref_frame_buffer(
    VP10_COMP *cpi, MV_REFERENCE_FRAME ref_frame) {
  VP10_COMMON *const cm = &cpi->common;
  const int buf_idx = get_ref_frame_buf_idx(cpi, ref_frame);
  return
      buf_idx != INVALID_IDX ? &cm->buffer_pool->frame_bufs[buf_idx].buf : NULL;
}

static INLINE const YV12_BUFFER_CONFIG *get_upsampled_ref(
    VP10_COMP *cpi, const MV_REFERENCE_FRAME ref_frame) {
  // Use up-sampled reference frames.
  const int buf_idx =
      cpi->upsampled_ref_idx[get_ref_frame_map_idx(cpi, ref_frame)];
  return &cpi->upsampled_ref_bufs[buf_idx].buf;
}

#if CONFIG_EXT_REFS
static INLINE int enc_is_ref_frame_buf(VP10_COMP *cpi,
                                       RefCntBuffer *frame_buf) {
  MV_REFERENCE_FRAME ref_frame;
  VP10_COMMON *const cm = &cpi->common;
  for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
    const int buf_idx = get_ref_frame_buf_idx(cpi, ref_frame);
    if (buf_idx == INVALID_IDX) continue;
    if (frame_buf == &cm->buffer_pool->frame_bufs[buf_idx]) break;
  }
  return (ref_frame <= ALTREF_FRAME);
}
#endif  // CONFIG_EXT_REFS

static INLINE unsigned int get_token_alloc(int mb_rows, int mb_cols) {
  // TODO(JBB): double check we can't exceed this token count if we have a
  // 32x32 transform crossing a boundary at a multiple of 16.
  // mb_rows, cols are in units of 16 pixels. We assume 3 planes all at full
  // resolution. We assume up to 1 token per pixel, and then allow
  // a head room of 1 EOSB token per 8x8 block per plane.
  return mb_rows * mb_cols * (16 * 16 + 4) * 3;
}

// Get the allocated token size for a tile. It does the same calculation as in
// the frame token allocation.
static INLINE unsigned int allocated_tokens(TileInfo tile) {
  int tile_mb_rows = (tile.mi_row_end - tile.mi_row_start + 1) >> 1;
  int tile_mb_cols = (tile.mi_col_end - tile.mi_col_start + 1) >> 1;

  return get_token_alloc(tile_mb_rows, tile_mb_cols);
}

void vp10_alloc_compressor_data(VP10_COMP *cpi);

void vp10_scale_references(VP10_COMP *cpi);

void vp10_update_reference_frames(VP10_COMP *cpi);

void vp10_set_high_precision_mv(VP10_COMP *cpi, int allow_high_precision_mv);

YV12_BUFFER_CONFIG *vp10_scale_if_required_fast(VP10_COMMON *cm,
                                               YV12_BUFFER_CONFIG *unscaled,
                                               YV12_BUFFER_CONFIG *scaled);

YV12_BUFFER_CONFIG *vp10_scale_if_required(VP10_COMMON *cm,
                                          YV12_BUFFER_CONFIG *unscaled,
                                          YV12_BUFFER_CONFIG *scaled);

void vp10_apply_encoding_flags(VP10_COMP *cpi, vpx_enc_frame_flags_t flags);

static INLINE int is_altref_enabled(const VP10_COMP *const cpi) {
  return cpi->oxcf.mode != REALTIME && cpi->oxcf.lag_in_frames > 0 &&
         cpi->oxcf.enable_auto_arf;
}

// TODO(zoeliu): To set up cpi->oxcf.enable_auto_brf
#if 0 && CONFIG_EXT_REFS
static INLINE int is_bwdref_enabled(const VP10_COMP *const cpi) {
  // NOTE(zoeliu): The enabling of bi-predictive frames depends on the use of
  //               alt_ref, and now will be off when the alt_ref interval is
  //               not sufficiently large.
  return is_altref_enabled(cpi) && cpi->oxcf.enable_auto_brf;
}
#endif  // CONFIG_EXT_REFS

static INLINE void set_ref_ptrs(VP10_COMMON *cm, MACROBLOCKD *xd,
                                MV_REFERENCE_FRAME ref0,
                                MV_REFERENCE_FRAME ref1) {
  xd->block_refs[0] = &cm->frame_refs[ref0 >= LAST_FRAME ? ref0 - LAST_FRAME
                                                         : 0];
  xd->block_refs[1] = &cm->frame_refs[ref1 >= LAST_FRAME ? ref1 - LAST_FRAME
                                                         : 0];
}

static INLINE int get_chessboard_index(const int frame_index) {
  return frame_index & 0x1;
}

static INLINE int *cond_cost_list(const struct VP10_COMP *cpi, int *cost_list) {
  return cpi->sf.mv.subpel_search_method != SUBPEL_TREE ? cost_list : NULL;
}

void vp10_new_framerate(VP10_COMP *cpi, double framerate);

#define LAYER_IDS_TO_IDX(sl, tl, num_tl) ((sl) * (num_tl) + (tl))

// Update up-sampled reference frame index.
static INLINE void uref_cnt_fb(EncRefCntBuffer *ubufs, int *uidx,
                               int new_uidx) {
  const int ref_index = *uidx;

  if (ref_index >= 0 && ubufs[ref_index].ref_count > 0)
    ubufs[ref_index].ref_count--;

  *uidx = new_uidx;
  ubufs[new_uidx].ref_count++;
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_ENCODER_H_
