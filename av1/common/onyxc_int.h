/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AV1_COMMON_ONYXC_INT_H_
#define AV1_COMMON_ONYXC_INT_H_

#include "./aom_config.h"
#include "./av1_rtcd.h"
#include "aom/internal/aom_codec_internal.h"
#include "aom_util/aom_thread.h"
#include "av1/common/alloccommon.h"
#include "av1/common/entropy.h"
#include "av1/common/entropymode.h"
#include "av1/common/entropymv.h"
#include "av1/common/frame_buffers.h"
#include "av1/common/loopfilter.h"
#include "av1/common/mv.h"
#include "av1/common/quant_common.h"
#if CONFIG_LOOP_RESTORATION
#include "av1/common/restoration.h"
#endif  // CONFIG_LOOP_RESTORATION
#include "av1/common/tile_common.h"
#include "av1/common/odintrin.h"
#if CONFIG_PVQ
#include "av1/common/pvq.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define REF_FRAMES_LOG2 3
#define REF_FRAMES (1 << REF_FRAMES_LOG2)

// 4 scratch frames for the new frames to support a maximum of 4 cores decoding
// in parallel, 3 for scaled references on the encoder.
// TODO(hkuang): Add ondemand frame buffers instead of hardcoding the number
// of framebuffers.
// TODO(jkoleszar): These 3 extra references could probably come from the
// normal reference pool.
#define FRAME_BUFFERS (REF_FRAMES + 7)

#if CONFIG_EXT_REFS
#define FRAME_CONTEXTS_LOG2 3
#else
#define FRAME_CONTEXTS_LOG2 2
#endif

#define FRAME_CONTEXTS (1 << FRAME_CONTEXTS_LOG2)

#define NUM_PING_PONG_BUFFERS 2

typedef enum {
  SINGLE_REFERENCE = 0,
  COMPOUND_REFERENCE = 1,
  REFERENCE_MODE_SELECT = 2,
  REFERENCE_MODES = 3,
} REFERENCE_MODE;

typedef enum {
  RESET_FRAME_CONTEXT_NONE = 0,
  RESET_FRAME_CONTEXT_CURRENT = 1,
  RESET_FRAME_CONTEXT_ALL = 2,
} RESET_FRAME_CONTEXT_MODE;

typedef enum {
  /**
   * Update frame context to values resulting from forward probability
   * updates signaled in the frame header
   */
  REFRESH_FRAME_CONTEXT_FORWARD,
  /**
   * Update frame context to values resulting from backward probability
   * updates based on entropy/counts in the decoded frame
   */
  REFRESH_FRAME_CONTEXT_BACKWARD,
} REFRESH_FRAME_CONTEXT_MODE;

typedef struct {
  int_mv mv[2];
#if CONFIG_REF_MV
  int_mv pred_mv[2];
#endif
  MV_REFERENCE_FRAME ref_frame[2];
} MV_REF;

typedef struct {
  int ref_count;
  MV_REF *mvs;
  int mi_rows;
  int mi_cols;
  aom_codec_frame_buffer_t raw_frame_buffer;
  YV12_BUFFER_CONFIG buf;

  // The Following variables will only be used in frame parallel decode.

  // frame_worker_owner indicates which FrameWorker owns this buffer. NULL means
  // that no FrameWorker owns, or is decoding, this buffer.
  AVxWorker *frame_worker_owner;

  // row and col indicate which position frame has been decoded to in real
  // pixel unit. They are reset to -1 when decoding begins and set to INT_MAX
  // when the frame is fully decoded.
  int row;
  int col;
} RefCntBuffer;

typedef struct BufferPool {
// Protect BufferPool from being accessed by several FrameWorkers at
// the same time during frame parallel decode.
// TODO(hkuang): Try to use atomic variable instead of locking the whole pool.
#if CONFIG_MULTITHREAD
  pthread_mutex_t pool_mutex;
#endif

  // Private data associated with the frame buffer callbacks.
  void *cb_priv;

  aom_get_frame_buffer_cb_fn_t get_fb_cb;
  aom_release_frame_buffer_cb_fn_t release_fb_cb;

  RefCntBuffer frame_bufs[FRAME_BUFFERS];

  // Frame buffers allocated internally by the codec.
  InternalFrameBufferList int_frame_buffers;
} BufferPool;

typedef struct AV1Common {
  struct aom_internal_error_info error;
  aom_color_space_t color_space;
  int color_range;
  int width;
  int height;
  int render_width;
  int render_height;
  int last_width;
  int last_height;

  // TODO(jkoleszar): this implies chroma ss right now, but could vary per
  // plane. Revisit as part of the future change to YV12_BUFFER_CONFIG to
  // support additional planes.
  int subsampling_x;
  int subsampling_y;

#if CONFIG_AOM_HIGHBITDEPTH
  // Marks if we need to use 16bit frame buffers (1: yes, 0: no).
  int use_highbitdepth;
#endif
#if CONFIG_CLPF
  // Two bits are used to signal the strength for all blocks and the
  // valid values are:
  // 0: no filtering
  // 1: strength = 1
  // 2: strength = 2
  // 3: strength = 4
  int clpf_strength_y;
  int clpf_strength_u;
  int clpf_strength_v;

  // If clpf_strength_y is not 0, another two bits are used to signal
  // the filter block size.  The valid values for clfp_size are:
  // 0: no block signalling
  // 1: 32x32
  // 2: 64x64
  // 3: 128x128
  CLPF_BLOCK_SIZE clpf_size;

  // Buffer for storing whether to filter individual blocks.
  int8_t *clpf_blocks;
  int clpf_stride;
#endif

  YV12_BUFFER_CONFIG *frame_to_show;
  RefCntBuffer *prev_frame;

  // TODO(hkuang): Combine this with cur_buf in macroblockd.
  RefCntBuffer *cur_frame;

  int ref_frame_map[REF_FRAMES]; /* maps fb_idx to reference slot */

  // Prepare ref_frame_map for the next frame.
  // Only used in frame parallel decode.
  int next_ref_frame_map[REF_FRAMES];

  // TODO(jkoleszar): could expand active_ref_idx to 4, with 0 as intra, and
  // roll new_fb_idx into it.

  // Each Inter frame can reference INTER_REFS_PER_FRAME buffers
  RefBuffer frame_refs[INTER_REFS_PER_FRAME];

  int new_fb_idx;

  FRAME_TYPE last_frame_type; /* last frame's frame type for motion search.*/
  FRAME_TYPE frame_type;

  int show_frame;
  int last_show_frame;
  int show_existing_frame;
#if CONFIG_EXT_REFS
  // Flag for a frame used as a reference - not written to the bitstream
  int is_reference_frame;
#endif  // CONFIG_EXT_REFS

  // Flag signaling that the frame is encoded using only INTRA modes.
  uint8_t intra_only;
  uint8_t last_intra_only;

  int allow_high_precision_mv;

#if CONFIG_PALETTE
  int allow_screen_content_tools;
#endif  // CONFIG_PALETTE

  // Flag signaling which frame contexts should be reset to default values.
  RESET_FRAME_CONTEXT_MODE reset_frame_context;

  // MBs, mb_rows/cols is in 16-pixel units; mi_rows/cols is in
  // MODE_INFO (8-pixel) units.
  int MBs;
  int mb_rows, mi_rows;
  int mb_cols, mi_cols;
  int mi_stride;

  /* profile settings */
  TX_MODE tx_mode;

  int base_qindex;
  int y_dc_delta_q;
  int uv_dc_delta_q;
  int uv_ac_delta_q;
  int16_t y_dequant[MAX_SEGMENTS][2];
  int16_t uv_dequant[MAX_SEGMENTS][2];

#if CONFIG_AOM_QM
  // Global quant matrix tables
  qm_val_t *giqmatrix[NUM_QM_LEVELS][2][2][TX_SIZES];
  qm_val_t *gqmatrix[NUM_QM_LEVELS][2][2][TX_SIZES];

  // Local quant matrix tables for each frame
  qm_val_t *y_iqmatrix[MAX_SEGMENTS][2][TX_SIZES];
  qm_val_t *uv_iqmatrix[MAX_SEGMENTS][2][TX_SIZES];
  // Encoder
  qm_val_t *y_qmatrix[MAX_SEGMENTS][2][TX_SIZES];
  qm_val_t *uv_qmatrix[MAX_SEGMENTS][2][TX_SIZES];

  int using_qmatrix;
  int min_qmlevel;
  int max_qmlevel;
#endif
#if CONFIG_NEW_QUANT
  dequant_val_type_nuq y_dequant_nuq[MAX_SEGMENTS][QUANT_PROFILES][COEF_BANDS];
  dequant_val_type_nuq uv_dequant_nuq[MAX_SEGMENTS][QUANT_PROFILES][COEF_BANDS];
#endif

  /* We allocate a MODE_INFO struct for each macroblock, together with
     an extra row on top and column on the left to simplify prediction. */
  int mi_alloc_size;
  MODE_INFO *mip; /* Base of allocated array */
  MODE_INFO *mi;  /* Corresponds to upper left visible macroblock */

  // TODO(agrange): Move prev_mi into encoder structure.
  // prev_mip and prev_mi will only be allocated in encoder.
  MODE_INFO *prev_mip; /* MODE_INFO array 'mip' from last decoded frame */
  MODE_INFO *prev_mi;  /* 'mi' from last frame (points into prev_mip) */

  // Separate mi functions between encoder and decoder.
  int (*alloc_mi)(struct AV1Common *cm, int mi_size);
  void (*free_mi)(struct AV1Common *cm);
  void (*setup_mi)(struct AV1Common *cm);

  // Grid of pointers to 8x8 MODE_INFO structs.  Any 8x8 not in the visible
  // area will be NULL.
  MODE_INFO **mi_grid_base;
  MODE_INFO **mi_grid_visible;
  MODE_INFO **prev_mi_grid_base;
  MODE_INFO **prev_mi_grid_visible;

  // Whether to use previous frame's motion vectors for prediction.
  int use_prev_frame_mvs;

  // Persistent mb segment id map used in prediction.
  int seg_map_idx;
  int prev_seg_map_idx;

  uint8_t *seg_map_array[NUM_PING_PONG_BUFFERS];
  uint8_t *last_frame_seg_map;
  uint8_t *current_frame_seg_map;
  int seg_map_alloc_size;

  InterpFilter interp_filter;

  loop_filter_info_n lf_info;
#if CONFIG_LOOP_RESTORATION
  RestorationInfo rst_info;
  RestorationInternal rst_internal;
#endif  // CONFIG_LOOP_RESTORATION

  // Flag signaling how frame contexts should be updated at the end of
  // a frame decode
  REFRESH_FRAME_CONTEXT_MODE refresh_frame_context;

  int ref_frame_sign_bias[TOTAL_REFS_PER_FRAME]; /* Two state 0, 1 */

  struct loopfilter lf;
  struct segmentation seg;

  int frame_parallel_decode;  // frame-based threading.

// Context probabilities for reference frame prediction
#if CONFIG_EXT_REFS
  MV_REFERENCE_FRAME comp_fwd_ref[FWD_REFS];
  MV_REFERENCE_FRAME comp_bwd_ref[BWD_REFS];
#else
  MV_REFERENCE_FRAME comp_fixed_ref;
  MV_REFERENCE_FRAME comp_var_ref[COMP_REFS];
#endif  // CONFIG_EXT_REFS
  REFERENCE_MODE reference_mode;

  FRAME_CONTEXT *fc;              /* this frame entropy */
  FRAME_CONTEXT *frame_contexts;  // FRAME_CONTEXTS
  unsigned int frame_context_idx; /* Context to use/update */
  FRAME_COUNTS counts;

#if CONFIG_ENTROPY
  // The initial probabilities for a frame, before any subframe backward update,
  // and after forward update.
  av1_coeff_probs_model starting_coef_probs[TX_SIZES][PLANE_TYPES];
  // Number of subframe backward updates already done
  uint8_t coef_probs_update_idx;
  // Signal if the backward update is subframe or end-of-frame
  uint8_t partial_prob_update;
  // Frame level flag to turn on/off subframe backward update
  uint8_t do_subframe_update;
#endif  // CONFIG_ENTROPY

  unsigned int current_video_frame;
  BITSTREAM_PROFILE profile;

  // AOM_BITS_8 in profile 0 or 1, AOM_BITS_10 or AOM_BITS_12 in profile 2 or 3.
  aom_bit_depth_t bit_depth;
  aom_bit_depth_t dequant_bit_depth;  // bit_depth of current dequantizer

  int error_resilient_mode;

#if !CONFIG_EXT_TILE
  int log2_tile_cols, log2_tile_rows;
#endif  // !CONFIG_EXT_TILE
  int tile_cols, tile_rows;
  int tile_width, tile_height;  // In MI units

  int byte_alignment;
  int skip_loop_filter;

  // Private data associated with the frame buffer callbacks.
  void *cb_priv;
  aom_get_frame_buffer_cb_fn_t get_fb_cb;
  aom_release_frame_buffer_cb_fn_t release_fb_cb;

  // Handles memory for the codec.
  InternalFrameBufferList int_frame_buffers;

  // External BufferPool passed from outside.
  BufferPool *buffer_pool;

  PARTITION_CONTEXT *above_seg_context;
  ENTROPY_CONTEXT *above_context[MAX_MB_PLANE];
#if CONFIG_VAR_TX
  TXFM_CONTEXT *above_txfm_context;
  TXFM_CONTEXT left_txfm_context[MAX_MIB_SIZE];
#endif
  int above_context_alloc_cols;

  // scratch memory for intraonly/keyframe forward updates from default tables
  // - this is intentionally not placed in FRAME_CONTEXT since it's reset upon
  // each keyframe and not used afterwards
  aom_prob kf_y_prob[INTRA_MODES][INTRA_MODES][INTRA_MODES - 1];
#if CONFIG_DAALA_EC
  aom_cdf_prob kf_y_cdf[INTRA_MODES][INTRA_MODES][INTRA_MODES];
#endif
#if CONFIG_GLOBAL_MOTION
  Global_Motion_Params global_motion[TOTAL_REFS_PER_FRAME];
#endif

  BLOCK_SIZE sb_size;  // Size of the superblock used for this frame
  int mib_size;        // Size of the superblock in units of MI blocks
  int mib_size_log2;   // Log 2 of above.
#if CONFIG_DERING
  int dering_level;
#endif

#if CONFIG_DELTA_Q
  int delta_q_present_flag;
  // Resolution of delta quant
  int delta_q_res;
#endif
#if CONFIG_TILE_GROUPS
  int num_tg;
#endif
} AV1_COMMON;

// TODO(hkuang): Don't need to lock the whole pool after implementing atomic
// frame reference count.
static void lock_buffer_pool(BufferPool *const pool) {
#if CONFIG_MULTITHREAD
  pthread_mutex_lock(&pool->pool_mutex);
#else
  (void)pool;
#endif
}

static void unlock_buffer_pool(BufferPool *const pool) {
#if CONFIG_MULTITHREAD
  pthread_mutex_unlock(&pool->pool_mutex);
#else
  (void)pool;
#endif
}

static INLINE YV12_BUFFER_CONFIG *get_ref_frame(AV1_COMMON *cm, int index) {
  if (index < 0 || index >= REF_FRAMES) return NULL;
  if (cm->ref_frame_map[index] < 0) return NULL;
  assert(cm->ref_frame_map[index] < FRAME_BUFFERS);
  return &cm->buffer_pool->frame_bufs[cm->ref_frame_map[index]].buf;
}

static INLINE YV12_BUFFER_CONFIG *get_frame_new_buffer(
    const AV1_COMMON *const cm) {
  return &cm->buffer_pool->frame_bufs[cm->new_fb_idx].buf;
}

static INLINE int get_free_fb(AV1_COMMON *cm) {
  RefCntBuffer *const frame_bufs = cm->buffer_pool->frame_bufs;
  int i;

  lock_buffer_pool(cm->buffer_pool);
  for (i = 0; i < FRAME_BUFFERS; ++i)
    if (frame_bufs[i].ref_count == 0) break;

  if (i != FRAME_BUFFERS) {
    frame_bufs[i].ref_count = 1;
  } else {
    // Reset i to be INVALID_IDX to indicate no free buffer found.
    i = INVALID_IDX;
  }

  unlock_buffer_pool(cm->buffer_pool);
  return i;
}

static INLINE void ref_cnt_fb(RefCntBuffer *bufs, int *idx, int new_idx) {
  const int ref_index = *idx;

  if (ref_index >= 0 && bufs[ref_index].ref_count > 0)
    bufs[ref_index].ref_count--;

  *idx = new_idx;

  bufs[new_idx].ref_count++;
}

static INLINE int mi_cols_aligned_to_sb(const AV1_COMMON *cm) {
  return ALIGN_POWER_OF_TWO(cm->mi_cols, cm->mib_size_log2);
}

static INLINE int mi_rows_aligned_to_sb(const AV1_COMMON *cm) {
  return ALIGN_POWER_OF_TWO(cm->mi_rows, cm->mib_size_log2);
}

static INLINE int frame_is_intra_only(const AV1_COMMON *const cm) {
  return cm->frame_type == KEY_FRAME || cm->intra_only;
}

static INLINE void av1_init_macroblockd(AV1_COMMON *cm, MACROBLOCKD *xd,
#if CONFIG_PVQ
                                        tran_low_t *pvq_ref_coeff,
#endif
                                        tran_low_t *dqcoeff) {
  int i;
  for (i = 0; i < MAX_MB_PLANE; ++i) {
    xd->plane[i].dqcoeff = dqcoeff;
#if CONFIG_PVQ
    xd->plane[i].pvq_ref_coeff = pvq_ref_coeff;
#endif
    xd->above_context[i] = cm->above_context[i];
    if (xd->plane[i].plane_type == PLANE_TYPE_Y) {
      memcpy(xd->plane[i].seg_dequant, cm->y_dequant, sizeof(cm->y_dequant));
#if CONFIG_AOM_QM
      memcpy(xd->plane[i].seg_iqmatrix, cm->y_iqmatrix, sizeof(cm->y_iqmatrix));
#endif

#if CONFIG_NEW_QUANT
      memcpy(xd->plane[i].seg_dequant_nuq, cm->y_dequant_nuq,
             sizeof(cm->y_dequant_nuq));
#endif
    } else {
      memcpy(xd->plane[i].seg_dequant, cm->uv_dequant, sizeof(cm->uv_dequant));
#if CONFIG_AOM_QM
      memcpy(xd->plane[i].seg_iqmatrix, cm->uv_iqmatrix,
             sizeof(cm->uv_iqmatrix));
#endif
#if CONFIG_NEW_QUANT
      memcpy(xd->plane[i].seg_dequant_nuq, cm->uv_dequant_nuq,
             sizeof(cm->uv_dequant_nuq));
#endif
    }
    xd->fc = cm->fc;
  }
  xd->above_seg_context = cm->above_seg_context;
#if CONFIG_VAR_TX
  xd->above_txfm_context = cm->above_txfm_context;
#endif
  xd->mi_stride = cm->mi_stride;
  xd->error_info = &cm->error;
}

static INLINE void set_skip_context(MACROBLOCKD *xd, int mi_row, int mi_col) {
  const int above_idx = mi_col * 2;
  const int left_idx = (mi_row * 2) & MAX_MIB_MASK_2;
  int i;
  for (i = 0; i < MAX_MB_PLANE; ++i) {
    struct macroblockd_plane *const pd = &xd->plane[i];
    pd->above_context = &xd->above_context[i][above_idx >> pd->subsampling_x];
    pd->left_context = &xd->left_context[i][left_idx >> pd->subsampling_y];
  }
}

static INLINE int calc_mi_size(int len) {
  // len is in mi units.
  return len + MAX_MIB_SIZE;
}

static INLINE void set_plane_n4(MACROBLOCKD *const xd, int bw, int bh, int bwl,
                                int bhl) {
  int i;
  for (i = 0; i < MAX_MB_PLANE; i++) {
    xd->plane[i].n4_w = (bw << 1) >> xd->plane[i].subsampling_x;
    xd->plane[i].n4_h = (bh << 1) >> xd->plane[i].subsampling_y;
    xd->plane[i].n4_wl = bwl - xd->plane[i].subsampling_x;
    xd->plane[i].n4_hl = bhl - xd->plane[i].subsampling_y;
    xd->plane[i].width = xd->plane[i].n4_w * 4;
    xd->plane[i].height = xd->plane[i].n4_h * 4;
  }
}

static INLINE void set_mi_row_col(MACROBLOCKD *xd, const TileInfo *const tile,
                                  int mi_row, int bh, int mi_col, int bw,
                                  int mi_rows, int mi_cols) {
  xd->mb_to_top_edge = -((mi_row * MI_SIZE) * 8);
  xd->mb_to_bottom_edge = ((mi_rows - bh - mi_row) * MI_SIZE) * 8;
  xd->mb_to_left_edge = -((mi_col * MI_SIZE) * 8);
  xd->mb_to_right_edge = ((mi_cols - bw - mi_col) * MI_SIZE) * 8;

  // Are edges available for intra prediction?
  xd->up_available = (mi_row > tile->mi_row_start);
  xd->left_available = (mi_col > tile->mi_col_start);
  if (xd->up_available) {
    xd->above_mi = xd->mi[-xd->mi_stride];
    // above_mi may be NULL in encoder's first pass.
    xd->above_mbmi = xd->above_mi ? &xd->above_mi->mbmi : NULL;
  } else {
    xd->above_mi = NULL;
    xd->above_mbmi = NULL;
  }

  if (xd->left_available) {
    xd->left_mi = xd->mi[-1];
    // left_mi may be NULL in encoder's first pass.
    xd->left_mbmi = xd->left_mi ? &xd->left_mi->mbmi : NULL;
  } else {
    xd->left_mi = NULL;
    xd->left_mbmi = NULL;
  }

  xd->n8_h = bh;
  xd->n8_w = bw;
#if CONFIG_REF_MV
  xd->is_sec_rect = 0;
  if (xd->n8_w < xd->n8_h)
    if (mi_col & (xd->n8_h - 1)) xd->is_sec_rect = 1;

  if (xd->n8_w > xd->n8_h)
    if (mi_row & (xd->n8_w - 1)) xd->is_sec_rect = 1;
#endif
}

static INLINE const aom_prob *get_y_mode_probs(const AV1_COMMON *cm,
                                               const MODE_INFO *mi,
                                               const MODE_INFO *above_mi,
                                               const MODE_INFO *left_mi,
                                               int block) {
  const PREDICTION_MODE above = av1_above_block_mode(mi, above_mi, block);
  const PREDICTION_MODE left = av1_left_block_mode(mi, left_mi, block);
  return cm->kf_y_prob[above][left];
}

#if CONFIG_DAALA_EC
static INLINE aom_cdf_prob *get_y_mode_cdf(AV1_COMMON *cm, const MODE_INFO *mi,
                                           const MODE_INFO *above_mi,
                                           const MODE_INFO *left_mi,
                                           int block) {
  const PREDICTION_MODE above = av1_above_block_mode(mi, above_mi, block);
  const PREDICTION_MODE left = av1_left_block_mode(mi, left_mi, block);
  return cm->kf_y_cdf[above][left];
}
#endif

static INLINE void update_partition_context(MACROBLOCKD *xd, int mi_row,
                                            int mi_col, BLOCK_SIZE subsize,
                                            BLOCK_SIZE bsize) {
  PARTITION_CONTEXT *const above_ctx = xd->above_seg_context + mi_col;
  PARTITION_CONTEXT *const left_ctx =
      xd->left_seg_context + (mi_row & MAX_MIB_MASK);

#if CONFIG_EXT_PARTITION_TYPES
  const int bw = num_8x8_blocks_wide_lookup[bsize];
  const int bh = num_8x8_blocks_high_lookup[bsize];
  memset(above_ctx, partition_context_lookup[subsize].above, bw);
  memset(left_ctx, partition_context_lookup[subsize].left, bh);
#else
  // num_4x4_blocks_wide_lookup[bsize] / 2
  const int bs = num_8x8_blocks_wide_lookup[bsize];

  // update the partition context at the end notes. set partition bits
  // of block sizes larger than the current one to be one, and partition
  // bits of smaller block sizes to be zero.
  memset(above_ctx, partition_context_lookup[subsize].above, bs);
  memset(left_ctx, partition_context_lookup[subsize].left, bs);
#endif  // CONFIG_EXT_PARTITION_TYPES
}

#if CONFIG_EXT_PARTITION_TYPES
static INLINE void update_ext_partition_context(MACROBLOCKD *xd, int mi_row,
                                                int mi_col, BLOCK_SIZE subsize,
                                                BLOCK_SIZE bsize,
                                                PARTITION_TYPE partition) {
  if (bsize >= BLOCK_8X8) {
    const int bsl = b_width_log2_lookup[bsize], hbs = (1 << bsl) / 4;
    BLOCK_SIZE bsize2 = get_subsize(bsize, PARTITION_SPLIT);
    switch (partition) {
      case PARTITION_SPLIT:
        if (bsize != BLOCK_8X8) break;
      case PARTITION_NONE:
      case PARTITION_HORZ:
      case PARTITION_VERT:
        update_partition_context(xd, mi_row, mi_col, subsize, bsize);
        break;
      case PARTITION_HORZ_A:
        update_partition_context(xd, mi_row, mi_col, bsize2, subsize);
        update_partition_context(xd, mi_row + hbs, mi_col, subsize, subsize);
        break;
      case PARTITION_HORZ_B:
        update_partition_context(xd, mi_row, mi_col, subsize, subsize);
        update_partition_context(xd, mi_row + hbs, mi_col, bsize2, subsize);
        break;
      case PARTITION_VERT_A:
        update_partition_context(xd, mi_row, mi_col, bsize2, subsize);
        update_partition_context(xd, mi_row, mi_col + hbs, subsize, subsize);
        break;
      case PARTITION_VERT_B:
        update_partition_context(xd, mi_row, mi_col, subsize, subsize);
        update_partition_context(xd, mi_row, mi_col + hbs, bsize2, subsize);
        break;
      default: assert(0 && "Invalid partition type");
    }
  }
}
#endif  // CONFIG_EXT_PARTITION_TYPES

static INLINE int partition_plane_context(const MACROBLOCKD *xd, int mi_row,
                                          int mi_col, BLOCK_SIZE bsize) {
  const PARTITION_CONTEXT *above_ctx = xd->above_seg_context + mi_col;
  const PARTITION_CONTEXT *left_ctx =
      xd->left_seg_context + (mi_row & MAX_MIB_MASK);
  const int bsl = mi_width_log2_lookup[bsize];
  int above = (*above_ctx >> bsl) & 1, left = (*left_ctx >> bsl) & 1;

  assert(b_width_log2_lookup[bsize] == b_height_log2_lookup[bsize]);
  assert(bsl >= 0);

  return (left * 2 + above) + bsl * PARTITION_PLOFFSET;
}

static INLINE int max_block_wide(const MACROBLOCKD *xd, const BLOCK_SIZE bsize,
                                 const int plane) {
  int max_blocks_wide = block_size_wide[bsize];
  const struct macroblockd_plane *const pd = &xd->plane[plane];

  if (xd->mb_to_right_edge < 0)
    max_blocks_wide += xd->mb_to_right_edge >> (3 + pd->subsampling_x);

  // Scale the width in the transform block unit.
  return max_blocks_wide >> tx_size_wide_log2[0];
}

static INLINE int max_block_high(const MACROBLOCKD *xd, const BLOCK_SIZE bsize,
                                 const int plane) {
  int max_blocks_high = block_size_high[bsize];
  const struct macroblockd_plane *const pd = &xd->plane[plane];

  if (xd->mb_to_bottom_edge < 0)
    max_blocks_high += xd->mb_to_bottom_edge >> (3 + pd->subsampling_y);

  // Scale the width in the transform block unit.
  return max_blocks_high >> tx_size_wide_log2[0];
}

static INLINE void av1_zero_above_context(AV1_COMMON *const cm,
                                          int mi_col_start, int mi_col_end) {
  const int width = mi_col_end - mi_col_start;

  const int offset_y = 2 * mi_col_start;
  const int width_y = 2 * width;
  const int offset_uv = offset_y >> cm->subsampling_x;
  const int width_uv = width_y >> cm->subsampling_x;

  av1_zero_array(cm->above_context[0] + offset_y, width_y);
  av1_zero_array(cm->above_context[1] + offset_uv, width_uv);
  av1_zero_array(cm->above_context[2] + offset_uv, width_uv);

  av1_zero_array(cm->above_seg_context + mi_col_start, width);

#if CONFIG_VAR_TX
  av1_zero_array(cm->above_txfm_context + mi_col_start, width);
#endif  // CONFIG_VAR_TX
}

static INLINE void av1_zero_left_context(MACROBLOCKD *const xd) {
  av1_zero(xd->left_context);
  av1_zero(xd->left_seg_context);
#if CONFIG_VAR_TX
  av1_zero(xd->left_txfm_context_buffer);
#endif
}

#if CONFIG_VAR_TX
static INLINE TX_SIZE get_min_tx_size(const TX_SIZE tx_size) {
  if (tx_size >= TX_SIZES_ALL) assert(0);
  return txsize_sqr_map[tx_size];
}

static INLINE void set_txfm_ctx(TXFM_CONTEXT *txfm_ctx, uint8_t txs, int len) {
  int i;
  for (i = 0; i < len; ++i) txfm_ctx[i] = txs;
}

static INLINE void set_txfm_ctxs(TX_SIZE tx_size, int n8_w, int n8_h,
                                 const MACROBLOCKD *xd) {
  uint8_t bw = tx_size_wide[tx_size];
  uint8_t bh = tx_size_high[tx_size];
  set_txfm_ctx(xd->above_txfm_context, bw, n8_w);
  set_txfm_ctx(xd->left_txfm_context, bh, n8_h);
}

static INLINE void txfm_partition_update(TXFM_CONTEXT *above_ctx,
                                         TXFM_CONTEXT *left_ctx,
                                         TX_SIZE tx_size) {
  BLOCK_SIZE bsize = txsize_to_bsize[tx_size];
  int bh = num_8x8_blocks_high_lookup[bsize];
  int bw = num_8x8_blocks_wide_lookup[bsize];
  uint8_t txw = tx_size_wide[tx_size];
  uint8_t txh = tx_size_high[tx_size];
  int i;
  for (i = 0; i < bh; ++i) left_ctx[i] = txh;
  for (i = 0; i < bw; ++i) above_ctx[i] = txw;
}

static INLINE int txfm_partition_context(TXFM_CONTEXT *above_ctx,
                                         TXFM_CONTEXT *left_ctx,
                                         const BLOCK_SIZE bsize,
                                         const TX_SIZE tx_size) {
  const uint8_t txw = tx_size_wide[tx_size];
  const uint8_t txh = tx_size_high[tx_size];
  const int above = *above_ctx < txw;
  const int left = *left_ctx < txh;
  TX_SIZE max_tx_size = max_txsize_lookup[bsize];
  int category = 15;

  if (max_tx_size == TX_32X32) {
    if (tx_size == TX_32X32)
      category = 0;
    else
      category = 1;
  } else if (max_tx_size == TX_16X16) {
    if (tx_size == TX_16X16)
      category = 2;
    else
      category = 3;
  } else if (max_tx_size == TX_8X8) {
    category = 4;
  }

  if (category == 15) return category;

  return category * 3 + above + left;
}
#endif

static INLINE PARTITION_TYPE get_partition(const AV1_COMMON *const cm,
                                           const int mi_row, const int mi_col,
                                           const BLOCK_SIZE bsize) {
  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) {
    return PARTITION_INVALID;
  } else {
    const int offset = mi_row * cm->mi_stride + mi_col;
    MODE_INFO **mi = cm->mi_grid_visible + offset;
    const MB_MODE_INFO *const mbmi = &mi[0]->mbmi;
    const int bsl = b_width_log2_lookup[bsize];
    const PARTITION_TYPE partition = partition_lookup[bsl][mbmi->sb_type];
#if !CONFIG_EXT_PARTITION_TYPES
    return partition;
#else
    const int hbs = num_8x8_blocks_wide_lookup[bsize] / 2;

    assert(cm->mi_grid_visible[offset] == &cm->mi[offset]);

    if (partition != PARTITION_NONE && bsize > BLOCK_8X8 &&
        mi_row + hbs < cm->mi_rows && mi_col + hbs < cm->mi_cols) {
      const BLOCK_SIZE h = get_subsize(bsize, PARTITION_HORZ_A);
      const BLOCK_SIZE v = get_subsize(bsize, PARTITION_VERT_A);
      const MB_MODE_INFO *const mbmi_right = &mi[hbs]->mbmi;
      const MB_MODE_INFO *const mbmi_below = &mi[hbs * cm->mi_stride]->mbmi;
      if (mbmi->sb_type == h) {
        return mbmi_below->sb_type == h ? PARTITION_HORZ : PARTITION_HORZ_B;
      } else if (mbmi->sb_type == v) {
        return mbmi_right->sb_type == v ? PARTITION_VERT : PARTITION_VERT_B;
      } else if (mbmi_below->sb_type == h) {
        return PARTITION_HORZ_A;
      } else if (mbmi_right->sb_type == v) {
        return PARTITION_VERT_A;
      } else {
        return PARTITION_SPLIT;
      }
    }

    return partition;
#endif  // !CONFIG_EXT_PARTITION_TYPES
  }
}

static INLINE void set_sb_size(AV1_COMMON *const cm, const BLOCK_SIZE sb_size) {
  cm->sb_size = sb_size;
  cm->mib_size = num_8x8_blocks_wide_lookup[cm->sb_size];
  cm->mib_size_log2 = mi_width_log2_lookup[cm->sb_size];
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_ONYXC_INT_H_
