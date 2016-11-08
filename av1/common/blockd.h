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

#ifndef AV1_COMMON_BLOCKD_H_
#define AV1_COMMON_BLOCKD_H_

#include "./aom_config.h"

#include "aom_dsp/aom_dsp_common.h"
#include "aom_ports/mem.h"
#include "aom_scale/yv12config.h"

#include "av1/common/common_data.h"
#include "av1/common/quant_common.h"
#include "av1/common/entropy.h"
#include "av1/common/entropymode.h"
#include "av1/common/mv.h"
#include "av1/common/scale.h"
#include "av1/common/seg_common.h"
#include "av1/common/tile_common.h"
#if CONFIG_PVQ
#include "av1/common/pvq.h"
#include "av1/common/pvq_state.h"
#include "av1/decoder/decint.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_MB_PLANE 3

typedef enum {
  KEY_FRAME = 0,
  INTER_FRAME = 1,
  FRAME_TYPES,
} FRAME_TYPE;

#if CONFIG_EXT_INTERP && SUPPORT_NONINTERPOLATING_FILTERS
#define IsInterpolatingFilter(filter) (av1_is_interpolating_filter(filter))
#else
#define IsInterpolatingFilter(filter) (1)
#endif  // CONFIG_EXT_INTERP && SUPPORT_NONINTERPOLATING_FILTERS

static INLINE int is_inter_mode(PREDICTION_MODE mode) {
#if CONFIG_EXT_INTER
  return mode >= NEARESTMV && mode <= NEW_NEWMV;
#else
  return mode >= NEARESTMV && mode <= NEWMV;
#endif  // CONFIG_EXT_INTER
}

#if CONFIG_PVQ
typedef struct PVQ_INFO {
  int theta[PVQ_MAX_PARTITIONS];
  int max_theta[PVQ_MAX_PARTITIONS];
  int qg[PVQ_MAX_PARTITIONS];
  int k[PVQ_MAX_PARTITIONS];
  od_coeff y[OD_TXSIZE_MAX * OD_TXSIZE_MAX];
  int nb_bands;
  int off[PVQ_MAX_PARTITIONS];
  int size[PVQ_MAX_PARTITIONS];
  int skip_rest;
  int skip_dir;
  int bs;           // log of the block size minus two,
                    // i.e. equivalent to aom's TX_SIZE
  int ac_dc_coded;  // block skip info, indicating whether DC/AC is coded.
                    // bit0: DC coded, bit1 : AC coded (1 means coded)
  tran_low_t dq_dc_residue;
} PVQ_INFO;

typedef struct PVQ_QUEUE {
  PVQ_INFO *buf;  // buffer for pvq info, stored in encoding order
  int curr_pos;   // curr position to write PVQ_INFO
  int buf_len;    // allocated buffer length
  int last_pos;   // last written position of PVQ_INFO in a tile
} PVQ_QUEUE;
#endif

#if CONFIG_EXT_INTER
static INLINE int is_inter_singleref_mode(PREDICTION_MODE mode) {
  return mode >= NEARESTMV && mode <= NEWFROMNEARMV;
}

static INLINE int is_inter_compound_mode(PREDICTION_MODE mode) {
  return mode >= NEAREST_NEARESTMV && mode <= NEW_NEWMV;
}

static INLINE PREDICTION_MODE compound_ref0_mode(PREDICTION_MODE mode) {
  static PREDICTION_MODE lut[MB_MODE_COUNT] = {
    MB_MODE_COUNT,  // DC_PRED            0
    MB_MODE_COUNT,  // V_PRED             1
    MB_MODE_COUNT,  // H_PRED             2
    MB_MODE_COUNT,  // D45_PRED           3
    MB_MODE_COUNT,  // D135_PRED          4
    MB_MODE_COUNT,  // D117_PRED          5
    MB_MODE_COUNT,  // D153_PRED          6
    MB_MODE_COUNT,  // D207_PRED          7
    MB_MODE_COUNT,  // D63_PRED           8
    MB_MODE_COUNT,  // TM_PRED            9
    MB_MODE_COUNT,  // NEARESTMV         10
    MB_MODE_COUNT,  // NEARMV            11
    MB_MODE_COUNT,  // ZEROMV            12
    MB_MODE_COUNT,  // NEWMV             13
    MB_MODE_COUNT,  // NEWFROMNEARMV     14
    NEARESTMV,      // NEAREST_NEARESTMV 15
    NEARESTMV,      // NEAREST_NEARMV    16
    NEARMV,         // NEAR_NEARESTMV    17
    NEARMV,         // NEAR_NEARMV       18
    NEARESTMV,      // NEAREST_NEWMV     19
    NEWMV,          // NEW_NEARESTMV     20
    NEARMV,         // NEAR_NEWMV        21
    NEWMV,          // NEW_NEARMV        22
    ZEROMV,         // ZERO_ZEROMV       23
    NEWMV,          // NEW_NEWMV         24
  };
  assert(is_inter_compound_mode(mode));
  return lut[mode];
}

static INLINE PREDICTION_MODE compound_ref1_mode(PREDICTION_MODE mode) {
  static PREDICTION_MODE lut[MB_MODE_COUNT] = {
    MB_MODE_COUNT,  // DC_PRED            0
    MB_MODE_COUNT,  // V_PRED             1
    MB_MODE_COUNT,  // H_PRED             2
    MB_MODE_COUNT,  // D45_PRED           3
    MB_MODE_COUNT,  // D135_PRED          4
    MB_MODE_COUNT,  // D117_PRED          5
    MB_MODE_COUNT,  // D153_PRED          6
    MB_MODE_COUNT,  // D207_PRED          7
    MB_MODE_COUNT,  // D63_PRED           8
    MB_MODE_COUNT,  // TM_PRED            9
    MB_MODE_COUNT,  // NEARESTMV         10
    MB_MODE_COUNT,  // NEARMV            11
    MB_MODE_COUNT,  // ZEROMV            12
    MB_MODE_COUNT,  // NEWMV             13
    MB_MODE_COUNT,  // NEWFROMNEARMV     14
    NEARESTMV,      // NEAREST_NEARESTMV 15
    NEARMV,         // NEAREST_NEARMV    16
    NEARESTMV,      // NEAR_NEARESTMV    17
    NEARMV,         // NEAR_NEARMV       18
    NEWMV,          // NEAREST_NEWMV     19
    NEARESTMV,      // NEW_NEARESTMV     20
    NEWMV,          // NEAR_NEWMV        21
    NEARMV,         // NEW_NEARMV        22
    ZEROMV,         // ZERO_ZEROMV       23
    NEWMV,          // NEW_NEWMV         24
  };
  assert(is_inter_compound_mode(mode));
  return lut[mode];
}

static INLINE int have_newmv_in_inter_mode(PREDICTION_MODE mode) {
  return (mode == NEWMV || mode == NEWFROMNEARMV || mode == NEW_NEWMV ||
          mode == NEAREST_NEWMV || mode == NEW_NEARESTMV ||
          mode == NEAR_NEWMV || mode == NEW_NEARMV);
}
#else

static INLINE int have_newmv_in_inter_mode(PREDICTION_MODE mode) {
  return (mode == NEWMV);
}
#endif  // CONFIG_EXT_INTER

/* For keyframes, intra block modes are predicted by the (already decoded)
   modes for the Y blocks to the left and above us; for interframes, there
   is a single probability table. */

typedef struct {
  PREDICTION_MODE as_mode;
  int_mv as_mv[2];  // first, second inter predictor motion vectors
#if CONFIG_REF_MV
  int_mv pred_mv[2];
#endif
#if CONFIG_EXT_INTER
  int_mv ref_mv[2];
#endif  // CONFIG_EXT_INTER
} b_mode_info;

typedef int8_t MV_REFERENCE_FRAME;

#if CONFIG_PALETTE
typedef struct {
  // Number of base colors for Y (0) and UV (1)
  uint8_t palette_size[2];
// Value of base colors for Y, U, and V
#if CONFIG_AOM_HIGHBITDEPTH
  uint16_t palette_colors[3 * PALETTE_MAX_SIZE];
#else
  uint8_t palette_colors[3 * PALETTE_MAX_SIZE];
#endif  // CONFIG_AOM_HIGHBITDEPTH
  // Only used by encoder to store the color index of the top left pixel.
  // TODO(huisu): move this to encoder
  uint8_t palette_first_color_idx[2];
} PALETTE_MODE_INFO;
#endif  // CONFIG_PALETTE

#if CONFIG_FILTER_INTRA
typedef struct {
  // 1: an ext intra mode is used; 0: otherwise.
  uint8_t use_filter_intra_mode[PLANE_TYPES];
  FILTER_INTRA_MODE filter_intra_mode[PLANE_TYPES];
} FILTER_INTRA_MODE_INFO;
#endif  // CONFIG_FILTER_INTRA

#if CONFIG_VAR_TX
#define TXB_COEFF_COST_MAP_SIZE (2 * MAX_MIB_SIZE)

// TODO(angiebird): Merge RD_COST and RD_STATS
typedef struct RD_STATS {
  int rate;
  int64_t dist;
  int64_t sse;
  int skip;
#if CONFIG_RD_DEBUG
  int txb_coeff_cost[MAX_MB_PLANE];
  int txb_coeff_cost_map[MAX_MB_PLANE][TXB_COEFF_COST_MAP_SIZE]
                        [TXB_COEFF_COST_MAP_SIZE];
#endif
} RD_STATS;
#endif  // CONFIG_VAR_TX

// This structure now relates to 8x8 block regions.
typedef struct {
  // Common for both INTER and INTRA blocks
  BLOCK_SIZE sb_type;
  PREDICTION_MODE mode;
  TX_SIZE tx_size;
#if CONFIG_VAR_TX
  // TODO(jingning): This effectively assigned a separate entry for each
  // 8x8 block. Apparently it takes much more space than needed.
  TX_SIZE inter_tx_size[MAX_MIB_SIZE][MAX_MIB_SIZE];
  TX_SIZE min_tx_size;
#endif
  int8_t skip;
  int8_t segment_id;
#if CONFIG_SUPERTX
  // Minimum of all segment IDs under the current supertx block.
  int8_t segment_id_supertx;
#endif                      // CONFIG_SUPERTX
  int8_t seg_id_predicted;  // valid only when temporal_update is enabled

  // Only for INTRA blocks
  PREDICTION_MODE uv_mode;
#if CONFIG_PALETTE
  PALETTE_MODE_INFO palette_mode_info;
#endif  // CONFIG_PALETTE

// Only for INTER blocks
#if CONFIG_DUAL_FILTER
  InterpFilter interp_filter[4];
#else
  InterpFilter interp_filter;
#endif
  MV_REFERENCE_FRAME ref_frame[2];
  TX_TYPE tx_type;

#if CONFIG_FILTER_INTRA
  FILTER_INTRA_MODE_INFO filter_intra_mode_info;
#endif  // CONFIG_FILTER_INTRA
#if CONFIG_EXT_INTRA
  int8_t angle_delta[2];
  // To-Do (huisu): this may be replaced by interp_filter
  INTRA_FILTER intra_filter;
#endif  // CONFIG_EXT_INTRA

#if CONFIG_EXT_INTER
  INTERINTRA_MODE interintra_mode;
  // TODO(debargha): Consolidate these flags
  int use_wedge_interintra;
  int interintra_wedge_index;
  int interintra_wedge_sign;
  int use_wedge_interinter;
  int interinter_wedge_index;
  int interinter_wedge_sign;
#endif  // CONFIG_EXT_INTER
  MOTION_MODE motion_mode;
  int_mv mv[2];
  int_mv pred_mv[2];
#if CONFIG_REF_MV
  uint8_t ref_mv_idx;
#endif
#if CONFIG_EXT_PARTITION_TYPES
  PARTITION_TYPE partition;
#endif
#if CONFIG_NEW_QUANT
  int dq_off_index;
  int send_dq_bit;
#endif  // CONFIG_NEW_QUANT
  /* deringing gain *per-superblock* */
  int8_t dering_gain;
#if CONFIG_DELTA_Q
  int current_q_index;
#endif
#if CONFIG_RD_DEBUG
  RD_STATS rd_stats;
  int mi_row;
  int mi_col;
#endif
} MB_MODE_INFO;

typedef struct MODE_INFO {
  MB_MODE_INFO mbmi;
  b_mode_info bmi[4];
} MODE_INFO;

static INLINE PREDICTION_MODE get_y_mode(const MODE_INFO *mi, int block) {
  return mi->mbmi.sb_type < BLOCK_8X8 ? mi->bmi[block].as_mode : mi->mbmi.mode;
}

static INLINE int is_inter_block(const MB_MODE_INFO *mbmi) {
  return mbmi->ref_frame[0] > INTRA_FRAME;
}

static INLINE int has_second_ref(const MB_MODE_INFO *mbmi) {
  return mbmi->ref_frame[1] > INTRA_FRAME;
}

PREDICTION_MODE av1_left_block_mode(const MODE_INFO *cur_mi,
                                    const MODE_INFO *left_mi, int b);

PREDICTION_MODE av1_above_block_mode(const MODE_INFO *cur_mi,
                                     const MODE_INFO *above_mi, int b);

enum mv_precision { MV_PRECISION_Q3, MV_PRECISION_Q4 };

struct buf_2d {
  uint8_t *buf;
  uint8_t *buf0;
  int width;
  int height;
  int stride;
};

typedef struct macroblockd_plane {
  tran_low_t *dqcoeff;
  PLANE_TYPE plane_type;
  int subsampling_x;
  int subsampling_y;
  struct buf_2d dst;
  struct buf_2d pre[2];
  ENTROPY_CONTEXT *above_context;
  ENTROPY_CONTEXT *left_context;
  int16_t seg_dequant[MAX_SEGMENTS][2];
#if CONFIG_NEW_QUANT
  dequant_val_type_nuq seg_dequant_nuq[MAX_SEGMENTS][QUANT_PROFILES]
                                      [COEF_BANDS];
#endif
#if CONFIG_PALETTE
  uint8_t *color_index_map;
#endif  // CONFIG_PALETTE

  // number of 4x4s in current block
  uint16_t n4_w, n4_h;
  // log2 of n4_w, n4_h
  uint8_t n4_wl, n4_hl;
  // block size in pixels
  uint8_t width, height;

#if CONFIG_AOM_QM
  const qm_val_t *seg_iqmatrix[MAX_SEGMENTS][2][TX_SIZES];
#endif
  // encoder
  const int16_t *dequant;
#if CONFIG_NEW_QUANT
  const dequant_val_type_nuq *dequant_val_nuq[QUANT_PROFILES];
#endif  // CONFIG_NEW_QUANT
#if CONFIG_AOM_QM
  const qm_val_t *seg_qmatrix[MAX_SEGMENTS][2][TX_SIZES];
#endif

#if CONFIG_PVQ
  DECLARE_ALIGNED(16, int16_t, pred[MAX_SB_SQUARE]);
  // PVQ: forward transformed predicted image, a reference for PVQ.
  tran_low_t *pvq_ref_coeff;
#endif
} MACROBLOCKD_PLANE;

#define BLOCK_OFFSET(x, i) ((x) + (i)*16)

typedef struct RefBuffer {
  // TODO(dkovalev): idx is not really required and should be removed, now it
  // is used in av1_onyxd_if.c
  int idx;
  YV12_BUFFER_CONFIG *buf;
  struct scale_factors sf;
} RefBuffer;

typedef struct macroblockd {
  struct macroblockd_plane plane[MAX_MB_PLANE];
  uint8_t bmode_blocks_wl;
  uint8_t bmode_blocks_hl;

  FRAME_COUNTS *counts;
  TileInfo tile;

  int mi_stride;

  MODE_INFO **mi;
  MODE_INFO *left_mi;
  MODE_INFO *above_mi;
  MB_MODE_INFO *left_mbmi;
  MB_MODE_INFO *above_mbmi;

  int up_available;
  int left_available;

  const aom_prob (*partition_probs)[PARTITION_TYPES - 1];

  /* Distance of MB away from frame edges */
  int mb_to_left_edge;
  int mb_to_right_edge;
  int mb_to_top_edge;
  int mb_to_bottom_edge;

  FRAME_CONTEXT *fc;

  /* pointers to reference frames */
  const RefBuffer *block_refs[2];

  /* pointer to current frame */
  const YV12_BUFFER_CONFIG *cur_buf;

  ENTROPY_CONTEXT *above_context[MAX_MB_PLANE];
  ENTROPY_CONTEXT left_context[MAX_MB_PLANE][2 * MAX_MIB_SIZE];

  PARTITION_CONTEXT *above_seg_context;
  PARTITION_CONTEXT left_seg_context[MAX_MIB_SIZE];

#if CONFIG_VAR_TX
  TXFM_CONTEXT *above_txfm_context;
  TXFM_CONTEXT *left_txfm_context;
  TXFM_CONTEXT left_txfm_context_buffer[MAX_MIB_SIZE];

  TX_SIZE max_tx_size;
#if CONFIG_SUPERTX
  TX_SIZE supertx_size;
#endif
#endif

  // dimension in the unit of 8x8 block of the current block
  uint8_t n8_w, n8_h;

#if CONFIG_REF_MV
  uint8_t ref_mv_count[MODE_CTX_REF_FRAMES];
  CANDIDATE_MV ref_mv_stack[MODE_CTX_REF_FRAMES][MAX_REF_MV_STACK_SIZE];
  uint8_t is_sec_rect;
#endif

#if CONFIG_PVQ
  daala_dec_ctx daala_dec;
#endif
#if CONFIG_AOM_HIGHBITDEPTH
  /* Bit depth: 8, 10, 12 */
  int bd;
#endif

  int qindex[MAX_SEGMENTS];
  int lossless[MAX_SEGMENTS];
  int corrupted;

  struct aom_internal_error_info *error_info;
#if CONFIG_GLOBAL_MOTION
  Global_Motion_Params *global_motion;
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_DELTA_Q
  int prev_qindex;
  int delta_qindex;
  int current_qindex;
#endif
} MACROBLOCKD;

static INLINE BLOCK_SIZE get_subsize(BLOCK_SIZE bsize,
                                     PARTITION_TYPE partition) {
  if (partition == PARTITION_INVALID)
    return BLOCK_INVALID;
  else
    return subsize_lookup[partition][bsize];
}

static const TX_TYPE intra_mode_to_tx_type_context[INTRA_MODES] = {
  DCT_DCT,    // DC
  ADST_DCT,   // V
  DCT_ADST,   // H
  DCT_DCT,    // D45
  ADST_ADST,  // D135
  ADST_DCT,   // D117
  DCT_ADST,   // D153
  DCT_ADST,   // D207
  ADST_DCT,   // D63
  ADST_ADST,  // TM
};

#if CONFIG_SUPERTX
static INLINE int supertx_enabled(const MB_MODE_INFO *mbmi) {
  return (int)txsize_sqr_map[mbmi->tx_size] >
         AOMMIN(b_width_log2_lookup[mbmi->sb_type],
                b_height_log2_lookup[mbmi->sb_type]);
}
#endif  // CONFIG_SUPERTX

#if CONFIG_EXT_TX
#define ALLOW_INTRA_EXT_TX 1

static const int num_ext_tx_set_inter[EXT_TX_SETS_INTER] = { 1, 16, 12, 2 };
static const int num_ext_tx_set_intra[EXT_TX_SETS_INTRA] = { 1, 7, 5 };

static INLINE int get_ext_tx_set(TX_SIZE tx_size, BLOCK_SIZE bs, int is_inter) {
  tx_size = txsize_sqr_map[tx_size];
  if (tx_size > TX_32X32 || bs < BLOCK_8X8) return 0;
  if (tx_size == TX_32X32) return is_inter ? 3 : 0;
  return (tx_size == TX_16X16 ? 2 : 1);
}

static const int use_intra_ext_tx_for_txsize[EXT_TX_SETS_INTRA][TX_SIZES] = {
  { 0, 0, 0, 0 },  // unused
  { 1, 1, 0, 0 },
  { 0, 0, 1, 0 },
};

static const int use_inter_ext_tx_for_txsize[EXT_TX_SETS_INTER][TX_SIZES] = {
  { 0, 0, 0, 0 },  // unused
  { 1, 1, 0, 0 },
  { 0, 0, 1, 0 },
  { 0, 0, 0, 1 },
};

// Transform types used in each intra set
static const int ext_tx_used_intra[EXT_TX_SETS_INTRA][TX_TYPES] = {
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
};

// Transform types used in each inter set
static const int ext_tx_used_inter[EXT_TX_SETS_INTER][TX_TYPES] = {
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 },
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
};

// 1D Transforms used in inter set, this needs to be changed if
// ext_tx_used_inter is changed
static const int ext_tx_used_inter_1D[EXT_TX_SETS_INTER][TX_TYPES_1D] = {
  { 1, 0, 0, 0 }, { 1, 1, 1, 1 }, { 1, 1, 1, 1 }, { 1, 0, 0, 1 },
};

static INLINE int get_ext_tx_types(TX_SIZE tx_size, BLOCK_SIZE bs,
                                   int is_inter) {
  const int set = get_ext_tx_set(tx_size, bs, is_inter);
  return is_inter ? num_ext_tx_set_inter[set] : num_ext_tx_set_intra[set];
}

#if CONFIG_RECT_TX
static INLINE int is_rect_tx_allowed_bsize(BLOCK_SIZE bsize) {
  static const char LUT[BLOCK_SIZES] = {
    0,  // BLOCK_4X4
    1,  // BLOCK_4X8
    1,  // BLOCK_8X4
    0,  // BLOCK_8X8
    1,  // BLOCK_8X16
    1,  // BLOCK_16X8
    0,  // BLOCK_16X16
    1,  // BLOCK_16X32
    1,  // BLOCK_32X16
    0,  // BLOCK_32X32
    0,  // BLOCK_32X64
    0,  // BLOCK_64X32
    0,  // BLOCK_64X64
#if CONFIG_EXT_PARTITION
    0,  // BLOCK_64X128
    0,  // BLOCK_128X64
    0,  // BLOCK_128X128
#endif  // CONFIG_EXT_PARTITION
  };

  return LUT[bsize];
}

static INLINE int is_rect_tx_allowed(const MACROBLOCKD *xd,
                                     const MB_MODE_INFO *mbmi) {
  return is_inter_block(mbmi) && is_rect_tx_allowed_bsize(mbmi->sb_type) &&
         !xd->lossless[mbmi->segment_id];
}

static INLINE int is_rect_tx(TX_SIZE tx_size) { return tx_size >= TX_SIZES; }
#endif  // CONFIG_RECT_TX
#endif  // CONFIG_EXT_TX

static INLINE TX_SIZE tx_size_from_tx_mode(BLOCK_SIZE bsize, TX_MODE tx_mode,
                                           int is_inter) {
  const TX_SIZE largest_tx_size = tx_mode_to_biggest_tx_size[tx_mode];
  const TX_SIZE max_tx_size = max_txsize_lookup[bsize];

#if CONFIG_EXT_TX && CONFIG_RECT_TX
  if (!is_inter) {
    return AOMMIN(max_tx_size, largest_tx_size);
  } else {
    const TX_SIZE max_rect_tx_size = max_txsize_rect_lookup[bsize];
    if (txsize_sqr_up_map[max_rect_tx_size] <= largest_tx_size) {
      return max_rect_tx_size;
    } else {
      return largest_tx_size;
    }
  }
#else
  (void)is_inter;
  return AOMMIN(max_tx_size, largest_tx_size);
#endif  // CONFIG_EXT_TX && CONFIG_RECT_TX
}

#if CONFIG_FILTER_INTRA
static const TX_TYPE filter_intra_mode_to_tx_type_lookup[FILTER_INTRA_MODES] = {
  DCT_DCT,    // FILTER_DC
  ADST_DCT,   // FILTER_V
  DCT_ADST,   // FILTER_H
  DCT_DCT,    // FILTER_D45
  ADST_ADST,  // FILTER_D135
  ADST_DCT,   // FILTER_D117
  DCT_ADST,   // FILTER_D153
  DCT_ADST,   // FILTER_D207
  ADST_DCT,   // FILTER_D63
  ADST_ADST,  // FILTER_TM
};
#endif  // CONFIG_FILTER_INTRA

#if CONFIG_EXT_INTRA
#define ANGLE_STEP 3
#define MAX_ANGLE_DELTAS 3
extern const int16_t dr_intra_derivative[90];
static const uint8_t mode_to_angle_map[INTRA_MODES] = {
  0, 90, 180, 45, 135, 111, 157, 203, 67, 0,
};

// Returns whether filter selection is needed for a given
// intra prediction angle.
int av1_is_intra_filter_switchable(int angle);
#endif  // CONFIG_EXT_INTRA

#if CONFIG_EXT_TILE
#define FIXED_TX_TYPE 1
#else
#define FIXED_TX_TYPE 0
#endif

static INLINE TX_TYPE get_default_tx_type(PLANE_TYPE plane_type,
                                          const MACROBLOCKD *xd, int block_idx,
                                          TX_SIZE tx_size) {
  const MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;

  if (is_inter_block(mbmi) || plane_type != PLANE_TYPE_Y ||
      xd->lossless[mbmi->segment_id] || tx_size >= TX_32X32)
    return DCT_DCT;

  return intra_mode_to_tx_type_context[plane_type == PLANE_TYPE_Y
                                           ? get_y_mode(xd->mi[0], block_idx)
                                           : mbmi->uv_mode];
}

static INLINE TX_TYPE get_tx_type(PLANE_TYPE plane_type, const MACROBLOCKD *xd,
                                  int block_idx, TX_SIZE tx_size) {
  const MODE_INFO *const mi = xd->mi[0];
  const MB_MODE_INFO *const mbmi = &mi->mbmi;

  if (FIXED_TX_TYPE)
    return get_default_tx_type(plane_type, xd, block_idx, tx_size);

#if CONFIG_EXT_INTRA || CONFIG_FILTER_INTRA
  if (!is_inter_block(mbmi)) {
#if CONFIG_FILTER_INTRA
    const int use_filter_intra_mode_info =
        mbmi->filter_intra_mode_info.use_filter_intra_mode[plane_type];
    const FILTER_INTRA_MODE filter_intra_mode =
        mbmi->filter_intra_mode_info.filter_intra_mode[plane_type];
#endif  // CONFIG_FILTER_INTRA
#if CONFIG_EXT_INTRA
    const PREDICTION_MODE mode = (plane_type == PLANE_TYPE_Y)
                                     ? get_y_mode(mi, block_idx)
                                     : mbmi->uv_mode;
#endif  // CONFIG_EXT_INTRA

    if (xd->lossless[mbmi->segment_id] || tx_size >= TX_32X32) return DCT_DCT;

#if CONFIG_EXT_TX
#if ALLOW_INTRA_EXT_TX
    if (mbmi->sb_type >= BLOCK_8X8 && plane_type == PLANE_TYPE_Y)
      return mbmi->tx_type;
#endif  // ALLOW_INTRA_EXT_TX
#endif  // CONFIG_EXT_TX

#if CONFIG_FILTER_INTRA
    if (use_filter_intra_mode_info)
      return filter_intra_mode_to_tx_type_lookup[filter_intra_mode];
#endif  // CONFIG_FILTER_INTRA
#if CONFIG_EXT_INTRA
    if (mode == DC_PRED) {
      return DCT_DCT;
    } else if (mode == TM_PRED) {
      return ADST_ADST;
    } else {
      int angle = mode_to_angle_map[mode];
      if (mbmi->sb_type >= BLOCK_8X8)
        angle += mbmi->angle_delta[plane_type] * ANGLE_STEP;
      assert(angle > 0 && angle < 270);
      if (angle == 135)
        return ADST_ADST;
      else if (angle < 45 || angle > 225)
        return DCT_DCT;
      else if (angle < 135)
        return ADST_DCT;
      else
        return DCT_ADST;
    }
#endif  // CONFIG_EXT_INTRA
  }
#endif  // CONFIG_EXT_INTRA || CONFIG_FILTER_INTRA

#if CONFIG_EXT_TX
#if EXT_TX_SIZES == 4
  if (xd->lossless[mbmi->segment_id] || txsize_sqr_map[tx_size] > TX_32X32 ||
      (txsize_sqr_map[tx_size] >= TX_32X32 && !is_inter_block(mbmi)))
#else
  if (xd->lossless[mbmi->segment_id] || txsize_sqr_map[tx_size] >= TX_32X32)
#endif
    return DCT_DCT;
  if (mbmi->sb_type >= BLOCK_8X8) {
    if (plane_type == PLANE_TYPE_Y) {
#if !ALLOW_INTRA_EXT_TX
      if (is_inter_block(mbmi))
#endif  // ALLOW_INTRA_EXT_TX
        return mbmi->tx_type;
    }
    if (is_inter_block(mbmi))
      // UV Inter only
      return (mbmi->tx_type == IDTX && txsize_sqr_map[tx_size] == TX_32X32)
                 ? DCT_DCT
                 : mbmi->tx_type;
  }

  // Sub8x8-Inter/Intra OR UV-Intra
  if (is_inter_block(mbmi))  // Sub8x8-Inter
    return DCT_DCT;
  else  // Sub8x8 Intra OR UV-Intra
    return intra_mode_to_tx_type_context[plane_type == PLANE_TYPE_Y
                                             ? get_y_mode(mi, block_idx)
                                             : mbmi->uv_mode];
#else   // CONFIG_EXT_TX
  (void)block_idx;
  if (plane_type != PLANE_TYPE_Y || xd->lossless[mbmi->segment_id] ||
      txsize_sqr_map[tx_size] >= TX_32X32)
    return DCT_DCT;
  return mbmi->tx_type;
#endif  // CONFIG_EXT_TX
}

void av1_setup_block_planes(MACROBLOCKD *xd, int ss_x, int ss_y);

static INLINE int tx_size_to_depth(const TX_SIZE tx_size) {
  return (int)(tx_size - TX_4X4);
}

static INLINE TX_SIZE depth_to_tx_size(const int depth) {
  return (TX_SIZE)(depth + TX_4X4);
}

static INLINE TX_SIZE get_uv_tx_size(const MB_MODE_INFO *mbmi,
                                     const struct macroblockd_plane *pd) {
  TX_SIZE uv_txsize;
#if CONFIG_SUPERTX
  if (supertx_enabled(mbmi))
    return uvsupertx_size_lookup[txsize_sqr_map[mbmi->tx_size]]
                                [pd->subsampling_x][pd->subsampling_y];
#endif  // CONFIG_SUPERTX
  uv_txsize = uv_txsize_lookup[mbmi->sb_type][mbmi->tx_size][pd->subsampling_x]
                              [pd->subsampling_y];
  assert(uv_txsize != TX_INVALID);
  return uv_txsize;
}

static INLINE BLOCK_SIZE
get_plane_block_size(BLOCK_SIZE bsize, const struct macroblockd_plane *pd) {
  return ss_size_lookup[bsize][pd->subsampling_x][pd->subsampling_y];
}

static INLINE void reset_skip_context(MACROBLOCKD *xd, BLOCK_SIZE bsize) {
  int i;
  for (i = 0; i < MAX_MB_PLANE; i++) {
    struct macroblockd_plane *const pd = &xd->plane[i];
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize, pd);
    memset(pd->above_context, 0,
           sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_wide_lookup[plane_bsize]);
    memset(pd->left_context, 0,
           sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_high_lookup[plane_bsize]);
  }
}

typedef void (*foreach_transformed_block_visitor)(int plane, int block,
                                                  int blk_row, int blk_col,
                                                  BLOCK_SIZE plane_bsize,
                                                  TX_SIZE tx_size, void *arg);

void av1_foreach_transformed_block_in_plane(
    const MACROBLOCKD *const xd, BLOCK_SIZE bsize, int plane,
    foreach_transformed_block_visitor visit, void *arg);

void av1_foreach_transformed_block(const MACROBLOCKD *const xd,
                                   BLOCK_SIZE bsize,
                                   foreach_transformed_block_visitor visit,
                                   void *arg);

void av1_set_contexts(const MACROBLOCKD *xd, struct macroblockd_plane *pd,
                      TX_SIZE tx_size, int has_eob, int aoff, int loff);

#if CONFIG_EXT_INTER
static INLINE int is_interintra_allowed_bsize(const BLOCK_SIZE bsize) {
  // TODO(debargha): Should this be bsize < BLOCK_LARGEST?
  return (bsize >= BLOCK_8X8) && (bsize < BLOCK_64X64);
}

static INLINE int is_interintra_allowed_mode(const PREDICTION_MODE mode) {
  return (mode >= NEARESTMV) && (mode <= NEWMV);
}

static INLINE int is_interintra_allowed_ref(const MV_REFERENCE_FRAME rf[2]) {
  return (rf[0] > INTRA_FRAME) && (rf[1] <= INTRA_FRAME);
}

static INLINE int is_interintra_allowed(const MB_MODE_INFO *mbmi) {
  return is_interintra_allowed_bsize(mbmi->sb_type) &&
         is_interintra_allowed_mode(mbmi->mode) &&
         is_interintra_allowed_ref(mbmi->ref_frame);
}

static INLINE int is_interintra_allowed_bsize_group(const int group) {
  int i;
  for (i = 0; i < BLOCK_SIZES; i++) {
    if (size_group_lookup[i] == group &&
        is_interintra_allowed_bsize((BLOCK_SIZE)i)) {
      return 1;
    }
  }
  return 0;
}

static INLINE int is_interintra_pred(const MB_MODE_INFO *mbmi) {
  return (mbmi->ref_frame[1] == INTRA_FRAME) && is_interintra_allowed(mbmi);
}
#endif  // CONFIG_EXT_INTER

#if CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
static INLINE int is_motion_variation_allowed_bsize(BLOCK_SIZE bsize) {
  return (bsize >= BLOCK_8X8);
}

static INLINE int is_motion_variation_allowed(const MB_MODE_INFO *mbmi) {
#if CONFIG_EXT_INTER
  return is_motion_variation_allowed_bsize(mbmi->sb_type) &&
         mbmi->ref_frame[1] != INTRA_FRAME;
#else
  return is_motion_variation_allowed_bsize(mbmi->sb_type);
#endif  // CONFIG_EXT_INTER
}

#if CONFIG_MOTION_VAR
static INLINE int is_neighbor_overlappable(const MB_MODE_INFO *mbmi) {
  return (is_inter_block(mbmi));
}
#endif  // CONFIG_MOTION_VAR
#endif  // CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_BLOCKD_H_
