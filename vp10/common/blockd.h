/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP10_COMMON_BLOCKD_H_
#define VP10_COMMON_BLOCKD_H_

#include "./vpx_config.h"

#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_ports/mem.h"
#include "vpx_scale/yv12config.h"

#include "vp10/common/common_data.h"
#include "vp10/common/quant_common.h"
#include "vp10/common/entropy.h"
#include "vp10/common/entropymode.h"
#include "vp10/common/mv.h"
#include "vp10/common/scale.h"
#include "vp10/common/seg_common.h"
#include "vp10/common/tile_common.h"

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_NEW_QUANT
#define QUANT_PROFILES 3

static INLINE int get_dq_profile_from_ctx(int q_ctx) {
  return VPXMIN(q_ctx, QUANT_PROFILES - 1);
}
#endif  // CONFIG_NEW_QUANT

#define MAX_MB_PLANE 3

typedef enum {
  KEY_FRAME = 0,
  INTER_FRAME = 1,
  FRAME_TYPES,
} FRAME_TYPE;

#if CONFIG_EXT_INTERP && SUPPORT_NONINTERPOLATING_FILTERS
#define IsInterpolatingFilter(filter)  (vp10_is_interpolating_filter(filter))
#else
#define IsInterpolatingFilter(filter)  (1)
#endif  // CONFIG_EXT_INTERP && SUPPORT_NONINTERPOLATING_FILTERS

static INLINE int is_inter_mode(PREDICTION_MODE mode) {
#if CONFIG_EXT_INTER
  return mode >= NEARESTMV && mode <= NEW_NEWMV;
#else
  return mode >= NEARESTMV && mode <= NEWMV;
#endif  // CONFIG_EXT_INTER
}

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
  return (mode == NEWMV || mode == NEWFROMNEARMV ||
          mode == NEW_NEWMV ||
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
  int_mv pred_mv_s8[2];
#endif
#if CONFIG_EXT_INTER
  int_mv ref_mv[2];
#endif  // CONFIG_EXT_INTER
} b_mode_info;

typedef int8_t MV_REFERENCE_FRAME;

typedef struct {
  // Number of base colors for Y (0) and UV (1)
  uint8_t palette_size[2];
  // Value of base colors for Y, U, and V
#if CONFIG_VP9_HIGHBITDEPTH
  uint16_t palette_colors[3 * PALETTE_MAX_SIZE];
#else
  uint8_t palette_colors[3 * PALETTE_MAX_SIZE];
#endif  // CONFIG_VP9_HIGHBITDEPTH
  // Only used by encoder to store the color index of the top left pixel.
  // TODO(huisu): move this to encoder
  uint8_t palette_first_color_idx[2];
} PALETTE_MODE_INFO;

#if CONFIG_EXT_INTRA
typedef struct {
  // 1: an ext intra mode is used; 0: otherwise.
  uint8_t use_ext_intra_mode[PLANE_TYPES];
  EXT_INTRA_MODE ext_intra_mode[PLANE_TYPES];
} EXT_INTRA_MODE_INFO;
#endif  // CONFIG_EXT_INTRA

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
#endif
  int8_t skip;
  int8_t has_no_coeffs;
  int8_t segment_id;
#if CONFIG_SUPERTX
  // Minimum of all segment IDs under the current supertx block.
  int8_t segment_id_supertx;
#endif  // CONFIG_SUPERTX
  int8_t seg_id_predicted;  // valid only when temporal_update is enabled

  // Only for INTRA blocks
  PREDICTION_MODE uv_mode;
  PALETTE_MODE_INFO palette_mode_info;

  // Only for INTER blocks
#if CONFIG_DUAL_FILTER
  INTERP_FILTER interp_filter[4];
#else
  INTERP_FILTER interp_filter;
#endif
  MV_REFERENCE_FRAME ref_frame[2];
  TX_TYPE tx_type;

#if CONFIG_EXT_INTRA
  EXT_INTRA_MODE_INFO ext_intra_mode_info;
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
  MOTION_VARIATION motion_variation;
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
} MB_MODE_INFO;

typedef struct MODE_INFO {
  MB_MODE_INFO mbmi;
  b_mode_info bmi[4];
} MODE_INFO;

static INLINE PREDICTION_MODE get_y_mode(const MODE_INFO *mi, int block) {
  return mi->mbmi.sb_type < BLOCK_8X8 ? mi->bmi[block].as_mode
                                      : mi->mbmi.mode;
}

static INLINE int is_inter_block(const MB_MODE_INFO *mbmi) {
  return mbmi->ref_frame[0] > INTRA_FRAME;
}

static INLINE int has_second_ref(const MB_MODE_INFO *mbmi) {
  return mbmi->ref_frame[1] > INTRA_FRAME;
}

PREDICTION_MODE vp10_left_block_mode(const MODE_INFO *cur_mi,
                                    const MODE_INFO *left_mi, int b);

PREDICTION_MODE vp10_above_block_mode(const MODE_INFO *cur_mi,
                                     const MODE_INFO *above_mi, int b);

enum mv_precision {
  MV_PRECISION_Q3,
  MV_PRECISION_Q4
};

struct buf_2d {
  uint8_t *buf;
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
  dequant_val_type_nuq
    seg_dequant_nuq[MAX_SEGMENTS][QUANT_PROFILES][COEF_BANDS];
#endif
  uint8_t *color_index_map;

  // number of 4x4s in current block
  uint16_t n4_w, n4_h;
  // log2 of n4_w, n4_h
  uint8_t n4_wl, n4_hl;

  // encoder
  const int16_t *dequant;
#if CONFIG_NEW_QUANT
  const dequant_val_type_nuq* dequant_val_nuq[QUANT_PROFILES];
#endif  // CONFIG_NEW_QUANT
} MACROBLOCKD_PLANE;

#define BLOCK_OFFSET(x, i) ((x) + (i) * 16)

typedef struct RefBuffer {
  // TODO(dkovalev): idx is not really required and should be removed, now it
  // is used in vp10_onyxd_if.c
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

  const vpx_prob (*partition_probs)[PARTITION_TYPES - 1];

  /* Distance of MB away from frame edges */
  int mb_to_left_edge;
  int mb_to_right_edge;
  int mb_to_top_edge;
  int mb_to_bottom_edge;

  FRAME_CONTEXT *fc;

  /* pointers to reference frames */
  RefBuffer *block_refs[2];

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

#if CONFIG_VP9_HIGHBITDEPTH
  /* Bit depth: 8, 10, 12 */
  int bd;
#endif

  int lossless[MAX_SEGMENTS];
  int corrupted;

  struct vpx_internal_error_info *error_info;
} MACROBLOCKD;

static INLINE BLOCK_SIZE get_subsize(BLOCK_SIZE bsize,
                                     PARTITION_TYPE partition) {
  if (partition == PARTITION_INVALID)
    return PARTITION_INVALID;
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
  return (int)mbmi->tx_size >
      VPXMIN(b_width_log2_lookup[mbmi->sb_type],
             b_height_log2_lookup[mbmi->sb_type]);
}
#endif  // CONFIG_SUPERTX

#if CONFIG_EXT_TX
#define ALLOW_INTRA_EXT_TX          1
// whether masked transforms are used for 32X32
#define USE_MSKTX_FOR_32X32         0
#define USE_REDUCED_TXSET_FOR_16X16 1

static const int num_ext_tx_set_inter[EXT_TX_SETS_INTER] = {
  1, 16, 12, 2
};
static const int num_ext_tx_set_intra[EXT_TX_SETS_INTRA] = {
  1, 7, 5
};

#if EXT_TX_SIZES == 4
static INLINE int get_ext_tx_set(TX_SIZE tx_size, BLOCK_SIZE bs,
                                 int is_inter) {
  if (tx_size > TX_32X32 || bs < BLOCK_8X8) return 0;
#if USE_REDUCED_TXSET_FOR_16X16
  if (tx_size == TX_32X32)
    return is_inter ? 3 - USE_MSKTX_FOR_32X32 : 0;
  return (tx_size == TX_16X16 ? 2 : 1);
#else
  if (tx_size == TX_32X32)
    return is_inter ? 3 - 2 * USE_MSKTX_FOR_32X32 : 0;
  return (tx_size == TX_16X16 && !is_inter ? 2 : 1);
#endif  // USE_REDUCED_TXSET_FOR_16X16
}

static const int use_intra_ext_tx_for_txsize[EXT_TX_SETS_INTRA][TX_SIZES] = {
  { 0, 0, 0, 0, },  // unused
  { 1, 1, 0, 0, },
  { 0, 0, 1, 0, },
};

static const int use_inter_ext_tx_for_txsize[EXT_TX_SETS_INTER][TX_SIZES] = {
  { 0, 0, 0, 0, },  // unused
  { 1, 1, (!USE_REDUCED_TXSET_FOR_16X16), USE_MSKTX_FOR_32X32, },
  { 0, 0, USE_REDUCED_TXSET_FOR_16X16, 0, },
  { 0, 0, 0, (!USE_MSKTX_FOR_32X32), },
};

#else  // EXT_TX_SIZES == 4

static INLINE int get_ext_tx_set(TX_SIZE tx_size, BLOCK_SIZE bs,
                                 int is_inter) {
  (void) is_inter;
  if (tx_size > TX_32X32 || bs < BLOCK_8X8) return 0;
  if (tx_size == TX_32X32) return 0;
#if USE_REDUCED_TXSET_FOR_16X16
  return (tx_size == TX_16X16 ? 2 : 1);
#else
  return (tx_size == TX_16X16 && !is_inter ? 2 : 1);
#endif  // USE_REDUCED_TXSET_FOR_16X16
}

static const int use_intra_ext_tx_for_txsize[EXT_TX_SETS_INTRA][TX_SIZES] = {
  { 0, 0, 0, 0, },  // unused
  { 1, 1, 0, 0, },
  { 0, 0, 1, 0, },
};

static const int use_inter_ext_tx_for_txsize[EXT_TX_SETS_INTER][TX_SIZES] = {
  { 0, 0, 0, 0, },  // unused
  { 1, 1, (!USE_REDUCED_TXSET_FOR_16X16), 0, },
  { 0, 0, USE_REDUCED_TXSET_FOR_16X16, 0, },
  { 0, 0, 0, 1, },
};
#endif  // EXT_TX_SIZES == 4

// Transform types used in each intra set
static const int ext_tx_used_intra[EXT_TX_SETS_INTRA][TX_TYPES] = {
  {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0},
  {1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
};

// Transform types used in each inter set
static const int ext_tx_used_inter[EXT_TX_SETS_INTER][TX_TYPES] = {
  {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0},
  {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
};

// 1D Transforms used in inter set, this needs to be changed if
// ext_tx_used_inter is changed
static const int ext_tx_used_inter_1D[EXT_TX_SETS_INTER][TX_TYPES_1D] = {
  {1, 0, 0, 0},
  {1, 1, 1, 1},
  {1, 1, 1, 1},
  {1, 0, 0, 1},
};

static INLINE int get_ext_tx_types(TX_SIZE tx_size, BLOCK_SIZE bs,
                                   int is_inter) {
  const int set = get_ext_tx_set(tx_size, bs, is_inter);
  return is_inter ? num_ext_tx_set_inter[set] : num_ext_tx_set_intra[set];
}
#endif  // CONFIG_EXT_TX

#if CONFIG_EXT_INTRA
#define ALLOW_FILTER_INTRA_MODES 1
#define ANGLE_STEP 3
#define MAX_ANGLE_DELTAS 3

extern const int16_t dr_intra_derivative[270][2];

static const uint8_t mode_to_angle_map[INTRA_MODES] = {
    0, 90, 180, 45, 135, 111, 157, 203, 67, 0,
};

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

int vp10_is_intra_filter_switchable(int angle);
#endif  // CONFIG_EXT_INTRA

#if CONFIG_EXT_TILE
#define FIXED_TX_TYPE 1
#else
#define FIXED_TX_TYPE 0
#endif

static INLINE TX_TYPE get_default_tx_type(PLANE_TYPE plane_type,
                                          const MACROBLOCKD *xd,
                                          int block_idx, TX_SIZE tx_size) {
  const MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;

  if (is_inter_block(mbmi) || plane_type != PLANE_TYPE_Y ||
      xd->lossless[mbmi->segment_id] || tx_size >= TX_32X32)
    return DCT_DCT;

  return intra_mode_to_tx_type_context[plane_type == PLANE_TYPE_Y ?
      get_y_mode(xd->mi[0], block_idx) : mbmi->uv_mode];
}

static INLINE TX_TYPE get_tx_type(PLANE_TYPE plane_type,
                                  const MACROBLOCKD *xd,
                                  int block_idx, TX_SIZE tx_size) {
  const MODE_INFO *const mi = xd->mi[0];
  const MB_MODE_INFO *const mbmi = &mi->mbmi;

  if (FIXED_TX_TYPE)
    return get_default_tx_type(plane_type, xd, block_idx, tx_size);

#if CONFIG_EXT_INTRA
  if (!is_inter_block(mbmi)) {
    const int use_ext_intra_mode_info =
        mbmi->ext_intra_mode_info.use_ext_intra_mode[plane_type];
    const EXT_INTRA_MODE ext_intra_mode =
        mbmi->ext_intra_mode_info.ext_intra_mode[plane_type];
    const PREDICTION_MODE mode = (plane_type == PLANE_TYPE_Y) ?
        get_y_mode(mi, block_idx) : mbmi->uv_mode;

    if (xd->lossless[mbmi->segment_id] || tx_size >= TX_32X32)
      return DCT_DCT;

#if CONFIG_EXT_TX
#if ALLOW_INTRA_EXT_TX
    if (mbmi->sb_type >= BLOCK_8X8 && plane_type == PLANE_TYPE_Y)
      return mbmi->tx_type;
#endif  // ALLOW_INTRA_EXT_TX
#endif  // CONFIG_EXT_TX

    if (use_ext_intra_mode_info)
      return filter_intra_mode_to_tx_type_lookup[ext_intra_mode];

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
  }
#endif  // CONFIG_EXT_INTRA

#if CONFIG_EXT_TX
#if EXT_TX_SIZES == 4
  if (xd->lossless[mbmi->segment_id] || tx_size > TX_32X32 ||
      (tx_size >= TX_32X32 && !is_inter_block(mbmi)))
#else
  if (xd->lossless[mbmi->segment_id] || tx_size >= TX_32X32)
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
      return (mbmi->tx_type == IDTX && tx_size == TX_32X32 ?
              DCT_DCT : mbmi->tx_type);
  }

  // Sub8x8-Inter/Intra OR UV-Intra
  if (is_inter_block(mbmi))  // Sub8x8-Inter
    return DCT_DCT;
  else  // Sub8x8 Intra OR UV-Intra
    return intra_mode_to_tx_type_context[plane_type == PLANE_TYPE_Y ?
        get_y_mode(mi, block_idx) : mbmi->uv_mode];
#else
  (void) block_idx;
  if (plane_type != PLANE_TYPE_Y || xd->lossless[mbmi->segment_id] ||
      tx_size >= TX_32X32)
    return DCT_DCT;
  return mbmi->tx_type;
#endif  // CONFIG_EXT_TX
}

void vp10_setup_block_planes(MACROBLOCKD *xd, int ss_x, int ss_y);

static INLINE TX_SIZE get_uv_tx_size_impl(TX_SIZE y_tx_size, BLOCK_SIZE bsize,
                                          int xss, int yss) {
  if (bsize < BLOCK_8X8) {
    return TX_4X4;
  } else {
    const BLOCK_SIZE plane_bsize = ss_size_lookup[bsize][xss][yss];
    return VPXMIN(y_tx_size, max_txsize_lookup[plane_bsize]);
  }
}

static INLINE TX_SIZE get_uv_tx_size(const MB_MODE_INFO *mbmi,
                                     const struct macroblockd_plane *pd) {
#if CONFIG_SUPERTX
  if (supertx_enabled(mbmi))
    return uvsupertx_size_lookup[mbmi->tx_size][pd->subsampling_x]
                                               [pd->subsampling_y];
#endif  // CONFIG_SUPERTX
  return get_uv_tx_size_impl(mbmi->tx_size, mbmi->sb_type, pd->subsampling_x,
                             pd->subsampling_y);
}

static INLINE BLOCK_SIZE get_plane_block_size(BLOCK_SIZE bsize,
    const struct macroblockd_plane *pd) {
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
                                                  TX_SIZE tx_size,
                                                  void *arg);

void vp10_foreach_transformed_block_in_plane(
    const MACROBLOCKD *const xd, BLOCK_SIZE bsize, int plane,
    foreach_transformed_block_visitor visit, void *arg);

void vp10_foreach_transformed_block(
    const MACROBLOCKD* const xd, BLOCK_SIZE bsize,
    foreach_transformed_block_visitor visit, void *arg);

void vp10_set_contexts(const MACROBLOCKD *xd, struct macroblockd_plane *pd,
                      BLOCK_SIZE plane_bsize, TX_SIZE tx_size, int has_eob,
                      int aoff, int loff);

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
  return is_interintra_allowed_bsize(mbmi->sb_type)
          && is_interintra_allowed_mode(mbmi->mode)
          && is_interintra_allowed_ref(mbmi->ref_frame);
}

static INLINE int is_interintra_allowed_bsize_group(const int group) {
  int i;
  for (i = 0; i < BLOCK_SIZES; i++) {
    if (size_group_lookup[i] == group &&
        is_interintra_allowed_bsize(i))
      return 1;
  }
  return 0;
}

static INLINE int is_interintra_pred(const MB_MODE_INFO *mbmi) {
  return (mbmi->ref_frame[1] == INTRA_FRAME) && is_interintra_allowed(mbmi);
}
#endif  // CONFIG_EXT_INTER

#if CONFIG_OBMC || CONFIG_WARPED_MOTION
static INLINE int is_motvar_allowed(const MB_MODE_INFO *mbmi) {
#if CONFIG_EXT_INTER
  return (mbmi->sb_type >= BLOCK_8X8 && mbmi->ref_frame[1] != INTRA_FRAME);
#else
  return (mbmi->sb_type >= BLOCK_8X8);
#endif  // CONFIG_EXT_INTER
}

#if CONFIG_OBMC
static INLINE int is_neighbor_overlappable(const MB_MODE_INFO *mbmi) {
  return (is_inter_block(mbmi));
}
#endif  // CONFIG_OBMC
#endif  // CONFIG_OBMC || CONFIG_WARPED_MOTION

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_COMMON_BLOCKD_H_
