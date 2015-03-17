/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP9_COMMON_VP9_BLOCKD_H_
#define VP9_COMMON_VP9_BLOCKD_H_

#include "./vpx_config.h"

#include "vpx_ports/mem.h"
#include "vpx_scale/yv12config.h"

#include "vp9/common/vp9_common_data.h"
#include "vp9/common/vp9_filter.h"
#include "vp9/common/vp9_mv.h"
#include "vp9/common/vp9_quant_common.h"
#include "vp9/common/vp9_scale.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BLOCK_SIZE_GROUPS 4
#define SKIP_CONTEXTS 3
#define INTER_MODE_CONTEXTS 7

#if CONFIG_COPY_MODE
#define COPY_MODE_CONTEXTS 5
#endif  // CONFIG_COPY_MODE

#if CONFIG_PALETTE
#define PALETTE_BUF_SIZE 16
#define PALETTE_MAX_SIZE 8
#define PALETTE_DELTA_BIT 0
#define PALETTE_COLOR_CONTEXTS 16
#endif  // CONFIG_PALETTE

/* Segment Feature Masks */
#define MAX_MV_REF_CANDIDATES 2

#define INTRA_INTER_CONTEXTS 4
#define COMP_INTER_CONTEXTS 5
#define REF_CONTEXTS 5

typedef enum {
  PLANE_TYPE_Y  = 0,
  PLANE_TYPE_UV = 1,
  PLANE_TYPES
} PLANE_TYPE;

#define MAX_MB_PLANE 3

typedef char ENTROPY_CONTEXT;

static INLINE int combine_entropy_contexts(ENTROPY_CONTEXT a,
                                           ENTROPY_CONTEXT b) {
  return (a != 0) + (b != 0);
}

typedef enum {
  KEY_FRAME = 0,
  INTER_FRAME = 1,
  FRAME_TYPES,
} FRAME_TYPE;

typedef enum {
  DC_PRED,         // Average of above and left pixels
  V_PRED,          // Vertical
  H_PRED,          // Horizontal
  D45_PRED,        // Directional 45  deg = round(arctan(1/1) * 180/pi)
  D135_PRED,       // Directional 135 deg = 180 - 45
  D117_PRED,       // Directional 117 deg = 180 - 63
  D153_PRED,       // Directional 153 deg = 180 - 27
  D207_PRED,       // Directional 207 deg = 180 + 27
  D63_PRED,        // Directional 63  deg = round(arctan(2/1) * 180/pi)
  TM_PRED,         // True-motion
  NEARESTMV,
  NEARMV,
  ZEROMV,
  NEWMV,
#if CONFIG_NEWMVREF
  NEAR_FORNEWMV,
#endif  // CONFIG_NEWMVREF
#if CONFIG_COMPOUND_MODES
  NEAREST_NEARESTMV,
  NEAREST_NEARMV,
  NEAR_NEARESTMV,
  NEAREST_NEWMV,
  NEW_NEARESTMV,
  NEAR_NEWMV,
  NEW_NEARMV,
  ZERO_ZEROMV,
  NEW_NEWMV,
#endif
  MB_MODE_COUNT
} PREDICTION_MODE;

#if CONFIG_COPY_MODE
typedef enum {
  NOREF,
  REF0,
  REF1,
  REF2,
  COPY_MODE_COUNT
} COPY_MODE;
#endif  // CONFIG_COPY_MODE

static INLINE int is_inter_mode(PREDICTION_MODE mode) {
#if CONFIG_NEWMVREF
  return mode >= NEARESTMV && mode <= NEAR_FORNEWMV;
#else
  return mode >= NEARESTMV && mode <= NEWMV;
#endif  // CONFIG_NEWMVREF
}

#if CONFIG_COMPOUND_MODES
static INLINE int is_inter_compound_mode(PREDICTION_MODE mode) {
  return mode >= NEAREST_NEARESTMV && mode <= NEW_NEWMV;
}
#endif

static INLINE int have_newmv_in_inter_mode(PREDICTION_MODE mode) {
#if CONFIG_COMPOUND_MODES
  return (mode == NEWMV ||
#if CONFIG_NEWMVREF
          mode == NEAR_FORNEWMV ||
#endif  // CONFIG_NEWMVREF
          mode == NEW_NEWMV ||
          mode == NEAREST_NEWMV ||
          mode == NEW_NEARESTMV ||
          mode == NEAR_NEWMV ||
          mode == NEW_NEARMV);
#else
#if CONFIG_NEWMVREF
  return (mode == NEWMV ||
          mode == NEAR_FORNEWMV);
#else
  return (mode == NEWMV);
#endif  // CONFIG_NEWMVREF
#endif  // CONFIG_COMPOUND_MODES
}

#define INTRA_MODES (TM_PRED + 1)

#if CONFIG_NEWMVREF
#define INTER_MODES (1 + NEAR_FORNEWMV - NEARESTMV)
#else
#define INTER_MODES (1 + NEWMV - NEARESTMV)
#endif  // CONFIG_NEWMVREF

#define INTER_OFFSET(mode) ((mode) - NEARESTMV)

#if CONFIG_COMPOUND_MODES

#define INTER_COMPOUND_MODES (1 + NEW_NEWMV - NEAREST_NEARESTMV)

#define INTER_COMPOUND_OFFSET(mode) ((mode) - NEAREST_NEARESTMV)

#endif

#if CONFIG_TX64X64
#define MAXTXLEN 64
#else
#define MAXTXLEN 32
#endif

/* For keyframes, intra block modes are predicted by the (already decoded)
   modes for the Y blocks to the left and above us; for interframes, there
   is a single probability table. */

typedef struct {
  PREDICTION_MODE as_mode;
  int_mv as_mv[2];  // first, second inter predictor motion vectors
#if CONFIG_NEWMVREF
  int_mv ref_mv[2];
#endif  // CONFIG_NEWMVREF
} b_mode_info;

// Note that the rate-distortion optimization loop, bit-stream writer, and
// decoder implementation modules critically rely on the enum entry values
// specified herein. They should be refactored concurrently.
typedef enum {
  NONE = -1,
  INTRA_FRAME = 0,
  LAST_FRAME = 1,
  GOLDEN_FRAME = 2,
  ALTREF_FRAME = 3,
  MAX_REF_FRAMES = 4
} MV_REFERENCE_FRAME;

// This structure now relates to 8x8 block regions.
typedef struct {
  // Common for both INTER and INTRA blocks
  BLOCK_SIZE sb_type;
  PREDICTION_MODE mode;
#if CONFIG_FILTERINTRA
  int filterbit, uv_filterbit;
#endif
  TX_SIZE tx_size;
  int8_t skip;
  int8_t segment_id;
  int8_t seg_id_predicted;  // valid only when temporal_update is enabled

  // Only for INTRA blocks
  PREDICTION_MODE uv_mode;

  // Only for INTER blocks
  MV_REFERENCE_FRAME ref_frame[2];
  int_mv mv[2];
  int_mv ref_mvs[MAX_REF_FRAMES][MAX_MV_REF_CANDIDATES];
  uint8_t mode_context[MAX_REF_FRAMES];
  INTERP_FILTER interp_filter;

#if CONFIG_EXT_TX
  EXT_TX_TYPE ext_txfrm;
#endif
#if CONFIG_TX_SKIP
  int tx_skip[PLANE_TYPES];
  int tx_skip_shift;
#endif  // CONFIG_TX_SKIP
#if CONFIG_COPY_MODE
  COPY_MODE copy_mode;
  int inter_ref_count;
#endif  // CONFIG_COPY_MODE
#if CONFIG_INTERINTRA
  PREDICTION_MODE interintra_mode;
  PREDICTION_MODE interintra_uv_mode;
#if CONFIG_WEDGE_PARTITION
  int use_wedge_interintra;
  int interintra_wedge_index;
  int interintra_uv_wedge_index;
#endif  // CONFIG_WEDGE_PARTITION
#endif  // CONFIG_INTERINTRA
#if CONFIG_WEDGE_PARTITION
  int use_wedge_interinter;
  int interinter_wedge_index;
#endif  // CONFIG_WEDGE_PARTITION
#if CONFIG_PALETTE
  int palette_enabled[2];
  int palette_size[2];
  int palette_indexed_size;
  int palette_literal_size;
  int current_palette_size;
  int palette_delta_bitdepth;
  uint8_t palette_colors[3 * PALETTE_MAX_SIZE];
  uint8_t palette_indexed_colors[PALETTE_MAX_SIZE];
  int8_t palette_color_delta[PALETTE_MAX_SIZE];
  uint8_t palette_literal_colors[PALETTE_MAX_SIZE];
  uint8_t *palette_color_map;
  uint8_t *palette_uv_color_map;
#endif  // CONFIG_PALETTE
} MB_MODE_INFO;

typedef struct MODE_INFO {
  struct MODE_INFO *src_mi;
  MB_MODE_INFO mbmi;
#if CONFIG_FILTERINTRA
  int b_filter_info[4];
#endif
  b_mode_info bmi[4];
} MODE_INFO;

static INLINE PREDICTION_MODE get_y_mode(const MODE_INFO *mi, int block) {
  return mi->mbmi.sb_type < BLOCK_8X8 ? mi->bmi[block].as_mode
                                      : mi->mbmi.mode;
}

#if CONFIG_FILTERINTRA
static INLINE int is_filter_allowed(PREDICTION_MODE mode) {
  (void)mode;
  return 1;
}

static INLINE int is_filter_enabled(TX_SIZE txsize) {
  return (txsize < TX_SIZES);
}
#endif

static INLINE int is_inter_block(const MB_MODE_INFO *mbmi) {
  return mbmi->ref_frame[0] > INTRA_FRAME;
}

static INLINE int has_second_ref(const MB_MODE_INFO *mbmi) {
  return mbmi->ref_frame[1] > INTRA_FRAME;
}

PREDICTION_MODE vp9_left_block_mode(const MODE_INFO *cur_mi,
                                    const MODE_INFO *left_mi, int b);

PREDICTION_MODE vp9_above_block_mode(const MODE_INFO *cur_mi,
                                     const MODE_INFO *above_mi, int b);

enum mv_precision {
  MV_PRECISION_Q3,
  MV_PRECISION_Q4
};

struct buf_2d {
  uint8_t *buf;
  int stride;
};

struct macroblockd_plane {
  tran_low_t *dqcoeff;
  PLANE_TYPE plane_type;
  int subsampling_x;
  int subsampling_y;
  struct buf_2d dst;
  struct buf_2d pre[2];
  const int16_t *dequant;
#if CONFIG_NEW_QUANT
  const dequant_val_type_nuq *dequant_val_nuq;
#endif
  ENTROPY_CONTEXT *above_context;
  ENTROPY_CONTEXT *left_context;
#if CONFIG_PALETTE
  uint8_t *color_index_map;
#endif
};

#define BLOCK_OFFSET(x, i) ((x) + (i) * 16)

typedef struct RefBuffer {
  // TODO(dkovalev): idx is not really required and should be removed, now it
  // is used in vp9_onyxd_if.c
  int idx;
  YV12_BUFFER_CONFIG *buf;
  struct scale_factors sf;
} RefBuffer;

typedef struct macroblockd {
  struct macroblockd_plane plane[MAX_MB_PLANE];

  int mi_stride;

  MODE_INFO *mi;

  int up_available;
  int left_available;

  /* Distance of MB away from frame edges */
  int mb_to_left_edge;
  int mb_to_right_edge;
  int mb_to_top_edge;
  int mb_to_bottom_edge;

  /* pointers to reference frames */
  RefBuffer *block_refs[2];

  /* pointer to current frame */
  const YV12_BUFFER_CONFIG *cur_buf;

  /* mc buffer */
  DECLARE_ALIGNED(16, uint8_t, mc_buf[80 * 2 * 80 * 2]);

#if CONFIG_VP9_HIGHBITDEPTH
  /* Bit depth: 8, 10, 12 */
  int bd;
  DECLARE_ALIGNED(16, uint16_t, mc_buf_high[80 * 2 * 80 * 2]);
#endif

  int lossless;

  int corrupted;

  DECLARE_ALIGNED(16, tran_low_t, dqcoeff[MAX_MB_PLANE][64 * 64]);
#if CONFIG_PALETTE
  DECLARE_ALIGNED(16, uint8_t, color_index_map[2][64 * 64]);
  DECLARE_ALIGNED(16, uint8_t, palette_map_buffer[64 * 64]);
#endif  // CONFIG_PALETTE

  ENTROPY_CONTEXT *above_context[MAX_MB_PLANE];
  ENTROPY_CONTEXT left_context[MAX_MB_PLANE][16];

  PARTITION_CONTEXT *above_seg_context;
  PARTITION_CONTEXT left_seg_context[8];
} MACROBLOCKD;

static INLINE BLOCK_SIZE get_subsize(BLOCK_SIZE bsize,
                                     PARTITION_TYPE partition) {
  return subsize_lookup[partition][bsize];
}

extern const TX_TYPE intra_mode_to_tx_type_lookup[INTRA_MODES];

#if CONFIG_SUPERTX

#define PARTITION_SUPERTX_CONTEXTS 2

#if CONFIG_TX64X64
#define MAX_SUPERTX_BLOCK_SIZE BLOCK_64X64
#else
#define MAX_SUPERTX_BLOCK_SIZE BLOCK_32X32
#endif  // CONFIG_TX64X64

static INLINE TX_SIZE bsize_to_tx_size(BLOCK_SIZE bsize) {
  const TX_SIZE bsize_to_tx_size_lookup[BLOCK_SIZES] = {
    TX_4X4, TX_4X4, TX_4X4,
    TX_8X8, TX_8X8, TX_8X8,
    TX_16X16, TX_16X16, TX_16X16,
    TX_32X32, TX_32X32, TX_32X32,
#if CONFIG_TX64X64
    TX_64X64
#else
    TX_32X32
#endif
  };
  return bsize_to_tx_size_lookup[bsize];
}

static INLINE int supertx_enabled(const MB_MODE_INFO *mbmi) {
  return (int)mbmi->tx_size >
         MIN(b_width_log2_lookup[mbmi->sb_type],
             b_height_log2_lookup[mbmi->sb_type]);
}
#endif  // CONFIG_SUPERTX

#if CONFIG_EXT_TX
static TX_TYPE ext_tx_to_txtype[EXT_TX_TYPES] = {
  DCT_DCT,
  ADST_ADST,
  FLIPADST_FLIPADST,
  ADST_FLIPADST,
  FLIPADST_ADST,
  ADST_DCT,
  DCT_ADST,
  FLIPADST_DCT,
  DCT_FLIPADST,
};
#endif  // CONFIG_EXT_TX

static INLINE TX_TYPE get_tx_type(PLANE_TYPE plane_type,
                                  const MACROBLOCKD *xd) {
  const MB_MODE_INFO *const mbmi = &xd->mi[0].src_mi->mbmi;

#if CONFIG_EXT_TX
  if (plane_type != PLANE_TYPE_Y || xd->lossless)
      return DCT_DCT;

  if (is_inter_block(mbmi)) {
    return ext_tx_to_txtype[mbmi->ext_txfrm];
  }
#else
  if (plane_type != PLANE_TYPE_Y || xd->lossless || is_inter_block(mbmi))
    return DCT_DCT;
#endif
  return intra_mode_to_tx_type_lookup[mbmi->mode];
}

static INLINE TX_TYPE get_tx_type_4x4(PLANE_TYPE plane_type,
                                      const MACROBLOCKD *xd, int ib) {
  const MODE_INFO *const mi = xd->mi[0].src_mi;

#if CONFIG_EXT_TX
  if (plane_type != PLANE_TYPE_Y || xd->lossless)
      return DCT_DCT;

  if (is_inter_block(&mi->mbmi)) {
    return ext_tx_to_txtype[mi->mbmi.ext_txfrm];
  }
#else
  if (plane_type != PLANE_TYPE_Y || xd->lossless || is_inter_block(&mi->mbmi))
    return DCT_DCT;
#endif

  return intra_mode_to_tx_type_lookup[get_y_mode(mi, ib)];
}

void vp9_setup_block_planes(MACROBLOCKD *xd, int ss_x, int ss_y);

static INLINE TX_SIZE get_uv_tx_size_impl(TX_SIZE y_tx_size, BLOCK_SIZE bsize,
                                          int xss, int yss) {
  if (bsize < BLOCK_8X8) {
    return TX_4X4;
  } else {
    const BLOCK_SIZE plane_bsize = ss_size_lookup[bsize][xss][yss];
    return MIN(y_tx_size, max_txsize_lookup[plane_bsize]);
  }
}

static INLINE TX_SIZE get_uv_tx_size(const MB_MODE_INFO *mbmi,
                                     const struct macroblockd_plane *pd) {
#if CONFIG_SUPERTX
  if (!supertx_enabled(mbmi)) {
    return get_uv_tx_size_impl(mbmi->tx_size, mbmi->sb_type, pd->subsampling_x,
                               pd->subsampling_y);
  } else {
    return uvsupertx_size_lookup[mbmi->tx_size][pd->subsampling_x]
                                               [pd->subsampling_y];
  }
#else
  return get_uv_tx_size_impl(mbmi->tx_size, mbmi->sb_type, pd->subsampling_x,
                             pd->subsampling_y);
#endif  // CONFIG_SUPERTX
}

static INLINE BLOCK_SIZE get_plane_block_size(BLOCK_SIZE bsize,
    const struct macroblockd_plane *pd) {
  return ss_size_lookup[bsize][pd->subsampling_x][pd->subsampling_y];
}

typedef void (*foreach_transformed_block_visitor)(int plane, int block,
                                                  BLOCK_SIZE plane_bsize,
                                                  TX_SIZE tx_size,
                                                  void *arg);

void vp9_foreach_transformed_block_in_plane(
    const MACROBLOCKD *const xd, BLOCK_SIZE bsize, int plane,
    foreach_transformed_block_visitor visit, void *arg);


void vp9_foreach_transformed_block(
    const MACROBLOCKD* const xd, BLOCK_SIZE bsize,
    foreach_transformed_block_visitor visit, void *arg);

static INLINE void txfrm_block_to_raster_xy(BLOCK_SIZE plane_bsize,
                                            TX_SIZE tx_size, int block,
                                            int *x, int *y) {
  const int bwl = b_width_log2_lookup[plane_bsize];
  const int tx_cols_log2 = bwl - tx_size;
  const int tx_cols = 1 << tx_cols_log2;
  const int raster_mb = block >> (tx_size << 1);
  *x = (raster_mb & (tx_cols - 1)) << tx_size;
  *y = (raster_mb >> tx_cols_log2) << tx_size;
}

void vp9_set_contexts(const MACROBLOCKD *xd, struct macroblockd_plane *pd,
                      BLOCK_SIZE plane_bsize, TX_SIZE tx_size, int has_eob,
                      int aoff, int loff);

#if CONFIG_INTERINTRA
static INLINE int is_interintra_allowed(BLOCK_SIZE sb_type) {
  return ((sb_type >= BLOCK_8X8) && (sb_type < BLOCK_64X64));
}
#endif  // CONFIG_INTERINTRA

#if CONFIG_WEDGE_PARTITION
#define WEDGE_BITS_SML   3
#define WEDGE_BITS_MED   4
#define WEDGE_BITS_BIG   5
#define WEDGE_NONE      -1

#define WEDGE_WEIGHT_BITS 6

static inline int get_wedge_bits(BLOCK_SIZE sb_type) {
  if (sb_type < BLOCK_8X8)
    return 0;
  if (sb_type <= BLOCK_8X8)
    return WEDGE_BITS_SML;
  else if (sb_type <= BLOCK_32X32)
    return WEDGE_BITS_MED;
  else
    return WEDGE_BITS_BIG;
}
#endif  // CONFIG_WEDGE_PARTITION

#if CONFIG_NEW_QUANT && CONFIG_TX_SKIP
static inline int is_rect_quant_used(const MB_MODE_INFO *mbmi,
                                     int plane) {
  return
      mbmi->tx_skip[plane != 0] &&
      ((plane == 0 && (mbmi->mode == V_PRED ||
                       mbmi->mode == H_PRED ||
                       mbmi->mode == TM_PRED)) ||
       (plane != 0 && (mbmi->uv_mode == V_PRED ||
                       mbmi->uv_mode == H_PRED ||
                       mbmi->uv_mode == TM_PRED)));
}
#endif  // CONFIG_NEW_QUANT && CONFIG_TX_SKIP

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_BLOCKD_H_
