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
#include "vpx_scale/yv12config.h"
#include "vp9/common/vp9_convolve.h"
#include "vp9/common/vp9_mv.h"
#include "vp9/common/vp9_treecoder.h"
#include "vpx_ports/mem.h"
#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_enums.h"

// #define MODE_STATS

#define MAX_MB_SEGMENTS     8
#define MB_SEG_TREE_PROBS   (MAX_MB_SEGMENTS-1)
#define PREDICTION_PROBS 3

#define DEFAULT_PRED_PROB_0 120
#define DEFAULT_PRED_PROB_1 80
#define DEFAULT_PRED_PROB_2 40

#define MBSKIP_CONTEXTS 3

#define MAX_REF_LF_DELTAS       4
#define MAX_MODE_LF_DELTAS      4

/* Segment Feature Masks */
#define SEGMENT_DELTADATA   0
#define SEGMENT_ABSDATA     1
#define MAX_MV_REFS 9
#define MAX_MV_REF_CANDIDATES 2

typedef enum {
  PLANE_TYPE_Y_WITH_DC,
  PLANE_TYPE_UV,
} PLANE_TYPE;

typedef char ENTROPY_CONTEXT;

typedef char PARTITION_CONTEXT;

static INLINE int combine_entropy_contexts(ENTROPY_CONTEXT a,
                                           ENTROPY_CONTEXT b) {
  return (a != 0) + (b != 0);
}

typedef enum {
  KEY_FRAME = 0,
  INTER_FRAME = 1
} FRAME_TYPE;

typedef enum {
  EIGHTTAP_SMOOTH,
  EIGHTTAP,
  EIGHTTAP_SHARP,
  BILINEAR,
  SWITCHABLE  /* should be the last one */
} INTERPOLATIONFILTERTYPE;

typedef enum {
  DC_PRED,         // Average of above and left pixels
  V_PRED,          // Vertical
  H_PRED,          // Horizontal
  D45_PRED,        // Directional 45  deg = round(arctan(1/1) * 180/pi)
  D135_PRED,       // Directional 135 deg = 180 - 45
  D117_PRED,       // Directional 117 deg = 180 - 63
  D153_PRED,       // Directional 153 deg = 180 - 27
  D27_PRED,        // Directional 27  deg = round(arctan(1/2) * 180/pi)
  D63_PRED,        // Directional 63  deg = round(arctan(2/1) * 180/pi)
  TM_PRED,         // True-motion
  I4X4_PRED,       // Each 4x4 subblock has its own mode
  NEARESTMV,
  NEARMV,
  ZEROMV,
  NEWMV,
  SPLITMV,
  MB_MODE_COUNT
} MB_PREDICTION_MODE;

static INLINE int is_inter_mode(MB_PREDICTION_MODE mode) {
  return mode >= NEARESTMV && mode <= SPLITMV;
}


// Segment level features.
typedef enum {
  SEG_LVL_ALT_Q = 0,               // Use alternate Quantizer ....
  SEG_LVL_ALT_LF = 1,              // Use alternate loop filter value...
  SEG_LVL_REF_FRAME = 2,           // Optional Segment reference frame
  SEG_LVL_SKIP = 3,                // Optional Segment (0,0) + skip mode
  SEG_LVL_MAX = 4                  // Number of MB level features supported
} SEG_LVL_FEATURES;

// Segment level features.
typedef enum {
  TX_4X4 = 0,                      // 4x4 dct transform
  TX_8X8 = 1,                      // 8x8 dct transform
  TX_16X16 = 2,                    // 16x16 dct transform
  TX_SIZE_MAX_MB = 3,              // Number of different transforms available
  TX_32X32 = TX_SIZE_MAX_MB,       // 32x32 dct transform
  TX_SIZE_MAX_SB,                  // Number of transforms available to SBs
} TX_SIZE;

typedef enum {
  DCT_DCT   = 0,                      // DCT  in both horizontal and vertical
  ADST_DCT  = 1,                      // ADST in vertical, DCT in horizontal
  DCT_ADST  = 2,                      // DCT  in vertical, ADST in horizontal
  ADST_ADST = 3                       // ADST in both directions
} TX_TYPE;

#define VP9_YMODES  (I4X4_PRED + 1)
#define VP9_UV_MODES (TM_PRED + 1)
#define VP9_I32X32_MODES (TM_PRED + 1)

#define VP9_MVREFS (1 + SPLITMV - NEARESTMV)

#define WHT_UPSCALE_FACTOR 2

typedef enum {
  B_DC_PRED,          /* average of above and left pixels */
  B_V_PRED,          /* vertical prediction */
  B_H_PRED,          /* horizontal prediction */
  B_D45_PRED,
  B_D135_PRED,
  B_D117_PRED,
  B_D153_PRED,
  B_D27_PRED,
  B_D63_PRED,
  B_TM_PRED,

  LEFT4X4,
  ABOVE4X4,
  ZERO4X4,
  NEW4X4,

  B_MODE_COUNT
} B_PREDICTION_MODE;

#define VP9_BINTRAMODES (LEFT4X4)
#define VP9_SUBMVREFS (1 + NEW4X4 - LEFT4X4)

#define VP9_KF_BINTRAMODES (VP9_BINTRAMODES)   /* 10 */
#define VP9_NKF_BINTRAMODES (VP9_BINTRAMODES)  /* 10 */

/* For keyframes, intra block modes are predicted by the (already decoded)
   modes for the Y blocks to the left and above us; for interframes, there
   is a single probability table. */

union b_mode_info {
  struct {
    B_PREDICTION_MODE first;
  } as_mode;
  int_mv as_mv[2];  // first, second inter predictor motion vectors
};

typedef enum {
  NONE = -1,
  INTRA_FRAME = 0,
  LAST_FRAME = 1,
  GOLDEN_FRAME = 2,
  ALTREF_FRAME = 3,
  MAX_REF_FRAMES = 4
} MV_REFERENCE_FRAME;

static INLINE int b_width_log2(BLOCK_SIZE_TYPE sb_type) {
  switch (sb_type) {
#if CONFIG_AB4X4
    case BLOCK_SIZE_SB4X8:
#endif
    case BLOCK_SIZE_AB4X4: return 0;
#if CONFIG_AB4X4
    case BLOCK_SIZE_SB8X4:
#endif
    case BLOCK_SIZE_SB8X8:
    case BLOCK_SIZE_SB8X16: return 1;
    case BLOCK_SIZE_SB16X8:
    case BLOCK_SIZE_MB16X16:
    case BLOCK_SIZE_SB16X32: return 2;
    case BLOCK_SIZE_SB32X16:
    case BLOCK_SIZE_SB32X32:
    case BLOCK_SIZE_SB32X64: return 3;
    case BLOCK_SIZE_SB64X32:
    case BLOCK_SIZE_SB64X64: return 4;
    default: assert(0);
  }
}

static INLINE int b_height_log2(BLOCK_SIZE_TYPE sb_type) {
  switch (sb_type) {
#if CONFIG_AB4X4
    case BLOCK_SIZE_SB8X4:
#endif
    case BLOCK_SIZE_AB4X4: return 0;
#if CONFIG_AB4X4
    case BLOCK_SIZE_SB4X8:
#endif
    case BLOCK_SIZE_SB8X8:
    case BLOCK_SIZE_SB16X8: return 1;
    case BLOCK_SIZE_SB8X16:
    case BLOCK_SIZE_MB16X16:
    case BLOCK_SIZE_SB32X16: return 2;
    case BLOCK_SIZE_SB16X32:
    case BLOCK_SIZE_SB32X32:
    case BLOCK_SIZE_SB64X32: return 3;
    case BLOCK_SIZE_SB32X64:
    case BLOCK_SIZE_SB64X64: return 4;
    default: assert(0);
  }
}

static INLINE int mi_width_log2(BLOCK_SIZE_TYPE sb_type) {
  int a = b_width_log2(sb_type) - 1;
  assert(a >= 0);
  return a;
}

static INLINE int mi_height_log2(BLOCK_SIZE_TYPE sb_type) {
  int a = b_height_log2(sb_type) - 1;
  assert(a >= 0);
  return a;
}

typedef struct {
  MB_PREDICTION_MODE mode, uv_mode;
  MV_REFERENCE_FRAME ref_frame, second_ref_frame;
  TX_SIZE txfm_size;
  int_mv mv[2]; // for each reference frame used
  int_mv ref_mvs[MAX_REF_FRAMES][MAX_MV_REF_CANDIDATES];
  int_mv best_mv, best_second_mv;

  int mb_mode_context[MAX_REF_FRAMES];

  unsigned char mb_skip_coeff;                                /* does this mb has coefficients at all, 1=no coefficients, 0=need decode tokens */
  unsigned char need_to_clamp_mvs;
  unsigned char need_to_clamp_secondmv;
  unsigned char segment_id;           // Segment id for current frame

  // Flags used for prediction status of various bistream signals
  unsigned char seg_id_predicted;
  unsigned char ref_predicted;

  // Indicates if the mb is part of the image (1) vs border (0)
  // This can be useful in determining whether the MB provides
  // a valid predictor
  unsigned char mb_in_image;

  INTERPOLATIONFILTERTYPE interp_filter;

  BLOCK_SIZE_TYPE sb_type;
} MB_MODE_INFO;

typedef struct {
  MB_MODE_INFO mbmi;
  union b_mode_info bmi[4];
} MODE_INFO;

struct scale_factors {
  int x_num;
  int x_den;
  int x_offset_q4;
  int x_step_q4;
  int y_num;
  int y_den;
  int y_offset_q4;
  int y_step_q4;

  int (*scale_value_x)(int val, const struct scale_factors *scale);
  int (*scale_value_y)(int val, const struct scale_factors *scale);
  void (*set_scaled_offsets)(struct scale_factors *scale, int row, int col);
  int_mv32 (*scale_motion_vector_q3_to_q4)(const int_mv *src_mv,
                                           const struct scale_factors *scale);
  int32_t (*scale_motion_vector_component_q4)(int mv_q4,
                                              int num,
                                              int den,
                                              int offset_q4);

  convolve_fn_t predict[2][2][2];  // horiz, vert, avg
};

enum { MAX_MB_PLANE = 3 };

struct buf_2d {
  uint8_t *buf;
  int stride;
};

struct macroblockd_plane {
  DECLARE_ALIGNED(16, int16_t,  qcoeff[64 * 64]);
  DECLARE_ALIGNED(16, int16_t,  dqcoeff[64 * 64]);
  DECLARE_ALIGNED(16, uint16_t, eobs[256]);
  DECLARE_ALIGNED(16, int16_t,  diff[64 * 64]);
  PLANE_TYPE plane_type;
  int subsampling_x;
  int subsampling_y;
  struct buf_2d dst;
  struct buf_2d pre[2];
  int16_t *dequant;
  ENTROPY_CONTEXT *above_context;
  ENTROPY_CONTEXT *left_context;
};

#define BLOCK_OFFSET(x, i, n) ((x) + (i) * (n))

typedef struct macroblockd {
  struct macroblockd_plane plane[MAX_MB_PLANE];

  struct scale_factors scale_factor[2];
  struct scale_factors scale_factor_uv[2];

  MODE_INFO *prev_mode_info_context;
  MODE_INFO *mode_info_context;
  int mode_info_stride;

  FRAME_TYPE frame_type;

  int up_available;
  int left_available;
  int right_available;

  // partition contexts
  PARTITION_CONTEXT *above_seg_context;
  PARTITION_CONTEXT *left_seg_context;

  /* 0 (disable) 1 (enable) segmentation */
  unsigned char segmentation_enabled;

  /* 0 (do not update) 1 (update) the macroblock segmentation map. */
  unsigned char update_mb_segmentation_map;

#if CONFIG_IMPLICIT_SEGMENTATION
  unsigned char allow_implicit_segment_update;
#endif

  /* 0 (do not update) 1 (update) the macroblock segmentation feature data. */
  unsigned char update_mb_segmentation_data;

  /* 0 (do not update) 1 (update) the macroblock segmentation feature data. */
  unsigned char mb_segment_abs_delta;

  /* Per frame flags that define which MB level features (such as quantizer or loop filter level) */
  /* are enabled and when enabled the proabilities used to decode the per MB flags in MB_MODE_INFO */

  // Probability Tree used to code Segment number
  vp9_prob mb_segment_tree_probs[MB_SEG_TREE_PROBS];

  // Segment features
  signed char segment_feature_data[MAX_MB_SEGMENTS][SEG_LVL_MAX];
  unsigned int segment_feature_mask[MAX_MB_SEGMENTS];

  /* mode_based Loop filter adjustment */
  unsigned char mode_ref_lf_delta_enabled;
  unsigned char mode_ref_lf_delta_update;

  /* Delta values have the range +/- MAX_LOOP_FILTER */
  /* 0 = Intra, Last, GF, ARF */
  signed char last_ref_lf_deltas[MAX_REF_LF_DELTAS];
  /* 0 = Intra, Last, GF, ARF */
  signed char ref_lf_deltas[MAX_REF_LF_DELTAS];
  /* 0 = I4X4_PRED, ZERO_MV, MV, SPLIT */
  signed char last_mode_lf_deltas[MAX_MODE_LF_DELTAS];
  /* 0 = I4X4_PRED, ZERO_MV, MV, SPLIT */
  signed char mode_lf_deltas[MAX_MODE_LF_DELTAS];

  /* Distance of MB away from frame edges */
  int mb_to_left_edge;
  int mb_to_right_edge;
  int mb_to_top_edge;
  int mb_to_bottom_edge;

  unsigned int frames_since_golden;
  unsigned int frames_till_alt_ref_frame;

  int lossless;
  /* Inverse transform function pointers. */
  void (*inv_txm4x4_1)(int16_t *input, int16_t *output, int pitch);
  void (*inv_txm4x4)(int16_t *input, int16_t *output, int pitch);
  void (*itxm_add)(int16_t *input, uint8_t *dest, int stride, int eob);
  void (*itxm_add_y_block)(int16_t *q, uint8_t *dst, int stride,
    struct macroblockd *xd);
  void (*itxm_add_uv_block)(int16_t *q, uint8_t *dst, int stride,
    uint16_t *eobs);

  struct subpix_fn_table  subpix;

  int allow_high_precision_mv;

  int corrupted;

  int sb_index;   // index of 32x32 block inside the 64x64 block
  int mb_index;   // index of 16x16 block inside the 32x32 block
  int b_index;    // index of 8x8 block inside the 16x16 block
#if CONFIG_AB4X4
  int ab_index;   // index of 4x4 block inside the 8x8 block
#endif
  int q_index;

} MACROBLOCKD;

static int *get_sb_index(MACROBLOCKD *xd, BLOCK_SIZE_TYPE subsize) {
  switch (subsize) {
    case BLOCK_SIZE_SB64X32:
    case BLOCK_SIZE_SB32X64:
    case BLOCK_SIZE_SB32X32:
      return &xd->sb_index;
    case BLOCK_SIZE_SB32X16:
    case BLOCK_SIZE_SB16X32:
    case BLOCK_SIZE_MB16X16:
      return &xd->mb_index;
    case BLOCK_SIZE_SB16X8:
    case BLOCK_SIZE_SB8X16:
    case BLOCK_SIZE_SB8X8:
      return &xd->b_index;
#if CONFIG_AB4X4
    case BLOCK_SIZE_SB8X4:
    case BLOCK_SIZE_SB4X8:
    case BLOCK_SIZE_AB4X4:
      return &xd->ab_index;
#endif
    default:
      assert(0);
      return NULL;
  }
}

static INLINE void update_partition_context(MACROBLOCKD *xd,
                                            BLOCK_SIZE_TYPE sb_type,
                                            BLOCK_SIZE_TYPE sb_size) {
  int bsl = mi_width_log2(sb_size), bs = 1 << bsl;
  int bwl = mi_width_log2(sb_type);
  int bhl = mi_height_log2(sb_type);
  int boffset = mi_width_log2(BLOCK_SIZE_SB64X64) - bsl;
  int i;
  // skip macroblock partition
  if (bsl == 0)
    return;

  // update the partition context at the end notes. set partition bits
  // of block sizes larger than the current one to be one, and partition
  // bits of smaller block sizes to be zero.
  if ((bwl == bsl) && (bhl == bsl)) {
    for (i = 0; i < bs; i++)
      xd->left_seg_context[i] = ~(0xf << boffset);
    for (i = 0; i < bs; i++)
      xd->above_seg_context[i] = ~(0xf << boffset);
  } else if ((bwl == bsl) && (bhl < bsl)) {
    for (i = 0; i < bs; i++)
      xd->left_seg_context[i] = ~(0xe << boffset);
    for (i = 0; i < bs; i++)
      xd->above_seg_context[i] = ~(0xf << boffset);
  }  else if ((bwl < bsl) && (bhl == bsl)) {
    for (i = 0; i < bs; i++)
      xd->left_seg_context[i] = ~(0xf << boffset);
    for (i = 0; i < bs; i++)
      xd->above_seg_context[i] = ~(0xe << boffset);
  } else if ((bwl < bsl) && (bhl < bsl)) {
    for (i = 0; i < bs; i++)
      xd->left_seg_context[i] = ~(0xe << boffset);
    for (i = 0; i < bs; i++)
      xd->above_seg_context[i] = ~(0xe << boffset);
  } else {
    assert(0);
  }
}

static INLINE int partition_plane_context(MACROBLOCKD *xd,
                                          BLOCK_SIZE_TYPE sb_type) {
  int bsl = mi_width_log2(sb_type), bs = 1 << bsl;
  int above = 0, left = 0, i;
  int boffset = mi_width_log2(BLOCK_SIZE_SB64X64) - bsl;

  assert(mi_width_log2(sb_type) == mi_height_log2(sb_type));
  assert(bsl >= 0);
  assert(boffset >= 0);

  for (i = 0; i < bs; i++)
    above |= (xd->above_seg_context[i] & (1 << boffset));
  for (i = 0; i < bs; i++)
    left |= (xd->left_seg_context[i] & (1 << boffset));

  above = (above > 0);
  left  = (left > 0);

  return (left * 2 + above) + (bsl - 1) * PARTITION_PLOFFSET;
}

static BLOCK_SIZE_TYPE get_subsize(BLOCK_SIZE_TYPE bsize,
                                   PARTITION_TYPE partition) {
  BLOCK_SIZE_TYPE subsize;
  switch (partition) {
    case PARTITION_NONE:
      subsize = bsize;
      break;
    case PARTITION_HORZ:
      if (bsize == BLOCK_SIZE_SB64X64)
        subsize = BLOCK_SIZE_SB64X32;
      else if (bsize == BLOCK_SIZE_SB32X32)
        subsize = BLOCK_SIZE_SB32X16;
      else if (bsize == BLOCK_SIZE_MB16X16)
        subsize = BLOCK_SIZE_SB16X8;
      else
        assert(0);
      break;
    case PARTITION_VERT:
      if (bsize == BLOCK_SIZE_SB64X64)
        subsize = BLOCK_SIZE_SB32X64;
      else if (bsize == BLOCK_SIZE_SB32X32)
        subsize = BLOCK_SIZE_SB16X32;
      else if (bsize == BLOCK_SIZE_MB16X16)
        subsize = BLOCK_SIZE_SB8X16;
      else
        assert(0);
      break;
    case PARTITION_SPLIT:
      if (bsize == BLOCK_SIZE_SB64X64)
        subsize = BLOCK_SIZE_SB32X32;
      else if (bsize == BLOCK_SIZE_SB32X32)
        subsize = BLOCK_SIZE_MB16X16;
      else if (bsize == BLOCK_SIZE_MB16X16)
        subsize = BLOCK_SIZE_SB8X8;
      else
        assert(0);
      break;
    default:
      assert(0);
  }
  return subsize;
}

#define ACTIVE_HT   110                // quantization stepsize threshold

#define ACTIVE_HT8  300

#define ACTIVE_HT16 300

// convert MB_PREDICTION_MODE to B_PREDICTION_MODE
static B_PREDICTION_MODE pred_mode_conv(MB_PREDICTION_MODE mode) {
  switch (mode) {
    case DC_PRED: return B_DC_PRED;
    case V_PRED: return B_V_PRED;
    case H_PRED: return B_H_PRED;
    case TM_PRED: return B_TM_PRED;
    case D45_PRED: return B_D45_PRED;
    case D135_PRED: return B_D135_PRED;
    case D117_PRED: return B_D117_PRED;
    case D153_PRED: return B_D153_PRED;
    case D27_PRED: return B_D27_PRED;
    case D63_PRED: return B_D63_PRED;
    default:
       assert(0);
       return B_MODE_COUNT;  // Dummy value
  }
}

// transform mapping
static TX_TYPE txfm_map(B_PREDICTION_MODE bmode) {
  switch (bmode) {
    case B_TM_PRED :
    case B_D135_PRED :
      return ADST_ADST;

    case B_V_PRED :
    case B_D117_PRED :
      return ADST_DCT;

    case B_H_PRED :
    case B_D153_PRED :
    case B_D27_PRED :
      return DCT_ADST;

    default:
      return DCT_DCT;
  }
}

#define USE_ADST_FOR_I16X16_8X8   1
#define USE_ADST_FOR_I16X16_4X4   1
#define USE_ADST_FOR_I8X8_4X4     1
#define USE_ADST_PERIPHERY_ONLY   1
#define USE_ADST_FOR_SB           1
#define USE_ADST_FOR_REMOTE_EDGE  0

static TX_TYPE get_tx_type_4x4(const MACROBLOCKD *xd, int ib) {
  // TODO(debargha): explore different patterns for ADST usage when blocksize
  // is smaller than the prediction size
  TX_TYPE tx_type = DCT_DCT;
  const BLOCK_SIZE_TYPE sb_type = xd->mode_info_context->mbmi.sb_type;
  const int wb = b_width_log2(sb_type), hb = b_height_log2(sb_type);
#if !USE_ADST_FOR_SB
  if (sb_type > BLOCK_SIZE_MB16X16)
    return tx_type;
#endif
  if (ib >= (1 << (wb + hb)))  // no chroma adst
    return tx_type;
  if (xd->lossless)
    return DCT_DCT;
  if (xd->mode_info_context->mbmi.mode == I4X4_PRED &&
      xd->q_index < ACTIVE_HT) {
    tx_type = txfm_map(
        xd->mode_info_context->bmi[ib].as_mode.first);
  } else if (xd->mode_info_context->mbmi.mode <= TM_PRED &&
             xd->q_index < ACTIVE_HT) {
#if USE_ADST_FOR_I16X16_4X4
#if USE_ADST_PERIPHERY_ONLY
    const int hmax = 1 << wb;
    tx_type = txfm_map(pred_mode_conv(xd->mode_info_context->mbmi.mode));
#if USE_ADST_FOR_REMOTE_EDGE
    if ((ib & (hmax - 1)) != 0 && ib >= hmax)
      tx_type = DCT_DCT;
#else
    if (ib >= 1 && ib < hmax) {
      if (tx_type == ADST_ADST) tx_type = ADST_DCT;
      else if (tx_type == DCT_ADST) tx_type = DCT_DCT;
    } else if (ib >= 1 && (ib & (hmax - 1)) == 0) {
      if (tx_type == ADST_ADST) tx_type = DCT_ADST;
      else if (tx_type == ADST_DCT) tx_type = DCT_DCT;
    } else if (ib != 0) {
      tx_type = DCT_DCT;
    }
#endif
#else
    // Use ADST
    tx_type = txfm_map(pred_mode_conv(xd->mode_info_context->mbmi.mode));
#endif
#else
    // Use 2D DCT
    tx_type = DCT_DCT;
#endif
  }
  return tx_type;
}

static TX_TYPE get_tx_type_8x8(const MACROBLOCKD *xd, int ib) {
  // TODO(debargha): explore different patterns for ADST usage when blocksize
  // is smaller than the prediction size
  TX_TYPE tx_type = DCT_DCT;
  const BLOCK_SIZE_TYPE sb_type = xd->mode_info_context->mbmi.sb_type;
  const int wb = b_width_log2(sb_type), hb = b_height_log2(sb_type);
#if !USE_ADST_FOR_SB
  if (sb_type > BLOCK_SIZE_MB16X16)
    return tx_type;
#endif
  if (ib >= (1 << (wb + hb)))  // no chroma adst
    return tx_type;
  if (xd->mode_info_context->mbmi.mode <= TM_PRED &&
      xd->q_index < ACTIVE_HT8) {
#if USE_ADST_FOR_I16X16_8X8
#if USE_ADST_PERIPHERY_ONLY
    const int hmax = 1 << wb;
    tx_type = txfm_map(pred_mode_conv(xd->mode_info_context->mbmi.mode));
#if USE_ADST_FOR_REMOTE_EDGE
    if ((ib & (hmax - 1)) != 0 && ib >= hmax)
      tx_type = DCT_DCT;
#else
    if (ib >= 1 && ib < hmax) {
      if (tx_type == ADST_ADST) tx_type = ADST_DCT;
      else if (tx_type == DCT_ADST) tx_type = DCT_DCT;
    } else if (ib >= 1 && (ib & (hmax - 1)) == 0) {
      if (tx_type == ADST_ADST) tx_type = DCT_ADST;
      else if (tx_type == ADST_DCT) tx_type = DCT_DCT;
    } else if (ib != 0) {
      tx_type = DCT_DCT;
    }
#endif
#else
    // Use ADST
    tx_type = txfm_map(pred_mode_conv(xd->mode_info_context->mbmi.mode));
#endif
#else
    // Use 2D DCT
    tx_type = DCT_DCT;
#endif
  }
  return tx_type;
}

static TX_TYPE get_tx_type_16x16(const MACROBLOCKD *xd, int ib) {
  TX_TYPE tx_type = DCT_DCT;
  const BLOCK_SIZE_TYPE sb_type = xd->mode_info_context->mbmi.sb_type;
  const int wb = b_width_log2(sb_type), hb = b_height_log2(sb_type);
#if !USE_ADST_FOR_SB
  if (sb_type > BLOCK_SIZE_MB16X16)
    return tx_type;
#endif
  if (ib >= (1 << (wb + hb)))
    return tx_type;
  if (xd->mode_info_context->mbmi.mode <= TM_PRED &&
      xd->q_index < ACTIVE_HT16) {
    tx_type = txfm_map(pred_mode_conv(xd->mode_info_context->mbmi.mode));
#if USE_ADST_PERIPHERY_ONLY
    if (sb_type > BLOCK_SIZE_MB16X16) {
      const int hmax = 1 << wb;
#if USE_ADST_FOR_REMOTE_EDGE
      if ((ib & (hmax - 1)) != 0 && ib >= hmax)
        tx_type = DCT_DCT;
#else
      if (ib >= 1 && ib < hmax) {
        if (tx_type == ADST_ADST) tx_type = ADST_DCT;
        else if (tx_type == DCT_ADST) tx_type = DCT_DCT;
      } else if (ib >= 1 && (ib & (hmax - 1)) == 0) {
        if (tx_type == ADST_ADST) tx_type = DCT_ADST;
        else if (tx_type == ADST_DCT) tx_type = DCT_DCT;
      } else if (ib != 0) {
        tx_type = DCT_DCT;
      }
#endif
    }
#endif
  }
  return tx_type;
}

void vp9_setup_block_dptrs(MACROBLOCKD *xd,
                           int subsampling_x, int subsampling_y);

static TX_SIZE get_uv_tx_size(const MACROBLOCKD *xd) {
  MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
  const TX_SIZE size = mbmi->txfm_size;

  switch (mbmi->sb_type) {
    case BLOCK_SIZE_SB64X64:
      return size;
    case BLOCK_SIZE_SB64X32:
    case BLOCK_SIZE_SB32X64:
    case BLOCK_SIZE_SB32X32:
      if (size == TX_32X32)
        return TX_16X16;
      else
        return size;
    case BLOCK_SIZE_SB32X16:
    case BLOCK_SIZE_SB16X32:
    case BLOCK_SIZE_MB16X16:
      if (size == TX_16X16)
        return TX_8X8;
      else
        return size;
    default:
      return TX_4X4;
  }

  return size;
}

struct plane_block_idx {
  int plane;
  int block;
};

// TODO(jkoleszar): returning a struct so it can be used in a const context,
// expect to refactor this further later.
static INLINE struct plane_block_idx plane_block_idx(int y_blocks,
                                                     int b_idx) {
  const int v_offset = y_blocks * 5 / 4;
  struct plane_block_idx res;

  if (b_idx < y_blocks) {
    res.plane = 0;
    res.block = b_idx;
  } else if (b_idx < v_offset) {
    res.plane = 1;
    res.block = b_idx - y_blocks;
  } else {
    assert(b_idx < y_blocks * 3 / 2);
    res.plane = 2;
    res.block = b_idx - v_offset;
  }
  return res;
}

typedef void (*foreach_transformed_block_visitor)(int plane, int block,
                                                  BLOCK_SIZE_TYPE bsize,
                                                  int ss_txfrm_size,
                                                  void *arg);
static INLINE void foreach_transformed_block_in_plane(
    const MACROBLOCKD* const xd, BLOCK_SIZE_TYPE bsize, int plane,
    foreach_transformed_block_visitor visit, void *arg) {
  const int bw = b_width_log2(bsize), bh = b_height_log2(bsize);

  // block and transform sizes, in number of 4x4 blocks log 2 ("*_b")
  // 4x4=0, 8x8=2, 16x16=4, 32x32=6, 64x64=8
  // transform size varies per plane, look it up in a common way.
  const TX_SIZE tx_size = plane ? get_uv_tx_size(xd)
                                : xd->mode_info_context->mbmi.txfm_size;
  const int block_size_b = bw + bh;
  const int txfrm_size_b = tx_size * 2;

  // subsampled size of the block
  const int ss_sum = xd->plane[plane].subsampling_x +
                     xd->plane[plane].subsampling_y;
  const int ss_block_size = block_size_b - ss_sum;

  const int step = 1 << txfrm_size_b;

  int i;

  assert(txfrm_size_b <= block_size_b);
  assert(txfrm_size_b <= ss_block_size);
  for (i = 0; i < (1 << ss_block_size); i += step) {
    visit(plane, i, bsize, txfrm_size_b, arg);
  }
}

static INLINE void foreach_transformed_block(
    const MACROBLOCKD* const xd, BLOCK_SIZE_TYPE bsize,
    foreach_transformed_block_visitor visit, void *arg) {
  int plane;

  for (plane = 0; plane < MAX_MB_PLANE; plane++) {
    foreach_transformed_block_in_plane(xd, bsize, plane,
                                       visit, arg);
  }
}

static INLINE void foreach_transformed_block_uv(
    const MACROBLOCKD* const xd, BLOCK_SIZE_TYPE bsize,
    foreach_transformed_block_visitor visit, void *arg) {
  int plane;

  for (plane = 1; plane < MAX_MB_PLANE; plane++) {
    foreach_transformed_block_in_plane(xd, bsize, plane,
                                       visit, arg);
  }
}

// TODO(jkoleszar): In principle, pred_w, pred_h are unnecessary, as we could
// calculate the subsampled BLOCK_SIZE_TYPE, but that type isn't defined for
// sizes smaller than 16x16 yet.
typedef void (*foreach_predicted_block_visitor)(int plane, int block,
                                                BLOCK_SIZE_TYPE bsize,
                                                int pred_w, int pred_h,
                                                void *arg);
static INLINE void foreach_predicted_block_in_plane(
    const MACROBLOCKD* const xd, BLOCK_SIZE_TYPE bsize, int plane,
    foreach_predicted_block_visitor visit, void *arg) {
  int i, x, y;
  const MB_PREDICTION_MODE mode = xd->mode_info_context->mbmi.mode;

  // block sizes in number of 4x4 blocks log 2 ("*_b")
  // 4x4=0, 8x8=2, 16x16=4, 32x32=6, 64x64=8
  // subsampled size of the block
  const int bw = b_width_log2(bsize) - xd->plane[plane].subsampling_x;
  const int bh = b_height_log2(bsize) - xd->plane[plane].subsampling_y;

  // size of the predictor to use.
  int pred_w, pred_h;

  if (mode == SPLITMV) {
    pred_w = 0;
    pred_h = 0;
  } else {
    pred_w = bw;
    pred_h = bh;
  }
  assert(pred_w <= bw);
  assert(pred_h <= bh);

  // visit each subblock in raster order
  i = 0;
  for (y = 0; y < 1 << bh; y += 1 << pred_h) {
    for (x = 0; x < 1 << bw; x += 1 << pred_w) {
      visit(plane, i, bsize, pred_w, pred_h, arg);
      i += 1 << pred_w;
    }
    i -= 1 << bw;
    i += 1 << (bw + pred_h);
  }
}
static INLINE void foreach_predicted_block(
    const MACROBLOCKD* const xd, BLOCK_SIZE_TYPE bsize,
    foreach_predicted_block_visitor visit, void *arg) {
  int plane;

  for (plane = 0; plane < MAX_MB_PLANE; plane++) {
    foreach_predicted_block_in_plane(xd, bsize, plane, visit, arg);
  }
}
static INLINE void foreach_predicted_block_uv(
    const MACROBLOCKD* const xd, BLOCK_SIZE_TYPE bsize,
    foreach_predicted_block_visitor visit, void *arg) {
  int plane;

  for (plane = 1; plane < MAX_MB_PLANE; plane++) {
    foreach_predicted_block_in_plane(xd, bsize, plane, visit, arg);
  }
}
static int raster_block_offset(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize,
                               int plane, int block, int stride) {
  const int bw = b_width_log2(bsize) - xd->plane[plane].subsampling_x;
  const int y = 4 * (block >> bw), x = 4 * (block & ((1 << bw) - 1));
  return y * stride + x;
}
static int16_t* raster_block_offset_int16(MACROBLOCKD *xd,
                                         BLOCK_SIZE_TYPE bsize,
                                         int plane, int block, int16_t *base) {
  const int bw = b_width_log2(bsize) - xd->plane[plane].subsampling_x;
  const int stride = 4 << bw;
  return base + raster_block_offset(xd, bsize, plane, block, stride);
}
static uint8_t* raster_block_offset_uint8(MACROBLOCKD *xd,
                                         BLOCK_SIZE_TYPE bsize,
                                         int plane, int block,
                                         uint8_t *base, int stride) {
  return base + raster_block_offset(xd, bsize, plane, block, stride);
}

static int txfrm_block_to_raster_block(MACROBLOCKD *xd,
                                       BLOCK_SIZE_TYPE bsize,
                                       int plane, int block,
                                       int ss_txfrm_size) {
  const int bwl = b_width_log2(bsize) - xd->plane[plane].subsampling_x;
  const int txwl = ss_txfrm_size / 2;
  const int tx_cols_lg2 = bwl - txwl;
  const int tx_cols = 1 << tx_cols_lg2;
  const int raster_mb = block >> ss_txfrm_size;
  const int x = (raster_mb & (tx_cols - 1)) << (txwl);
  const int y = raster_mb >> tx_cols_lg2 << (txwl);
  return x + (y << bwl);
}

static void txfrm_block_to_raster_xy(MACROBLOCKD *xd,
                                     BLOCK_SIZE_TYPE bsize,
                                     int plane, int block,
                                     int ss_txfrm_size,
                                     int *x, int *y) {
  const int bwl = b_width_log2(bsize) - xd->plane[plane].subsampling_x;
  const int txwl = ss_txfrm_size / 2;
  const int tx_cols_lg2 = bwl - txwl;
  const int tx_cols = 1 << tx_cols_lg2;
  const int raster_mb = block >> ss_txfrm_size;
  *x = (raster_mb & (tx_cols - 1)) << (txwl);
  *y = raster_mb >> tx_cols_lg2 << (txwl);
}

#endif  // VP9_COMMON_VP9_BLOCKD_H_
