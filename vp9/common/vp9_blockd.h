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

#define BLOCK_SIZE_GROUPS   4
#define MAX_MB_SEGMENTS     8
#define MB_SEG_TREE_PROBS   (MAX_MB_SEGMENTS-1)

#define PREDICTION_PROBS 3

#define MBSKIP_CONTEXTS 3

#define MAX_REF_LF_DELTAS       4
#define MAX_MODE_LF_DELTAS      2

/* Segment Feature Masks */
#define SEGMENT_DELTADATA   0
#define SEGMENT_ABSDATA     1
#define MAX_MV_REF_CANDIDATES 2

#define INTRA_INTER_CONTEXTS 4
#define COMP_INTER_CONTEXTS 5
#define REF_CONTEXTS 5

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
  INTER_FRAME = 1,
  NUM_FRAME_TYPES,
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
  NEARESTMV,
  NEARMV,
  ZEROMV,
  NEWMV,
  MB_MODE_COUNT
} MB_PREDICTION_MODE;

static INLINE int is_inter_mode(MB_PREDICTION_MODE mode) {
  return mode >= NEARESTMV && mode <= NEWMV;
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
  TX_32X32 = 3,                    // 32x32 dct transform
  TX_SIZE_MAX_SB,                  // Number of transforms available to SBs
} TX_SIZE;

typedef enum {
  DCT_DCT   = 0,                      // DCT  in both horizontal and vertical
  ADST_DCT  = 1,                      // ADST in vertical, DCT in horizontal
  DCT_ADST  = 2,                      // DCT  in vertical, ADST in horizontal
  ADST_ADST = 3                       // ADST in both directions
} TX_TYPE;

#define VP9_INTRA_MODES (TM_PRED + 1)

#define VP9_INTER_MODES (1 + NEWMV - NEARESTMV)

#define WHT_UPSCALE_FACTOR 2

#define TX_SIZE_PROBS  6  // (TX_SIZE_MAX_SB * (TX_SIZE_MAX_SB - 1) / 2)

#define get_tx_probs(c, b) ((b) < BLOCK_SIZE_MB16X16 ? \
                            (c)->fc.tx_probs_8x8p :    \
                            (b) < BLOCK_SIZE_SB32X32 ? \
                            (c)->fc.tx_probs_16x16p : (c)->fc.tx_probs_32x32p)

/* For keyframes, intra block modes are predicted by the (already decoded)
   modes for the Y blocks to the left and above us; for interframes, there
   is a single probability table. */

union b_mode_info {
  struct {
    MB_PREDICTION_MODE first;
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
    case BLOCK_SIZE_SB4X8:
    case BLOCK_SIZE_AB4X4: return 0;
    case BLOCK_SIZE_SB8X4:
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
      return -1;
  }
}

static INLINE int b_height_log2(BLOCK_SIZE_TYPE sb_type) {
  switch (sb_type) {
    case BLOCK_SIZE_SB8X4:
    case BLOCK_SIZE_AB4X4: return 0;
    case BLOCK_SIZE_SB4X8:
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
      return -1;
  }
}

static INLINE int mi_width_log2(BLOCK_SIZE_TYPE sb_type) {
  int a = b_width_log2(sb_type) - 1;
  // align 4x4 block to mode_info
  if (a < 0)
    a = 0;
  assert(a >= 0);
  return a;
}

static INLINE int mi_height_log2(BLOCK_SIZE_TYPE sb_type) {
  int a = b_height_log2(sb_type) - 1;
  if (a < 0)
    a = 0;
  assert(a >= 0);
  return a;
}

typedef struct {
  MB_PREDICTION_MODE mode, uv_mode;
  MV_REFERENCE_FRAME ref_frame[2];
  TX_SIZE txfm_size;
  int_mv mv[2]; // for each reference frame used
  int_mv ref_mvs[MAX_REF_FRAMES][MAX_MV_REF_CANDIDATES];
  int_mv best_mv, best_second_mv;

  int mb_mode_context[MAX_REF_FRAMES];

  unsigned char mb_skip_coeff;                                /* does this mb has coefficients at all, 1=no coefficients, 0=need decode tokens */
  unsigned char segment_id;           // Segment id for current frame

  // Flags used for prediction status of various bistream signals
  unsigned char seg_id_predicted;

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

#define VP9_REF_SCALE_SHIFT 14
#define VP9_REF_NO_SCALE 16384

struct scale_factors {
  int x_scale_fp;   // horizontal fixed point scale factor
  int y_scale_fp;   // vertical fixed point scale factor
  int x_offset_q4;
  int x_step_q4;
  int y_offset_q4;
  int y_step_q4;

  int (*scale_value_x)(int val, const struct scale_factors *scale);
  int (*scale_value_y)(int val, const struct scale_factors *scale);
  void (*set_scaled_offsets)(struct scale_factors *scale, int row, int col);
  int_mv32 (*scale_mv_q3_to_q4)(const int_mv *src_mv,
                                const struct scale_factors *scale);
  int32_t (*scale_mv_component_q4)(int mv_q4, int scale_fp, int offset_q4);

  convolve_fn_t predict[2][2][2];  // horiz, vert, avg
};

#if CONFIG_ALPHA
enum { MAX_MB_PLANE = 4 };
#else
enum { MAX_MB_PLANE = 3 };
#endif

struct buf_2d {
  uint8_t *buf;
  int stride;
};

struct macroblockd_plane {
  DECLARE_ALIGNED(16, int16_t,  qcoeff[64 * 64]);
  DECLARE_ALIGNED(16, int16_t,  dqcoeff[64 * 64]);
  DECLARE_ALIGNED(16, uint16_t, eobs[256]);
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

  /* 0 (do not update) 1 (update) the macroblock segmentation feature data. */
  unsigned char update_mb_segmentation_data;

  /* 0 (do not update) 1 (update) the macroblock segmentation feature data. */
  unsigned char mb_segment_abs_delta;

  /* Per frame flags that define which MB level features (such as quantizer or loop filter level) */
  /* are enabled and when enabled the proabilities used to decode the per MB flags in MB_MODE_INFO */

  // Probability Tree used to code Segment number
  vp9_prob mb_segment_tree_probs[MB_SEG_TREE_PROBS];

  // Segment features
  int16_t segment_feature_data[MAX_MB_SEGMENTS][SEG_LVL_MAX];
  unsigned int segment_feature_mask[MAX_MB_SEGMENTS];

  /* mode_based Loop filter adjustment */
  unsigned char mode_ref_lf_delta_enabled;
  unsigned char mode_ref_lf_delta_update;

  /* Delta values have the range +/- MAX_LOOP_FILTER */
  /* 0 = Intra, Last, GF, ARF */
  signed char last_ref_lf_deltas[MAX_REF_LF_DELTAS];
  /* 0 = Intra, Last, GF, ARF */
  signed char ref_lf_deltas[MAX_REF_LF_DELTAS];
  /* 0 = ZERO_MV, MV */
  signed char last_mode_lf_deltas[MAX_MODE_LF_DELTAS];
  /* 0 = ZERO_MV, MV */
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
  void (*inv_txm4x4_1_add)(int16_t *input, uint8_t *dest, int stride);
  void (*inv_txm4x4_add)(int16_t *input, uint8_t *dest, int stride);
  void (*itxm_add)(int16_t *input, uint8_t *dest, int stride, int eob);

  struct subpix_fn_table  subpix;

  int allow_high_precision_mv;

  int corrupted;

  int sb_index;   // index of 32x32 block inside the 64x64 block
  int mb_index;   // index of 16x16 block inside the 32x32 block
  int b_index;    // index of 8x8 block inside the 16x16 block
  int ab_index;   // index of 4x4 block inside the 8x8 block
  int q_index;

} MACROBLOCKD;

static int *get_sb_index(MACROBLOCKD *xd, BLOCK_SIZE_TYPE subsize) {
  switch (subsize) {
    case BLOCK_SIZE_SB64X64:
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
    case BLOCK_SIZE_SB8X4:
    case BLOCK_SIZE_SB4X8:
    case BLOCK_SIZE_AB4X4:
      return &xd->ab_index;
    default:
      assert(0);
      return NULL;
  }
}

static INLINE void update_partition_context(MACROBLOCKD *xd,
                                            BLOCK_SIZE_TYPE sb_type,
                                            BLOCK_SIZE_TYPE sb_size) {
  int bsl = b_width_log2(sb_size), bs = (1 << bsl) / 2;
  int bwl = b_width_log2(sb_type);
  int bhl = b_height_log2(sb_type);
  int boffset = b_width_log2(BLOCK_SIZE_SB64X64) - bsl;
  int i;

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

  return (left * 2 + above) + bsl * PARTITION_PLOFFSET;
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
      else if (bsize == BLOCK_SIZE_SB8X8)
        subsize = BLOCK_SIZE_SB8X4;
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
      else if (bsize == BLOCK_SIZE_SB8X8)
        subsize = BLOCK_SIZE_SB4X8;
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
      else if (bsize == BLOCK_SIZE_SB8X8)
        subsize = BLOCK_SIZE_AB4X4;
      else
        assert(0);
      break;
    default:
      assert(0);
  }
  return subsize;
}

// transform mapping
static TX_TYPE txfm_map(MB_PREDICTION_MODE bmode) {
  switch (bmode) {
    case TM_PRED :
    case D135_PRED :
      return ADST_ADST;

    case V_PRED :
    case D117_PRED :
    case D63_PRED:
      return ADST_DCT;

    case H_PRED :
    case D153_PRED :
    case D27_PRED :
      return DCT_ADST;

    default:
      return DCT_DCT;
  }
}

static TX_TYPE get_tx_type_4x4(const MACROBLOCKD *xd, int ib) {
  TX_TYPE tx_type;
  MODE_INFO *mi = xd->mode_info_context;
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  if (xd->lossless || mbmi->ref_frame[0] != INTRA_FRAME)
    return DCT_DCT;
  if (mbmi->sb_type < BLOCK_SIZE_SB8X8) {
    tx_type = txfm_map(mi->bmi[ib].as_mode.first);
  } else {
    assert(mbmi->mode <= TM_PRED);
    tx_type = txfm_map(mbmi->mode);
  }
  return tx_type;
}

static TX_TYPE get_tx_type_8x8(const MACROBLOCKD *xd, int ib) {
  TX_TYPE tx_type = DCT_DCT;
  if (xd->mode_info_context->mbmi.mode <= TM_PRED) {
    tx_type = txfm_map(xd->mode_info_context->mbmi.mode);
  }
  return tx_type;
}

static TX_TYPE get_tx_type_16x16(const MACROBLOCKD *xd, int ib) {
  TX_TYPE tx_type = DCT_DCT;
  if (xd->mode_info_context->mbmi.mode <= TM_PRED) {
    tx_type = txfm_map(xd->mode_info_context->mbmi.mode);
  }
  return tx_type;
}

void vp9_setup_block_dptrs(MACROBLOCKD *xd,
                           int subsampling_x, int subsampling_y);

static TX_SIZE get_uv_tx_size(const MB_MODE_INFO *mbmi) {
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

static INLINE int plane_block_width(BLOCK_SIZE_TYPE bsize,
                                    const struct macroblockd_plane* plane) {
  return 4 << (b_width_log2(bsize) - plane->subsampling_x);
}

static INLINE int plane_block_height(BLOCK_SIZE_TYPE bsize,
                                     const struct macroblockd_plane* plane) {
  return 4 << (b_height_log2(bsize) - plane->subsampling_y);
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
  const MB_MODE_INFO* mbmi = &xd->mode_info_context->mbmi;
  const TX_SIZE tx_size = plane ? get_uv_tx_size(mbmi)
                                : mbmi->txfm_size;
  const int block_size_b = bw + bh;
  const int txfrm_size_b = tx_size * 2;

  // subsampled size of the block
  const int ss_sum = xd->plane[plane].subsampling_x
      + xd->plane[plane].subsampling_y;
  const int ss_block_size = block_size_b - ss_sum;

  const int step = 1 << txfrm_size_b;

  int i;

  assert(txfrm_size_b <= block_size_b);
  assert(txfrm_size_b <= ss_block_size);

  // If mb_to_right_edge is < 0 we are in a situation in which
  // the current block size extends into the UMV and we won't
  // visit the sub blocks that are wholly within the UMV.
  if (xd->mb_to_right_edge < 0 || xd->mb_to_bottom_edge < 0) {
    int r, c;
    const int sw = bw - xd->plane[plane].subsampling_x;
    const int sh = bh - xd->plane[plane].subsampling_y;
    int max_blocks_wide = 1 << sw;
    int max_blocks_high = 1 << sh;

    // xd->mb_to_right_edge is in units of pixels * 8.  This converts
    // it to 4x4 block sizes.
    if (xd->mb_to_right_edge < 0)
      max_blocks_wide +=
          + (xd->mb_to_right_edge >> (5 + xd->plane[plane].subsampling_x));

    if (xd->mb_to_bottom_edge < 0)
      max_blocks_high +=
          + (xd->mb_to_bottom_edge >> (5 + xd->plane[plane].subsampling_y));

    i = 0;
    // Unlike the normal case - in here we have to keep track of the
    // row and column of the blocks we use so that we know if we are in
    // the unrestricted motion border..
    for (r = 0; r < (1 << sh); r += (1 << tx_size)) {
      for (c = 0; c < (1 << sw); c += (1 << tx_size)) {
        if (r < max_blocks_high && c < max_blocks_wide)
          visit(plane, i, bsize, txfrm_size_b, arg);
        i += step;
      }
    }
  } else {
    for (i = 0; i < (1 << ss_block_size); i += step) {
      visit(plane, i, bsize, txfrm_size_b, arg);
    }
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

  // block sizes in number of 4x4 blocks log 2 ("*_b")
  // 4x4=0, 8x8=2, 16x16=4, 32x32=6, 64x64=8
  // subsampled size of the block
  const int bwl = b_width_log2(bsize) - xd->plane[plane].subsampling_x;
  const int bhl = b_height_log2(bsize) - xd->plane[plane].subsampling_y;

  // size of the predictor to use.
  int pred_w, pred_h;

  if (xd->mode_info_context->mbmi.sb_type < BLOCK_SIZE_SB8X8) {
    assert(bsize == BLOCK_SIZE_SB8X8);
    pred_w = 0;
    pred_h = 0;
  } else {
    pred_w = bwl;
    pred_h = bhl;
  }
  assert(pred_w <= bwl);
  assert(pred_h <= bhl);

  // visit each subblock in raster order
  i = 0;
  for (y = 0; y < 1 << bhl; y += 1 << pred_h) {
    for (x = 0; x < 1 << bwl; x += 1 << pred_w) {
      visit(plane, i, bsize, pred_w, pred_h, arg);
      i += 1 << pred_w;
    }
    i += (1 << (bwl + pred_h)) - (1 << bwl);
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
  const int stride = plane_block_width(bsize, &xd->plane[plane]);
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

static void extend_for_intra(MACROBLOCKD* const xd, int plane, int block,
                             BLOCK_SIZE_TYPE bsize, int ss_txfrm_size) {
  const int bw = plane_block_width(bsize, &xd->plane[plane]);
  const int bh = plane_block_height(bsize, &xd->plane[plane]);
  int x, y;
  txfrm_block_to_raster_xy(xd, bsize, plane, block, ss_txfrm_size, &x, &y);
  x = x * 4 - 1;
  y = y * 4 - 1;
  // Copy a pixel into the umv if we are in a situation where the block size
  // extends into the UMV.
  // TODO(JBB): Should be able to do the full extend in place so we don't have
  // to do this multiple times.
  if (xd->mb_to_right_edge < 0) {
    int umv_border_start = bw
        + (xd->mb_to_right_edge >> (3 + xd->plane[plane].subsampling_x));

    if (x + bw > umv_border_start)
      vpx_memset(
          xd->plane[plane].dst.buf + y * xd->plane[plane].dst.stride
              + umv_border_start,
          *(xd->plane[plane].dst.buf + y * xd->plane[plane].dst.stride
              + umv_border_start - 1),
          bw);
  }
  if (xd->mb_to_bottom_edge < 0) {
    int umv_border_start = bh
        + (xd->mb_to_bottom_edge >> (3 + xd->plane[plane].subsampling_y));
    int i;
    uint8_t c = *(xd->plane[plane].dst.buf
        + (umv_border_start - 1) * xd->plane[plane].dst.stride + x);

    uint8_t *d = xd->plane[plane].dst.buf
        + umv_border_start * xd->plane[plane].dst.stride + x;

    if (y + bh > umv_border_start)
      for (i = 0; i < bh; i++, d += xd->plane[plane].dst.stride)
        *d = c;
  }
}
static void set_contexts_on_border(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize,
                                   int plane, int ss_tx_size, int eob, int aoff,
                                   int loff, ENTROPY_CONTEXT *A,
                                   ENTROPY_CONTEXT *L) {
  const int bw = b_width_log2(bsize), bh = b_height_log2(bsize);
  const int sw = bw - xd->plane[plane].subsampling_x;
  const int sh = bh - xd->plane[plane].subsampling_y;
  int mi_blocks_wide = 1 << sw;
  int mi_blocks_high = 1 << sh;
  int tx_size_in_blocks = (1 << ss_tx_size);
  int above_contexts = tx_size_in_blocks;
  int left_contexts = tx_size_in_blocks;
  int pt;

  // xd->mb_to_right_edge is in units of pixels * 8.  This converts
  // it to 4x4 block sizes.
  if (xd->mb_to_right_edge < 0) {
    mi_blocks_wide += (xd->mb_to_right_edge
        >> (5 + xd->plane[plane].subsampling_x));
  }

  // this code attempts to avoid copying into contexts that are outside
  // our border.  Any blocks that do are set to 0...
  if (above_contexts + aoff > mi_blocks_wide)
    above_contexts = mi_blocks_wide - aoff;

  if (xd->mb_to_bottom_edge < 0) {
    mi_blocks_high += (xd->mb_to_bottom_edge
        >> (5 + xd->plane[plane].subsampling_y));
  }
  if (left_contexts + loff > mi_blocks_high) {
    left_contexts = mi_blocks_high - loff;
  }

  for (pt = 0; pt < above_contexts; pt++)
    A[pt] = eob > 0;
  for (pt = above_contexts; pt < (1 << ss_tx_size); pt++)
    A[pt] = 0;
  for (pt = 0; pt < left_contexts; pt++)
    L[pt] = eob > 0;
  for (pt = left_contexts; pt < (1 << ss_tx_size); pt++)
    L[pt] = 0;
}


#endif  // VP9_COMMON_VP9_BLOCKD_H_
