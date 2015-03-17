/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_RD_H_
#define VP9_ENCODER_VP9_RD_H_

#include <limits.h>

#include "vp9/common/vp9_blockd.h"

#include "vp9/encoder/vp9_block.h"
#include "vp9/encoder/vp9_context_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RDDIV_BITS          7

#define RDCOST(RM, DM, R, D) \
  (((128 + ((int64_t)R) * (RM)) >> 8) + (D << DM))
#define QIDX_SKIP_THRESH     115

#define MV_COST_WEIGHT      108
#define MV_COST_WEIGHT_SUB  120

#define INVALID_MV 0x80008000

#if CONFIG_COMPOUND_MODES

#if CONFIG_INTERINTRA

#if CONFIG_NEWMVREF
#define MAX_MODES 55
#define INTERINTRA_START_MODE 43
#else
#define MAX_MODES 52
#define INTERINTRA_START_MODE 40
#endif  // CONFIG_NEWMVREF

#else  // CONFIG_INTERINTRA

#if CONFIG_NEWMVREF
#define MAX_MODES 43
#else
#define MAX_MODES 40
#endif  // CONFIG_NEWMVREF

#endif  // CONFIG_INTERINTRA

#else   // CONFIG_COMPOUND_MODES

#if CONFIG_INTERINTRA

#if CONFIG_NEWMVREF
#define MAX_MODES 47
#define INTERINTRA_START_MODE 35
#else
#define MAX_MODES 42
#define INTERINTRA_START_MODE 30
#endif  // CONFIG_NEWMVREF

#else  // CONFIG_INTERINTRA

#if CONFIG_NEWMVREF
#define MAX_MODES 35
#else
#define MAX_MODES 30
#endif  // CONFIG_NEWMVREF

#endif  // CONFIG_INTERINTRA

#endif  // CONFIG_COMPOUND_MODES

#define MAX_REFS  6

// This enumerator type needs to be kept aligned with the mode order in
// const MODE_DEFINITION vp9_mode_order[MAX_MODES] used in the rd code.
typedef enum {
  THR_NEARESTMV,
  THR_NEARESTA,
  THR_NEARESTG,

  THR_DC,

  THR_NEWMV,
  THR_NEWA,
  THR_NEWG,

#if CONFIG_NEWMVREF
  THR_NEAR_FORNEWMV,
  THR_NEAR_FORNEWA,
  THR_NEAR_FORNEWG,
#endif  // CONFIG_NEWMVREF

  THR_NEARMV,
  THR_NEARA,
  THR_NEARG,

  THR_ZEROMV,
  THR_ZEROG,
  THR_ZEROA,

#if CONFIG_COMPOUND_MODES
  THR_COMP_NEAREST_NEARESTLA,
  THR_COMP_NEAREST_NEARESTGA,
#else
  THR_COMP_NEARESTLA,
  THR_COMP_NEARESTGA,
#endif

  THR_TM,

#if CONFIG_COMPOUND_MODES
  THR_COMP_NEAR_NEARESTLA,
  THR_COMP_NEAR_NEARESTGA,
  THR_COMP_NEAREST_NEARLA,
  THR_COMP_NEAREST_NEARGA,

  THR_COMP_NEW_NEARESTLA,
  THR_COMP_NEW_NEARESTGA,
  THR_COMP_NEAREST_NEWLA,
  THR_COMP_NEAREST_NEWGA,

  THR_COMP_NEW_NEARLA,
  THR_COMP_NEW_NEARGA,
  THR_COMP_NEAR_NEWLA,
  THR_COMP_NEAR_NEWGA,

  THR_COMP_NEW_NEWLA,
  THR_COMP_NEW_NEWGA,

  THR_COMP_ZERO_ZEROLA,
  THR_COMP_ZERO_ZEROGA,
#else
  THR_COMP_NEARLA,
  THR_COMP_NEWLA,
  THR_COMP_NEARGA,
  THR_COMP_NEWGA,

#if CONFIG_NEWMVREF
  THR_COMP_NEAR_FORNEWLA,
  THR_COMP_NEAR_FORNEWGA,
#endif  // CONFIG_NEWMVREF

  THR_COMP_ZEROLA,
  THR_COMP_ZEROGA,
#endif  // CONFIG_COMPOUND_MODES

  THR_H_PRED,
  THR_V_PRED,
  THR_D135_PRED,
  THR_D207_PRED,
  THR_D153_PRED,
  THR_D63_PRED,
  THR_D117_PRED,
  THR_D45_PRED,

#if CONFIG_INTERINTRA
  THR_COMP_INTERINTRA_ZEROL,
  THR_COMP_INTERINTRA_NEARESTL,
  THR_COMP_INTERINTRA_NEARL,
  THR_COMP_INTERINTRA_NEWL,

  THR_COMP_INTERINTRA_ZEROG,
  THR_COMP_INTERINTRA_NEARESTG,
  THR_COMP_INTERINTRA_NEARG,
  THR_COMP_INTERINTRA_NEWG,

  THR_COMP_INTERINTRA_ZEROA,
  THR_COMP_INTERINTRA_NEARESTA,
  THR_COMP_INTERINTRA_NEARA,
  THR_COMP_INTERINTRA_NEWA,
#endif  // CONFIG_INTERINTRA
} THR_MODES;

typedef enum {
  THR_LAST,
  THR_GOLD,
  THR_ALTR,
  THR_COMP_LA,
  THR_COMP_GA,
  THR_INTRA,
} THR_MODES_SUB8X8;

typedef struct RD_OPT {
  // Thresh_mult is used to set a threshold for the rd score. A higher value
  // means that we will accept the best mode so far more often. This number
  // is used in combination with the current block size, and thresh_freq_fact
  // to pick a threshold.
  int thresh_mult[MAX_MODES];
  int thresh_mult_sub8x8[MAX_REFS];

  int threshes[MAX_SEGMENTS][BLOCK_SIZES][MAX_MODES];
  int thresh_freq_fact[BLOCK_SIZES][MAX_MODES];

  int mode_map[BLOCK_SIZES][MAX_MODES];

  int64_t comp_pred_diff[REFERENCE_MODES];
  int64_t prediction_type_threshes[MAX_REF_FRAMES][REFERENCE_MODES];
  int64_t tx_select_diff[TX_MODES];
  // TODO(agrange): can this overflow?
  int tx_select_threshes[MAX_REF_FRAMES][TX_MODES];

  int64_t filter_diff[SWITCHABLE_FILTER_CONTEXTS];
  int64_t filter_threshes[MAX_REF_FRAMES][SWITCHABLE_FILTER_CONTEXTS];
  int64_t filter_cache[SWITCHABLE_FILTER_CONTEXTS];
  int64_t mask_filter;

  int RDMULT;
  int RDDIV;
} RD_OPT;

typedef struct RD_COST {
  int rate;
  int64_t dist;
  int64_t rdcost;
} RD_COST;

// Reset the rate distortion cost values to maximum (invalid) value.
void vp9_rd_cost_reset(RD_COST *rd_cost);
// Initialize the rate distortion cost values to zero.
void vp9_rd_cost_init(RD_COST *rd_cost);

struct TileInfo;
struct VP9_COMP;
struct macroblock;

int vp9_compute_rd_mult(const struct VP9_COMP *cpi, int qindex);

void vp9_initialize_rd_consts(struct VP9_COMP *cpi);

void vp9_initialize_me_consts(struct VP9_COMP *cpi, int qindex);

void vp9_model_rd_from_var_lapndz(unsigned int var, unsigned int n,
                                  unsigned int qstep, int *rate,
                                  int64_t *dist);

int vp9_get_switchable_rate(const struct VP9_COMP *cpi);

const YV12_BUFFER_CONFIG *vp9_get_scaled_ref_frame(const struct VP9_COMP *cpi,
                                                   int ref_frame);

void vp9_init_me_luts();

void vp9_get_entropy_contexts(BLOCK_SIZE bsize, TX_SIZE tx_size,
                              const struct macroblockd_plane *pd,
                              ENTROPY_CONTEXT t_above[16],
                              ENTROPY_CONTEXT t_left[16]);

void vp9_set_rd_speed_thresholds(struct VP9_COMP *cpi);

void vp9_set_rd_speed_thresholds_sub8x8(struct VP9_COMP *cpi);

static INLINE int rd_less_than_thresh(int64_t best_rd, int thresh,
                                      int thresh_fact) {
    return best_rd < ((int64_t)thresh * thresh_fact >> 5) || thresh == INT_MAX;
}

void vp9_mv_pred(struct VP9_COMP *cpi, MACROBLOCK *x,
                 uint8_t *ref_y_buffer, int ref_y_stride,
                 int ref_frame, BLOCK_SIZE block_size);

void vp9_setup_pred_block(const MACROBLOCKD *xd,
                          struct buf_2d dst[MAX_MB_PLANE],
                          const YV12_BUFFER_CONFIG *src,
                          int mi_row, int mi_col,
                          const struct scale_factors *scale,
                          const struct scale_factors *scale_uv);

int vp9_get_intra_cost_penalty(int qindex, int qdelta,
                               vpx_bit_depth_t bit_depth);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_RD_H_
