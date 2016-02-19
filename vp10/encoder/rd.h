/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_RD_H_
#define VP10_ENCODER_RD_H_

#include <limits.h>

#include "vp10/common/blockd.h"

#include "vp10/encoder/block.h"
#include "vp10/encoder/context_tree.h"
#include "vp10/encoder/cost.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RDDIV_BITS          7
#define RD_EPB_SHIFT        6

#define RDCOST(RM, DM, R, D) \
  (ROUND_POWER_OF_TWO(((int64_t)R) * (RM), VP9_PROB_COST_SHIFT) + (D << DM))

#define RDCOST_DBL(RM, DM, R, D)                                   \
  (((((double)(R)) * (RM)) / (double)(1 << VP9_PROB_COST_SHIFT)) + \
   ((double)(D) * (1 << (DM))))

#define QIDX_SKIP_THRESH     115

#define MV_COST_WEIGHT      108
#define MV_COST_WEIGHT_SUB  120

#define INVALID_MV 0x80008000

#if CONFIG_EXT_REFS
#if CONFIG_EXT_INTER
#define MAX_MODES 85
#else
#define MAX_MODES 54
#endif  // CONFIG_EXT_INTER
#else
#if CONFIG_EXT_INTER
#define MAX_MODES 43
#else
#define MAX_MODES 30
#endif  // CONFIG_EXT_INTER
#endif  // CONFIG_EXT_REFS

#if CONFIG_EXT_REFS
#define MAX_REFS  12
#else
#define MAX_REFS  6
#endif  // CONFIG_EXT_REFS

#define RD_THRESH_MAX_FACT 64
#define RD_THRESH_INC      1

// This enumerator type needs to be kept aligned with the mode order in
// const MODE_DEFINITION vp10_mode_order[MAX_MODES] used in the rd code.
typedef enum {
  THR_NEARESTMV,
#if CONFIG_EXT_REFS
  THR_NEARESTL2,
  THR_NEARESTL3,
  THR_NEARESTL4,
#endif  // CONFIG_EXT_REFS
  THR_NEARESTA,
  THR_NEARESTG,

  THR_DC,

  THR_NEWMV,
#if CONFIG_EXT_REFS
  THR_NEWL2,
  THR_NEWL3,
  THR_NEWL4,
#endif  // CONFIG_EXT_REFS
  THR_NEWA,
  THR_NEWG,

  THR_NEARMV,
#if CONFIG_EXT_REFS
  THR_NEARL2,
  THR_NEARL3,
  THR_NEARL4,
#endif  // CONFIG_EXT_REFS
  THR_NEARA,
  THR_NEARG,

#if CONFIG_EXT_INTER
  THR_NEWFROMNEARMV,
#if CONFIG_EXT_REFS
  THR_NEWFROMNEARL2,
  THR_NEWFROMNEARL3,
  THR_NEWFROMNEARL4,
#endif  // CONFIG_EXT_REFS
  THR_NEWFROMNEARA,
  THR_NEWFROMNEARG,
#endif  // CONFIG_EXT_INTER

  THR_ZEROMV,
#if CONFIG_EXT_REFS
  THR_ZEROL2,
  THR_ZEROL3,
  THR_ZEROL4,
#endif  // CONFIG_EXT_REFS
  THR_ZEROG,
  THR_ZEROA,

#if CONFIG_EXT_INTER
  THR_COMP_NEAREST_NEARESTLA,
#if CONFIG_EXT_REFS
  THR_COMP_NEAREST_NEARESTL2A,
  THR_COMP_NEAREST_NEARESTL3A,
  THR_COMP_NEAREST_NEARESTL4A,
#endif  // CONFIG_EXT_REFS
  THR_COMP_NEAREST_NEARESTGA,
#else  // CONFIG_EXT_INTER
  THR_COMP_NEARESTLA,
#if CONFIG_EXT_REFS
  THR_COMP_NEARESTL2A,
  THR_COMP_NEARESTL3A,
  THR_COMP_NEARESTL4A,
#endif  // CONFIG_EXT_REFS
  THR_COMP_NEARESTGA,
#endif  // CONFIG_EXT_INTER

  THR_TM,

#if CONFIG_EXT_INTER
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

#if CONFIG_EXT_REFS
  THR_COMP_NEAR_NEARESTL2A,
  THR_COMP_NEAREST_NEARL2A,
  THR_COMP_NEW_NEARESTL2A,
  THR_COMP_NEAREST_NEWL2A,
  THR_COMP_NEW_NEARL2A,
  THR_COMP_NEAR_NEWL2A,
  THR_COMP_NEW_NEWL2A,
  THR_COMP_ZERO_ZEROL2A,

  THR_COMP_NEAR_NEARESTL3A,
  THR_COMP_NEAREST_NEARL3A,
  THR_COMP_NEW_NEARESTL3A,
  THR_COMP_NEAREST_NEWL3A,
  THR_COMP_NEW_NEARL3A,
  THR_COMP_NEAR_NEWL3A,
  THR_COMP_NEW_NEWL3A,
  THR_COMP_ZERO_ZEROL3A,

  THR_COMP_NEAR_NEARESTL4A,
  THR_COMP_NEAREST_NEARL4A,
  THR_COMP_NEW_NEARESTL4A,
  THR_COMP_NEAREST_NEWL4A,
  THR_COMP_NEW_NEARL4A,
  THR_COMP_NEAR_NEWL4A,
  THR_COMP_NEW_NEWL4A,
  THR_COMP_ZERO_ZEROL4A,
#endif  // CONFIG_EXT_REFS
#else
  THR_COMP_NEARLA,
  THR_COMP_NEWLA,
#if CONFIG_EXT_REFS
  THR_COMP_NEARL2A,
  THR_COMP_NEWL2A,
  THR_COMP_NEARL3A,
  THR_COMP_NEWL3A,
  THR_COMP_NEARL4A,
  THR_COMP_NEWL4A,
#endif  // CONFIG_EXT_REFS
  THR_COMP_NEARGA,
  THR_COMP_NEWGA,

  THR_COMP_ZEROLA,
#if CONFIG_EXT_REFS
  THR_COMP_ZEROL2A,
  THR_COMP_ZEROL3A,
  THR_COMP_ZEROL4A,
#endif  // CONFIG_EXT_REFS
  THR_COMP_ZEROGA,
#endif  // CONFIG_EXT_INTER

  THR_H_PRED,
  THR_V_PRED,
  THR_D135_PRED,
  THR_D207_PRED,
  THR_D153_PRED,
  THR_D63_PRED,
  THR_D117_PRED,
  THR_D45_PRED,
} THR_MODES;

typedef enum {
  THR_LAST,
#if CONFIG_EXT_REFS
  THR_LAST2,
  THR_LAST3,
  THR_LAST4,
#endif  // CONFIG_EXT_REFS
  THR_GOLD,
  THR_ALTR,
  THR_COMP_LA,
#if CONFIG_EXT_REFS
  THR_COMP_L2A,
  THR_COMP_L3A,
  THR_COMP_L4A,
#endif  // CONFIG_EXT_REFS
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

  int64_t prediction_type_threshes[MAX_REF_FRAMES][REFERENCE_MODES];

  int64_t filter_threshes[MAX_REF_FRAMES][SWITCHABLE_FILTER_CONTEXTS];

  int RDMULT;
  int RDDIV;
} RD_OPT;

typedef struct RD_COST {
  int rate;
  int64_t dist;
  int64_t rdcost;
} RD_COST;

// Reset the rate distortion cost values to maximum (invalid) value.
void vp10_rd_cost_reset(RD_COST *rd_cost);
// Initialize the rate distortion cost values to zero.
void vp10_rd_cost_init(RD_COST *rd_cost);

struct TileInfo;
struct TileDataEnc;
struct VP10_COMP;
struct macroblock;

int vp10_compute_rd_mult(const struct VP10_COMP *cpi, int qindex);

void vp10_initialize_rd_consts(struct VP10_COMP *cpi);

void vp10_initialize_me_consts(struct VP10_COMP *cpi,
                               MACROBLOCK *x, int qindex);

void vp10_model_rd_from_var_lapndz(unsigned int var, unsigned int n,
                                  unsigned int qstep, int *rate,
                                  int64_t *dist);

int vp10_get_switchable_rate(const struct VP10_COMP *cpi,
                            const MACROBLOCKD *const xd);

int vp10_raster_block_offset(BLOCK_SIZE plane_bsize,
                            int raster_block, int stride);

int16_t* vp10_raster_block_offset_int16(BLOCK_SIZE plane_bsize,
                                       int raster_block, int16_t *base);

YV12_BUFFER_CONFIG *vp10_get_scaled_ref_frame(const struct VP10_COMP *cpi,
                                             int ref_frame);

void vp10_init_me_luts(void);

#if CONFIG_REF_MV
void vp10_set_mvcost(MACROBLOCK *x, MV_REFERENCE_FRAME ref_frame);
#endif

void vp10_get_entropy_contexts(BLOCK_SIZE bsize, TX_SIZE tx_size,
                              const struct macroblockd_plane *pd,
                              ENTROPY_CONTEXT t_above[16],
                              ENTROPY_CONTEXT t_left[16]);

void vp10_set_rd_speed_thresholds(struct VP10_COMP *cpi);

void vp10_set_rd_speed_thresholds_sub8x8(struct VP10_COMP *cpi);

void vp10_update_rd_thresh_fact(int (*fact)[MAX_MODES], int rd_thresh,
                               int bsize, int best_mode_index);

static INLINE int rd_less_than_thresh(int64_t best_rd, int thresh,
                                      int thresh_fact) {
    return best_rd < ((int64_t)thresh * thresh_fact >> 5) || thresh == INT_MAX;
}

void vp10_mv_pred(struct VP10_COMP *cpi, MACROBLOCK *x,
                 uint8_t *ref_y_buffer, int ref_y_stride,
                 int ref_frame, BLOCK_SIZE block_size);

static INLINE void set_error_per_bit(MACROBLOCK *x, int rdmult) {
  x->errorperbit = rdmult >> RD_EPB_SHIFT;
  x->errorperbit += (x->errorperbit == 0);
}

void vp10_setup_pred_block(const MACROBLOCKD *xd,
                          struct buf_2d dst[MAX_MB_PLANE],
                          const YV12_BUFFER_CONFIG *src,
                          int mi_row, int mi_col,
                          const struct scale_factors *scale,
                          const struct scale_factors *scale_uv);

int vp10_get_intra_cost_penalty(int qindex, int qdelta,
                               vpx_bit_depth_t bit_depth);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_RD_H_
