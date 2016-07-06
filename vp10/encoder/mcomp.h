/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP10_ENCODER_MCOMP_H_
#define VP10_ENCODER_MCOMP_H_

#include "vp10/encoder/block.h"
#include "vpx_dsp/variance.h"

#ifdef __cplusplus
extern "C" {
#endif

// The maximum number of steps in a step search given the largest
// allowed initial step
#define MAX_MVSEARCH_STEPS 11
// Max full pel mv specified in the unit of full pixel
// Enable the use of motion vector in range [-1023, 1023].
#define MAX_FULL_PEL_VAL ((1 << (MAX_MVSEARCH_STEPS - 1)) - 1)
// Maximum size of the first step in full pel units
#define MAX_FIRST_STEP (1 << (MAX_MVSEARCH_STEPS-1))
// Allowed motion vector pixel distance outside image border
// for Block_16x16
#define BORDER_MV_PIXELS_B16 (16 + VP9_INTERP_EXTEND)

// motion search site
typedef struct search_site {
  MV mv;
  int offset;
} search_site;

typedef struct search_site_config {
  search_site ss[8 * MAX_MVSEARCH_STEPS + 1];
  int ss_count;
  int searches_per_step;
} search_site_config;

void vp10_init_dsmotion_compensation(search_site_config *cfg, int stride);
void vp10_init3smotion_compensation(search_site_config *cfg,  int stride);

void vp10_set_mv_search_range(MACROBLOCK *x, const MV *mv);
int vp10_mv_bit_cost(const MV *mv, const MV *ref,
                    const int *mvjcost, int *mvcost[2], int weight);

// Utility to compute variance + MV rate cost for a given MV
int vp10_get_mvpred_var(const MACROBLOCK *x,
                       const MV *best_mv, const MV *center_mv,
                       const vp10_variance_fn_ptr_t *vfp,
                       int use_mvcost);
int vp10_get_mvpred_av_var(const MACROBLOCK *x,
                          const MV *best_mv, const MV *center_mv,
                          const uint8_t *second_pred,
                          const vp10_variance_fn_ptr_t *vfp,
                          int use_mvcost);

struct VP10_COMP;
struct SPEED_FEATURES;

int vp10_init_search_range(int size);

int vp10_refining_search_sad(const struct macroblock *x,
                            struct mv *ref_mv,
                            int sad_per_bit, int distance,
                            const vp10_variance_fn_ptr_t *fn_ptr,
                            const struct mv *center_mv);

// Runs sequence of diamond searches in smaller steps for RD.
int vp10_full_pixel_diamond(const struct VP10_COMP *cpi, MACROBLOCK *x,
                           MV *mvp_full, int step_param,
                           int sadpb, int further_steps, int do_refine,
                           int *cost_list,
                           const vp10_variance_fn_ptr_t *fn_ptr,
                           const MV *ref_mv, MV *dst_mv);

// Perform integral projection based motion estimation.
unsigned int vp10_int_pro_motion_estimation(const struct VP10_COMP *cpi,
                                            MACROBLOCK *x,
                                            BLOCK_SIZE bsize,
                                            int mi_row, int mi_col);


int vp10_hex_search(MACROBLOCK *x,
                    MV *start_mv,
                    int search_param,
                    int sad_per_bit,
                    int do_init_search,
                    int *cost_list,
                    const vp10_variance_fn_ptr_t *vfp,
                    int use_mvcost,
                    const MV *center_mv);

typedef int (fractional_mv_step_fp) (
    MACROBLOCK *x,
    const MV *ref_mv,
    int allow_hp,
    int error_per_bit,
    const vp10_variance_fn_ptr_t *vfp,
    int forced_stop,  // 0 - full, 1 - qtr only, 2 - half only
    int iters_per_step,
    int *cost_list,
    int *mvjcost, int *mvcost[2],
    int *distortion, unsigned int *sse1,
    const uint8_t *second_pred,
    int w, int h, int use_upsampled_ref);

extern fractional_mv_step_fp vp10_find_best_sub_pixel_tree;
extern fractional_mv_step_fp vp10_find_best_sub_pixel_tree_pruned;
extern fractional_mv_step_fp vp10_find_best_sub_pixel_tree_pruned_more;
extern fractional_mv_step_fp vp10_find_best_sub_pixel_tree_pruned_evenmore;

typedef int (*vp10_full_search_fn_t)(const MACROBLOCK *x,
                                     const MV *ref_mv, int sad_per_bit,
                                     int distance,
                                     const vp10_variance_fn_ptr_t *fn_ptr,
                                     const MV *center_mv, MV *best_mv);

typedef int (*vp10_diamond_search_fn_t)(const MACROBLOCK *x,
                                        const search_site_config *cfg,
                                        MV *ref_mv, MV *best_mv,
                                        int search_param, int sad_per_bit,
                                        int *num00,
                                        const vp10_variance_fn_ptr_t *fn_ptr,
                                        const MV *center_mv);

int vp10_refining_search_8p_c(MACROBLOCK *x,
                              int error_per_bit,
                              int search_range,
                              const vp10_variance_fn_ptr_t *fn_ptr,
                              const MV *center_mv, const uint8_t *second_pred);

struct VP10_COMP;

int vp10_full_pixel_search(struct VP10_COMP *cpi, MACROBLOCK *x,
                           BLOCK_SIZE bsize, MV *mvp_full,
                           int step_param, int error_per_bit,
                           int *cost_list, const MV *ref_mv,
                           int var_max, int rd);

#if CONFIG_EXT_INTER
int vp10_find_best_masked_sub_pixel_tree(const MACROBLOCK *x,
                                         const uint8_t *mask, int mask_stride,
                                         MV *bestmv, const MV *ref_mv,
                                         int allow_hp,
                                         int error_per_bit,
                                         const vp10_variance_fn_ptr_t *vfp,
                                         int forced_stop,
                                         int iters_per_step,
                                         int *mvjcost, int *mvcost[2],
                                         int *distortion,
                                         unsigned int *sse1,
                                         int is_second);
int vp10_find_best_masked_sub_pixel_tree_up(struct VP10_COMP *cpi,
                                            MACROBLOCK *x,
                                            const uint8_t *mask,
                                            int mask_stride,
                                            int mi_row, int mi_col,
                                            MV *bestmv, const MV *ref_mv,
                                            int allow_hp,
                                            int error_per_bit,
                                            const vp10_variance_fn_ptr_t *vfp,
                                            int forced_stop,
                                            int iters_per_step,
                                            int *mvjcost, int *mvcost[2],
                                            int *distortion,
                                            unsigned int *sse1,
                                            int is_second,
                                            int use_upsampled_ref);
int vp10_masked_full_pixel_diamond(const struct VP10_COMP *cpi, MACROBLOCK *x,
                                   const uint8_t *mask, int mask_stride,
                                   MV *mvp_full, int step_param,
                                   int sadpb, int further_steps, int do_refine,
                                   const vp10_variance_fn_ptr_t *fn_ptr,
                                   const MV *ref_mv, MV *dst_mv,
                                   int is_second);
#endif  // CONFIG_EXT_INTER

#if CONFIG_OBMC
int vp10_obmc_full_pixel_diamond(const struct VP10_COMP *cpi, MACROBLOCK *x,
                                 const int32_t *wsrc,
                                 const int32_t *mask,
                                 MV *mvp_full, int step_param,
                                 int sadpb, int further_steps, int do_refine,
                                 const vp10_variance_fn_ptr_t *fn_ptr,
                                 const MV *ref_mv, MV *dst_mv,
                                 int is_second);
int vp10_find_best_obmc_sub_pixel_tree_up(struct VP10_COMP *cpi, MACROBLOCK *x,
                                          const int32_t *wsrc,
                                          const int32_t *mask,
                                          int mi_row, int mi_col,
                                          MV *bestmv, const MV *ref_mv,
                                          int allow_hp, int error_per_bit,
                                          const vp10_variance_fn_ptr_t *vfp,
                                          int forced_stop, int iters_per_step,
                                          int *mvjcost, int *mvcost[2],
                                          int *distortion, unsigned int *sse1,
                                          int is_second,
                                          int use_upsampled_ref);
#endif  // CONFIG_OBMC
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_MCOMP_H_
