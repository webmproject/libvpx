/*
 *  Copyright (c) 2019 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VP9_ENCODER_VP9_NON_GREEDY_MV_H_
#define VPX_VP9_ENCODER_VP9_NON_GREEDY_MV_H_

#include "vp9/common/vp9_enums.h"
#include "vp9/common/vp9_blockd.h"
#include "vpx_scale/yv12config.h"
#include "vpx_dsp/variance.h"

#ifdef __cplusplus
extern "C" {
#endif
#define NB_MVS_NUM 4
#define LOG2_PRECISION 20

int64_t vp9_nb_mvs_inconsistency(const MV *mv, const int_mv *nb_full_mvs,
                                 int mv_num);

void vp9_get_smooth_motion_field(const MV *search_mf, const int (*M)[4],
                                 int rows, int cols, float alpha, int num_iters,
                                 MV *smooth_mf);

void vp9_get_local_structure(const YV12_BUFFER_CONFIG *ref_frame,
                             const vp9_variance_fn_ptr_t *fn_ptr, int mi_rows,
                             int mi_cols, BLOCK_SIZE bsize, int (*M)[4]);

#ifdef __cplusplus
}  // extern "C"
#endif
#endif  // VPX_VP9_ENCODER_VP9_NON_GREEDY_MV_H_
