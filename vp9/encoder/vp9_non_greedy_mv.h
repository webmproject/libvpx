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

#ifdef __cplusplus
extern "C" {
#endif

#define NB_MVS_NUM 4
#define LOG2_PRECISION 20

int64_t vp9_nb_mvs_inconsistency(const MV *mv, const int_mv *nb_full_mvs,
                                 int mv_num);

#ifdef __cplusplus
}  // extern "C"
#endif
#endif  // VPX_VP9_ENCODER_VP9_NON_GREEDY_MV_H_
