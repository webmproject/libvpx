/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_PRED_COMMON_H_
#define VP9_COMMON_VP9_PRED_COMMON_H_

#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_onyxc_int.h"

int vp9_get_segment_id(VP9_COMMON *cm, const uint8_t *segment_ids,
                       BLOCK_SIZE_TYPE bsize, int mi_row, int mi_col);


static INLINE unsigned char vp9_get_pred_context_seg_id(const VP9_COMMON *cm,
                                                        const MACROBLOCKD *xd) {
  int pred_context;
  const MODE_INFO * const mi = xd->mode_info_context;
  const MODE_INFO * const above_mi = mi - cm->mode_info_stride;
  const MODE_INFO * const left_mi = mi - 1;
  pred_context = above_mi->mbmi.seg_id_predicted;
  if (xd->left_available)
    pred_context += left_mi->mbmi.seg_id_predicted;
  return pred_context;
}
static INLINE vp9_prob vp9_get_pred_prob_seg_id(const VP9_COMMON *cm,
                                                const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_seg_id(cm, xd);
  return cm->segment_pred_probs[pred_context];
}
static INLINE unsigned char vp9_get_pred_flag_seg_id(
    const MACROBLOCKD * const xd) {
  return xd->mode_info_context->mbmi.seg_id_predicted;
}

void vp9_set_pred_flag_seg_id(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize,
                       unsigned char pred_flag);

static INLINE unsigned char vp9_get_pred_context_mbskip(const VP9_COMMON *cm,
                                                        const MACROBLOCKD *xd) {
  int pred_context;
  const MODE_INFO * const mi = xd->mode_info_context;
  const MODE_INFO * const above_mi = mi - cm->mode_info_stride;
  const MODE_INFO * const left_mi = mi - 1;
  pred_context = above_mi->mbmi.mb_skip_coeff;
  if (xd->left_available)
    pred_context += left_mi->mbmi.mb_skip_coeff;
  return pred_context;
}
static INLINE vp9_prob vp9_get_pred_prob_mbskip(const VP9_COMMON *cm,
                                                const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_mbskip(cm, xd);
  return cm->fc.mbskip_probs[pred_context];
}
static INLINE unsigned char vp9_get_pred_flag_mbskip(
    const MACROBLOCKD * const xd) {
  return xd->mode_info_context->mbmi.mb_skip_coeff;
}
void vp9_set_pred_flag_mbskip(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize,
                              unsigned char pred_flag);

unsigned char vp9_get_pred_context_switchable_interp(const VP9_COMMON *cm,
                                                     const MACROBLOCKD *xd);

static INLINE const vp9_prob *vp9_get_pred_probs_switchable_interp(
    const VP9_COMMON *cm, const MACROBLOCKD * xd) {
  const int pred_context = vp9_get_pred_context_switchable_interp(cm, xd);
  return &cm->fc.switchable_interp_prob[pred_context][0];
}

unsigned char vp9_get_pred_context_intra_inter(const VP9_COMMON *cm,
                                               const MACROBLOCKD *xd);
static INLINE vp9_prob vp9_get_pred_prob_intra_inter(const VP9_COMMON *cm,
                                                     const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_intra_inter(cm, xd);

  return cm->fc.intra_inter_prob[pred_context];
}

unsigned char vp9_get_pred_context_comp_inter_inter(const VP9_COMMON *cm,
                                                    const MACROBLOCKD *xd);


static INLINE vp9_prob vp9_get_pred_prob_comp_inter_inter(const VP9_COMMON *cm,
                                                          const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_comp_inter_inter(cm, xd);
  return cm->fc.comp_inter_prob[pred_context];
}

unsigned char vp9_get_pred_context_comp_ref_p(const VP9_COMMON *cm,
                                              const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_comp_ref_p(const VP9_COMMON *cm,
                                                    const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_comp_ref_p(cm, xd);
  return cm->fc.comp_ref_prob[pred_context];
}

unsigned char vp9_get_pred_context_single_ref_p1(const VP9_COMMON *cm,
                                                 const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_single_ref_p1(const VP9_COMMON *cm,
                                                       const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_single_ref_p1(cm, xd);
  return cm->fc.single_ref_prob[pred_context][0];
}

unsigned char vp9_get_pred_context_single_ref_p2(const VP9_COMMON *cm,
                                                 const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_single_ref_p2(const VP9_COMMON *cm,
                                                       const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_single_ref_p2(cm, xd);
  return cm->fc.single_ref_prob[pred_context][1];
}

unsigned char vp9_get_pred_context_tx_size(const VP9_COMMON *cm,
                                           const MACROBLOCKD *xd);

static INLINE const vp9_prob *vp9_get_pred_probs_tx_size(const VP9_COMMON *cm,
                                                         const MACROBLOCKD * xd) {
  const MODE_INFO *const mi = xd->mode_info_context;
  const int pred_context = vp9_get_pred_context_tx_size(cm, xd);
  if (mi->mbmi.sb_type < BLOCK_SIZE_MB16X16)
    return cm->fc.tx_probs_8x8p[pred_context];
  else if (mi->mbmi.sb_type < BLOCK_SIZE_SB32X32)
    return cm->fc.tx_probs_16x16p[pred_context];
  else
    return cm->fc.tx_probs_32x32p[pred_context];
}

#endif  // VP9_COMMON_VP9_PRED_COMMON_H_
