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
                       BLOCK_SIZE bsize, int mi_row, int mi_col);


static INLINE int vp9_get_pred_context_seg_id(const MACROBLOCKD *xd) {
  const MODE_INFO * const above_mi = xd->mi_8x8[-xd->mode_info_stride];
  const MODE_INFO * const left_mi = xd->mi_8x8[-1];
  const int above_sip = above_mi ? above_mi->mbmi.seg_id_predicted : 0;
  const int left_sip = left_mi ? left_mi->mbmi.seg_id_predicted : 0;

  return above_sip + (xd->left_available ? left_sip : 0);
}

static INLINE vp9_prob vp9_get_pred_prob_seg_id(struct segmentation *seg,
                                                const MACROBLOCKD *xd) {
  return seg->pred_probs[vp9_get_pred_context_seg_id(xd)];
}

void vp9_set_pred_flag_seg_id(MACROBLOCKD *xd, uint8_t pred_flag);

static INLINE int vp9_get_pred_context_mbskip(const MACROBLOCKD *xd) {
  const MODE_INFO * const above_mi = xd->mi_8x8[-xd->mode_info_stride];
  const MODE_INFO * const left_mi = xd->mi_8x8[-1];
  const int above_skip_coeff = above_mi ? above_mi->mbmi.skip_coeff : 0;
  const int left_skip_coeff = left_mi ? left_mi->mbmi.skip_coeff : 0;

  return above_skip_coeff + (xd->left_available ? left_skip_coeff : 0);
}

static INLINE vp9_prob vp9_get_pred_prob_mbskip(const VP9_COMMON *cm,
                                                const MACROBLOCKD *xd) {
  return cm->fc.mbskip_probs[vp9_get_pred_context_mbskip(xd)];
}

static INLINE unsigned char vp9_get_pred_flag_mbskip(const MACROBLOCKD *xd) {
  return xd->this_mi->mbmi.skip_coeff;
}

void vp9_set_pred_flag_mbskip(MACROBLOCKD *xd, BLOCK_SIZE bsize,
                              uint8_t pred_flag);

unsigned char vp9_get_pred_context_switchable_interp(const MACROBLOCKD *xd);

unsigned char vp9_get_pred_context_intra_inter(const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_intra_inter(const VP9_COMMON *cm,
                                                     const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_intra_inter(xd);
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

unsigned char vp9_get_pred_context_single_ref_p1(const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_single_ref_p1(const VP9_COMMON *cm,
                                                       const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_single_ref_p1(xd);
  return cm->fc.single_ref_prob[pred_context][0];
}

unsigned char vp9_get_pred_context_single_ref_p2(const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_single_ref_p2(const VP9_COMMON *cm,
                                                       const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_single_ref_p2(xd);
  return cm->fc.single_ref_prob[pred_context][1];
}

unsigned char vp9_get_pred_context_tx_size(const MACROBLOCKD *xd);

static const vp9_prob *get_tx_probs(BLOCK_SIZE bsize, uint8_t context,
                                    const struct tx_probs *tx_probs) {
  if (bsize < BLOCK_16X16)
    return tx_probs->p8x8[context];
  else if (bsize < BLOCK_32X32)
    return tx_probs->p16x16[context];
  else
    return tx_probs->p32x32[context];
}

static const vp9_prob *get_tx_probs2(const MACROBLOCKD *xd,
                                     const struct tx_probs *tx_probs,
                                     const MODE_INFO *m) {
  const BLOCK_SIZE bsize = m->mbmi.sb_type;
  const int context = vp9_get_pred_context_tx_size(xd);
  return get_tx_probs(bsize, context, tx_probs);
}

static void update_tx_counts(BLOCK_SIZE bsize, uint8_t context,
                             TX_SIZE tx_size, struct tx_counts *tx_counts) {
  if (bsize >= BLOCK_32X32)
    tx_counts->p32x32[context][tx_size]++;
  else if (bsize >= BLOCK_16X16)
    tx_counts->p16x16[context][tx_size]++;
  else
    tx_counts->p8x8[context][tx_size]++;
}

#endif  // VP9_COMMON_VP9_PRED_COMMON_H_
