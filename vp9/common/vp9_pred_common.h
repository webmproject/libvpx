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

#ifdef __cplusplus
extern "C" {
#endif

static INLINE const MODE_INFO *get_above_mi(const MACROBLOCKD *const xd) {
  return xd->up_available ? xd->mi[-xd->mi_stride].src_mi : NULL;
}

static INLINE const MODE_INFO *get_left_mi(const MACROBLOCKD *const xd) {
  return xd->left_available ? xd->mi[-1].src_mi : NULL;
}

int vp9_get_segment_id(const VP9_COMMON *cm, const uint8_t *segment_ids,
                       BLOCK_SIZE bsize, int mi_row, int mi_col);

static INLINE int vp9_get_pred_context_seg_id(const MACROBLOCKD *xd) {
  const MODE_INFO *const above_mi = get_above_mi(xd);
  const MODE_INFO *const left_mi = get_left_mi(xd);
  const int above_sip = (above_mi != NULL) ?
                        above_mi->mbmi.seg_id_predicted : 0;
  const int left_sip = (left_mi != NULL) ? left_mi->mbmi.seg_id_predicted : 0;

  return above_sip + left_sip;
}

static INLINE vp9_prob vp9_get_pred_prob_seg_id(const struct segmentation *seg,
                                                const MACROBLOCKD *xd) {
  return seg->pred_probs[vp9_get_pred_context_seg_id(xd)];
}

static INLINE int vp9_get_skip_context(const MACROBLOCKD *xd) {
  const MODE_INFO *const above_mi = get_above_mi(xd);
  const MODE_INFO *const left_mi = get_left_mi(xd);
  const int above_skip = (above_mi != NULL) ? above_mi->mbmi.skip : 0;
  const int left_skip = (left_mi != NULL) ? left_mi->mbmi.skip : 0;
  return above_skip + left_skip;
}

static INLINE vp9_prob vp9_get_skip_prob(const VP9_COMMON *cm,
                                         const MACROBLOCKD *xd) {
  return cm->fc.skip_probs[vp9_get_skip_context(xd)];
}

#if CONFIG_SR_MODE
#include "vp9/common/vp9_sr_txfm.h"
static INLINE int vp9_get_sr_context(const MACROBLOCKD *xd,
                                     BLOCK_SIZE bsize) {
  TX_SIZE max_tx_size = max_txsize_lookup[bsize];
  int ctx;
  (void)xd;

  assert(max_tx_size >= MIN_SR_TX_SIZE &&
         max_tx_size <= MAX_SR_TX_SIZE);
  ctx = max_tx_size - MIN_SR_TX_SIZE;

  return ctx;
}

static INLINE vp9_prob vp9_get_sr_prob(const VP9_COMMON *cm,
                                       const MACROBLOCKD *xd,
                                       BLOCK_SIZE bsize) {
  int sr_ctx = vp9_get_sr_context(xd, bsize);
  assert(sr_ctx >= 0 && sr_ctx < SR_CONTEXTS);
  return cm->fc.sr_probs[sr_ctx];
}

#if SR_USE_MULTI_F
static INLINE vp9_prob vp9_get_sr_usfilter_context(const MACROBLOCKD *xd) {
  (void) xd;
  return 0;

  /*const MODE_INFO *const above_mi = get_above_mi(xd);
  const MODE_INFO *const left_mi = get_left_mi(xd);
  int above_sr_ver =
      (above_mi != NULL && above_mi->mbmi.sr && !above_mi->mbmi.skip) ?
      idx_to_v(above_mi->mbmi.us_filter_idx) : SR_USFILTER_NUM_D;
  int left_sr_hor =
      (left_mi != NULL && left_mi->mbmi.sr && !left_mi->mbmi.skip) ?
      idx_to_h(left_mi->mbmi.us_filter_idx) : SR_USFILTER_NUM_D;
  return above_sr_ver * 3 + left_sr_hor;*/
}

static INLINE const vp9_prob * vp9_get_sr_usfilter_prob(const VP9_COMMON *cm,
                                                const MACROBLOCKD *xd) {
  int sr_usfilter_ctx = vp9_get_sr_usfilter_context(xd);
  assert(sr_usfilter_ctx >= 0 && sr_usfilter_ctx < SR_USFILTER_CONTEXTS);
  return cm->fc.sr_usfilter_probs[sr_usfilter_ctx];
}
#endif  // SR_USE_MULTI_F
#endif  // CONFIG_SR_MODE

int vp9_get_pred_context_switchable_interp(const MACROBLOCKD *xd);

int vp9_get_intra_inter_context(const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_intra_inter_prob(const VP9_COMMON *cm,
                                                const MACROBLOCKD *xd) {
  return cm->fc.intra_inter_prob[vp9_get_intra_inter_context(xd)];
}

int vp9_get_reference_mode_context(const VP9_COMMON *cm, const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_reference_mode_prob(const VP9_COMMON *cm,
                                                   const MACROBLOCKD *xd) {
  return cm->fc.comp_inter_prob[vp9_get_reference_mode_context(cm, xd)];
}

int vp9_get_pred_context_comp_ref_p(const VP9_COMMON *cm,
                                    const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_comp_ref_p(const VP9_COMMON *cm,
                                                    const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_comp_ref_p(cm, xd);
  return cm->fc.comp_ref_probs[pred_context][0];
}

#if CONFIG_MULTI_REF
int vp9_get_pred_context_comp_ref_p1(const VP9_COMMON *cm,
                                     const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_comp_ref_p1(const VP9_COMMON *cm,
                                                     const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_comp_ref_p1(cm, xd);
  return cm->fc.comp_ref_probs[pred_context][1];
}

int vp9_get_pred_context_comp_ref_p2(const VP9_COMMON *cm,
                                     const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_comp_ref_p2(const VP9_COMMON *cm,
                                                     const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_comp_ref_p2(cm, xd);
  return cm->fc.comp_ref_probs[pred_context][2];
}

int vp9_get_pred_context_comp_ref_p3(const VP9_COMMON *cm,
                                     const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_comp_ref_p3(const VP9_COMMON *cm,
                                                     const MACROBLOCKD *xd) {
  const int pred_context = vp9_get_pred_context_comp_ref_p3(cm, xd);
  return cm->fc.comp_ref_probs[pred_context][3];
}
#endif  // CONFIG_MULTI_REF

int vp9_get_pred_context_single_ref_p1(const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_single_ref_p1(const VP9_COMMON *cm,
                                                       const MACROBLOCKD *xd) {
  return cm->fc.single_ref_probs[vp9_get_pred_context_single_ref_p1(xd)][0];
}

int vp9_get_pred_context_single_ref_p2(const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_single_ref_p2(const VP9_COMMON *cm,
                                                       const MACROBLOCKD *xd) {
  return cm->fc.single_ref_probs[vp9_get_pred_context_single_ref_p2(xd)][1];
}

#if CONFIG_MULTI_REF
int vp9_get_pred_context_single_ref_p3(const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_single_ref_p3(const VP9_COMMON *cm,
                                                       const MACROBLOCKD *xd) {
  return cm->fc.single_ref_probs[vp9_get_pred_context_single_ref_p3(xd)][2];
}

int vp9_get_pred_context_single_ref_p4(const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_single_ref_p4(const VP9_COMMON *cm,
                                                       const MACROBLOCKD *xd) {
  return cm->fc.single_ref_probs[vp9_get_pred_context_single_ref_p4(xd)][3];
}

int vp9_get_pred_context_single_ref_p5(const MACROBLOCKD *xd);

static INLINE vp9_prob vp9_get_pred_prob_single_ref_p5(const VP9_COMMON *cm,
                                                       const MACROBLOCKD *xd) {
  return cm->fc.single_ref_probs[vp9_get_pred_context_single_ref_p5(xd)][4];
}
#endif  // CONFIG_MULTI_REF

int vp9_get_tx_size_context(const MACROBLOCKD *xd);

static INLINE const vp9_prob *get_tx_probs(TX_SIZE max_tx_size, int ctx,
                                           const struct tx_probs *tx_probs) {
  switch (max_tx_size) {
    case TX_8X8:
      return tx_probs->p8x8[ctx];
    case TX_16X16:
      return tx_probs->p16x16[ctx];
    case TX_32X32:
      return tx_probs->p32x32[ctx];
#if CONFIG_TX64X64
    case TX_64X64:
      return tx_probs->p64x64[ctx];
#endif
    default:
      assert(0 && "Invalid max_tx_size.");
      return NULL;
  }
}

static INLINE const vp9_prob *get_tx_probs2(TX_SIZE max_tx_size,
                                            const MACROBLOCKD *xd,
                                            const struct tx_probs *tx_probs) {
  return get_tx_probs(max_tx_size, vp9_get_tx_size_context(xd), tx_probs);
}

static INLINE unsigned int *get_tx_counts(TX_SIZE max_tx_size, int ctx,
                                          struct tx_counts *tx_counts) {
  switch (max_tx_size) {
    case TX_8X8:
      return tx_counts->p8x8[ctx];
    case TX_16X16:
      return tx_counts->p16x16[ctx];
    case TX_32X32:
      return tx_counts->p32x32[ctx];
#if CONFIG_TX64X64
    case TX_64X64:
      return tx_counts->p64x64[ctx];
#endif
    default:
      assert(0 && "Invalid max_tx_size.");
      return NULL;
  }
}

#if CONFIG_SR_MODE
static INLINE unsigned int *get_real_tx_counts(TX_SIZE max_tx_size, int ctx,
                                          struct tx_counts *tx_counts) {
  switch (max_tx_size) {
    case TX_8X8:
      return tx_counts->real_p8x8[ctx];
    case TX_16X16:
      return tx_counts->real_p16x16[ctx];
    case TX_32X32:
      return tx_counts->real_p32x32[ctx];
#if CONFIG_TX64X64
    case TX_64X64:
      return tx_counts->real_p64x64[ctx];
#endif  // CONFIG_TX64X64
    default:
      assert(0 && "Invalid max_tx_size.");
      return NULL;
  }
}
#endif  // CONFIG_SR_MODE

#if CONFIG_COPY_MODE
int vp9_get_copy_mode_context(const MACROBLOCKD *xd);
#endif  // CONFIG_COPY_MODE

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_PRED_COMMON_H_
