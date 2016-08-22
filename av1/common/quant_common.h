/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_COMMON_QUANT_COMMON_H_
#define VP10_COMMON_QUANT_COMMON_H_

#include "aom/vpx_codec.h"
#include "av1/common/seg_common.h"
#include "av1/common/enums.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MINQ 0
#define MAXQ 255
#define QINDEX_RANGE (MAXQ - MINQ + 1)
#define QINDEX_BITS 8
#if CONFIG_AOM_QM
// Total number of QM sets stored
#define QM_LEVEL_BITS 4
#define NUM_QM_LEVELS (1 << QM_LEVEL_BITS)
/* Offset into the list of QMs. Actual number of levels used is
   (NUM_QM_LEVELS-AOM_QM_OFFSET)
   Lower value of AOM_QM_OFFSET implies more heavily weighted matrices.*/
#define DEFAULT_QM_FIRST (NUM_QM_LEVELS / 2)
#define DEFAULT_QM_LAST (NUM_QM_LEVELS - 1)
#endif

struct VP10Common;

int16_t vp10_dc_quant(int qindex, int delta, vpx_bit_depth_t bit_depth);
int16_t vp10_ac_quant(int qindex, int delta, vpx_bit_depth_t bit_depth);

int vp10_get_qindex(const struct segmentation *seg, int segment_id,
                    int base_qindex);
#if CONFIG_AOM_QM
// Reduce the large number of quantizers to a smaller number of levels for which
// different matrices may be defined
static inline int aom_get_qmlevel(int qindex, int first, int last) {
  int qmlevel = (qindex * (last + 1 - first) + QINDEX_RANGE / 2) / QINDEX_RANGE;
  qmlevel = VPXMIN(qmlevel + first, NUM_QM_LEVELS - 1);
  return qmlevel;
}
void aom_qm_init(struct VP10Common *cm);
qm_val_t *aom_iqmatrix(struct VP10Common *cm, int qindex, int comp,
                       int log2sizem2, int is_intra);
qm_val_t *aom_qmatrix(struct VP10Common *cm, int qindex, int comp,
                      int log2sizem2, int is_intra);
#endif

#if CONFIG_NEW_QUANT

#define QUANT_PROFILES 3
#define QUANT_RANGES 2
#define NUQ_KNOTS 3

typedef tran_low_t dequant_val_type_nuq[NUQ_KNOTS + 1];
typedef tran_low_t cuml_bins_type_nuq[NUQ_KNOTS];
void vp10_get_dequant_val_nuq(int q, int qindex, int band, tran_low_t *dq,
                              tran_low_t *cuml_bins, int dq_off_index);
tran_low_t vp10_dequant_abscoeff_nuq(int v, int q, const tran_low_t *dq);
tran_low_t vp10_dequant_coeff_nuq(int v, int q, const tran_low_t *dq);

static INLINE int get_dq_profile_from_ctx(int q_ctx) {
  return VPXMIN(q_ctx, QUANT_PROFILES - 1);
}
#endif  // CONFIG_NEW_QUANT

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_COMMON_QUANT_COMMON_H_
