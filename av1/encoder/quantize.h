/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AV1_ENCODER_QUANTIZE_H_
#define AV1_ENCODER_QUANTIZE_H_

#include "./aom_config.h"
#include "av1/common/quant_common.h"
#include "av1/common/scan.h"
#include "av1/encoder/block.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QUANT_PARAM { int log_scale; } QUANT_PARAM;

typedef void (*AV1_QUANT_FACADE)(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                                 const MACROBLOCK_PLANE *p,
                                 tran_low_t *qcoeff_ptr,
                                 const MACROBLOCKD_PLANE *pd,
                                 tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                                 const SCAN_ORDER *sc, const QUANT_PARAM *qparam
#if CONFIG_AOM_QM
                                 ,
                                 const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr
#endif
                                 );

typedef struct {
#if CONFIG_NEW_QUANT
  DECLARE_ALIGNED(16, tran_low_t, y_cuml_bins_nuq[QUANT_PROFILES][QINDEX_RANGE]
                                                 [COEF_BANDS][NUQ_KNOTS]);
  DECLARE_ALIGNED(16, tran_low_t, uv_cuml_bins_nuq[QUANT_PROFILES][QINDEX_RANGE]
                                                  [COEF_BANDS][NUQ_KNOTS]);
#endif  // CONFIG_NEW_QUANT
  // 0: dc 1: ac 2-8: ac repeated to SIMD width
  DECLARE_ALIGNED(16, int16_t, y_quant[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, y_quant_shift[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, y_zbin[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, y_round[QINDEX_RANGE][8]);

  // TODO(jingning): in progress of re-working the quantization. will decide
  // if we want to deprecate the current use of y_quant.
  DECLARE_ALIGNED(16, int16_t, y_quant_fp[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, uv_quant_fp[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, y_round_fp[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, uv_round_fp[QINDEX_RANGE][8]);

  DECLARE_ALIGNED(16, int16_t, uv_quant[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, uv_quant_shift[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, uv_zbin[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, uv_round[QINDEX_RANGE][8]);
} QUANTS;

struct AV1_COMP;
struct AV1Common;

void av1_frame_init_quantizer(struct AV1_COMP *cpi);

void av1_init_plane_quantizers(const struct AV1_COMP *cpi, MACROBLOCK *x,
                               int segment_id);

void av1_init_quantizer(struct AV1_COMP *cpi);

void av1_set_quantizer(struct AV1Common *cm, int q);

int av1_quantizer_to_qindex(int quantizer);

int av1_qindex_to_quantizer(int qindex);

void av1_quantize_skip(intptr_t n_coeffs, tran_low_t *qcoeff_ptr,
                       tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr);

void av1_quantize_fp_facade(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                            const MACROBLOCK_PLANE *p, tran_low_t *qcoeff_ptr,
                            const MACROBLOCKD_PLANE *pd,
                            tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                            const SCAN_ORDER *sc, const QUANT_PARAM *qparam
#if CONFIG_AOM_QM
                            ,
                            const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr
#endif
                            );

void av1_quantize_b_facade(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                           const MACROBLOCK_PLANE *p, tran_low_t *qcoeff_ptr,
                           const MACROBLOCKD_PLANE *pd, tran_low_t *dqcoeff_ptr,
                           uint16_t *eob_ptr, const SCAN_ORDER *sc,
                           const QUANT_PARAM *qparam
#if CONFIG_AOM_QM
                           ,
                           const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr
#endif
                           );

void av1_quantize_dc_facade(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                            const MACROBLOCK_PLANE *p, tran_low_t *qcoeff_ptr,
                            const MACROBLOCKD_PLANE *pd,
                            tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                            const SCAN_ORDER *sc, const QUANT_PARAM *qparam
#if CONFIG_AOM_QM
                            ,
                            const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr
#endif
                            );

#if CONFIG_NEW_QUANT
void quantize_dc_nuq(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                     int skip_block, const int16_t quant,
                     const int16_t quant_shift, const int16_t dequant,
                     const tran_low_t *cuml_bins_ptr,
                     const tran_low_t *dequant_val, tran_low_t *qcoeff_ptr,
                     tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr);
void quantize_dc_32x32_nuq(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                           int skip_block, const int16_t quant,
                           const int16_t quant_shift, const int16_t dequant,
                           const tran_low_t *cuml_bins_ptr,
                           const tran_low_t *dequant_val,
                           tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                           uint16_t *eob_ptr);
void quantize_dc_fp_nuq(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                        int skip_block, const int16_t quant,
                        const int16_t dequant, const tran_low_t *cuml_bins_ptr,
                        const tran_low_t *dequant_val, tran_low_t *qcoeff_ptr,
                        tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr);
void quantize_dc_32x32_fp_nuq(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                              int skip_block, const int16_t quant,
                              const int16_t dequant,
                              const tran_low_t *cuml_bins_ptr,
                              const tran_low_t *dequant_val,
                              tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                              uint16_t *eob_ptr);
#endif  // CONFIG_NEW_QUANT

#if CONFIG_AOM_HIGHBITDEPTH
void av1_highbd_quantize_fp_facade(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
    tran_low_t *qcoeff_ptr, const MACROBLOCKD_PLANE *pd,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr, const SCAN_ORDER *sc,
    const QUANT_PARAM *qparam
#if CONFIG_AOM_QM
    ,
    const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr
#endif
    );

void av1_highbd_quantize_b_facade(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
    tran_low_t *qcoeff_ptr, const MACROBLOCKD_PLANE *pd,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr, const SCAN_ORDER *sc,
    const QUANT_PARAM *qparam
#if CONFIG_AOM_QM
    ,
    const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr
#endif
    );

void av1_highbd_quantize_dc_facade(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
    tran_low_t *qcoeff_ptr, const MACROBLOCKD_PLANE *pd,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr, const SCAN_ORDER *sc,
    const QUANT_PARAM *qparam
#if CONFIG_AOM_QM
    ,
    const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr
#endif
    );

void av1_highbd_quantize_dc(const tran_low_t *coeff_ptr, int n_coeffs,
                            int skip_block, const int16_t *round_ptr,
                            const int16_t quant, tran_low_t *qcoeff_ptr,
                            tran_low_t *dqcoeff_ptr, const int16_t dequant_ptr,
                            uint16_t *eob_ptr,
#if CONFIG_AOM_QM
                            const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr,
#endif
                            const int log_scale);
#if CONFIG_NEW_QUANT
void highbd_quantize_dc_nuq(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                            int skip_block, const int16_t quant,
                            const int16_t quant_shift, const int16_t dequant,
                            const tran_low_t *cuml_bins_ptr,
                            const tran_low_t *dequant_val,
                            tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                            uint16_t *eob_ptr);
void highbd_quantize_dc_32x32_nuq(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block,
    const int16_t quant, const int16_t quant_shift, const int16_t dequant,
    const tran_low_t *cuml_bins_ptr, const tran_low_t *dequant_val,
    tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr);
void highbd_quantize_dc_fp_nuq(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                               int skip_block, const int16_t quant,
                               const int16_t dequant,
                               const tran_low_t *cuml_bins_ptr,
                               const tran_low_t *dequant_val,
                               tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                               uint16_t *eob_ptr);
void highbd_quantize_dc_32x32_fp_nuq(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block,
    const int16_t quant, const int16_t dequant, const tran_low_t *cuml_bins_ptr,
    const tran_low_t *dequant_val, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr);

#endif  // CONFIG_NEW_QUANT
#endif  // CONFIG_AOM_HIGHBITDEPTH

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_QUANTIZE_H_
