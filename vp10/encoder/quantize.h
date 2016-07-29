/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_QUANTIZE_H_
#define VP10_ENCODER_QUANTIZE_H_

#include "./vpx_config.h"
#include "vp10/common/scan.h"
#include "vp10/encoder/block.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QUANT_PARAM {
  int log_scale;
} QUANT_PARAM;

typedef void (*VP10_QUANT_FACADE)(const tran_low_t *coeff_ptr,
                                  intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
                                  tran_low_t *qcoeff_ptr,
                                  const MACROBLOCKD_PLANE *pd,
                                  tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                                  const scan_order *sc,
                                  const QUANT_PARAM *qparam);

typedef struct {
#if CONFIG_NEW_QUANT
  DECLARE_ALIGNED(16, tran_low_t,
                  y_cuml_bins_nuq[QUANT_PROFILES][QINDEX_RANGE][COEF_BANDS]
                               [NUQ_KNOTS]);
  DECLARE_ALIGNED(16, tran_low_t,
                  uv_cuml_bins_nuq[QUANT_PROFILES][QINDEX_RANGE][COEF_BANDS]
                                [NUQ_KNOTS]);
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

struct VP10_COMP;
struct VP10Common;

void vp10_frame_init_quantizer(struct VP10_COMP *cpi);

void vp10_init_plane_quantizers(const struct VP10_COMP *cpi, MACROBLOCK *x,
                                int segment_id);

void vp10_init_quantizer(struct VP10_COMP *cpi);

void vp10_set_quantizer(struct VP10Common *cm, int q);

int vp10_quantizer_to_qindex(int quantizer);

int vp10_qindex_to_quantizer(int qindex);

void vp10_quantize_skip(intptr_t n_coeffs, tran_low_t *qcoeff_ptr,
                        tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr);

void vp10_quantize_fp_facade(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                             const MACROBLOCK_PLANE *p, tran_low_t *qcoeff_ptr,
                             const MACROBLOCKD_PLANE *pd,
                             tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                             const scan_order *sc, const QUANT_PARAM *qparam);

void vp10_quantize_b_facade(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                            const MACROBLOCK_PLANE *p, tran_low_t *qcoeff_ptr,
                            const MACROBLOCKD_PLANE *pd,
                            tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                            const scan_order *sc, const QUANT_PARAM *qparam);

void vp10_quantize_dc_facade(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                             const MACROBLOCK_PLANE *p, tran_low_t *qcoeff_ptr,
                             const MACROBLOCKD_PLANE *pd,
                             tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                             const scan_order *sc, const QUANT_PARAM *qparam);

#if CONFIG_NEW_QUANT
void quantize_dc_nuq(const tran_low_t *coeff_ptr,
                     intptr_t n_coeffs,
                     int skip_block,
                     const int16_t quant,
                     const int16_t quant_shift,
                     const int16_t dequant,
                     const tran_low_t *cuml_bins_ptr,
                     const tran_low_t *dequant_val,
                     tran_low_t *qcoeff_ptr,
                     tran_low_t *dqcoeff_ptr,
                     uint16_t *eob_ptr);
void quantize_dc_32x32_nuq(const tran_low_t *coeff_ptr,
                           intptr_t n_coeffs,
                           int skip_block,
                           const int16_t quant,
                           const int16_t quant_shift,
                           const int16_t dequant,
                           const tran_low_t *cuml_bins_ptr,
                           const tran_low_t *dequant_val,
                           tran_low_t *qcoeff_ptr,
                           tran_low_t *dqcoeff_ptr,
                           uint16_t *eob_ptr);
void quantize_dc_fp_nuq(const tran_low_t *coeff_ptr,
                        intptr_t n_coeffs,
                        int skip_block,
                        const int16_t quant,
                        const int16_t dequant,
                        const tran_low_t *cuml_bins_ptr,
                        const tran_low_t *dequant_val,
                        tran_low_t *qcoeff_ptr,
                        tran_low_t *dqcoeff_ptr,
                        uint16_t *eob_ptr);
void quantize_dc_32x32_fp_nuq(const tran_low_t *coeff_ptr,
                              intptr_t n_coeffs,
                              int skip_block,
                              const int16_t quant,
                              const int16_t dequant,
                              const tran_low_t *cuml_bins_ptr,
                              const tran_low_t *dequant_val,
                              tran_low_t *qcoeff_ptr,
                              tran_low_t *dqcoeff_ptr,
                              uint16_t *eob_ptr);
#endif  // CONFIG_NEW_QUANT

#if CONFIG_VPX_HIGHBITDEPTH
void vp10_highbd_quantize_fp_facade(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
    tran_low_t *qcoeff_ptr, const MACROBLOCKD_PLANE *pd,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr, const scan_order *sc,
    const QUANT_PARAM *qparam);

void vp10_highbd_quantize_b_facade(const tran_low_t *coeff_ptr,
                                   intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
                                   tran_low_t *qcoeff_ptr,
                                   const MACROBLOCKD_PLANE *pd,
                                   tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                                   const scan_order *sc,
                                   const QUANT_PARAM *qparam);

void vp10_highbd_quantize_dc_facade(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
    tran_low_t *qcoeff_ptr, const MACROBLOCKD_PLANE *pd,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr, const scan_order *sc,
    const QUANT_PARAM *qparam);

void vp10_highbd_quantize_dc(const tran_low_t *coeff_ptr,
                            int n_coeffs, int skip_block,
                            const int16_t *round_ptr, const int16_t quant,
                            tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                            const int16_t dequant_ptr, uint16_t *eob_ptr,
                            const int log_scale);
#if CONFIG_NEW_QUANT
void highbd_quantize_dc_nuq(const tran_low_t *coeff_ptr,
                            intptr_t n_coeffs,
                            int skip_block,
                            const int16_t quant,
                            const int16_t quant_shift,
                            const int16_t dequant,
                            const tran_low_t *cuml_bins_ptr,
                            const tran_low_t *dequant_val,
                            tran_low_t *qcoeff_ptr,
                            tran_low_t *dqcoeff_ptr,
                            uint16_t *eob_ptr);
void highbd_quantize_dc_32x32_nuq(const tran_low_t *coeff_ptr,
                                  intptr_t n_coeffs,
                                  int skip_block,
                                  const int16_t quant,
                                  const int16_t quant_shift,
                                  const int16_t dequant,
                                  const tran_low_t *cuml_bins_ptr,
                                  const tran_low_t *dequant_val,
                                  tran_low_t *qcoeff_ptr,
                                  tran_low_t *dqcoeff_ptr,
                                  uint16_t *eob_ptr);
void highbd_quantize_dc_fp_nuq(const tran_low_t *coeff_ptr,
                               intptr_t n_coeffs,
                               int skip_block,
                               const int16_t quant,
                               const int16_t dequant,
                               const tran_low_t *cuml_bins_ptr,
                               const tran_low_t *dequant_val,
                               tran_low_t *qcoeff_ptr,
                               tran_low_t *dqcoeff_ptr,
                               uint16_t *eob_ptr);
void highbd_quantize_dc_32x32_fp_nuq(const tran_low_t *coeff_ptr,
                                     intptr_t n_coeffs,
                                     int skip_block,
                                     const int16_t quant,
                                     const int16_t dequant,
                                     const tran_low_t *cuml_bins_ptr,
                                     const tran_low_t *dequant_val,
                                     tran_low_t *qcoeff_ptr,
                                     tran_low_t *dqcoeff_ptr,
                                     uint16_t *eob_ptr);

#endif  // CONFIG_NEW_QUANT
#endif  // CONFIG_VPX_HIGHBITDEPTH

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_QUANTIZE_H_
