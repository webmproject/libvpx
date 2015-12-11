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

typedef void (*VP10_QUANT_FACADE)(const tran_low_t *coeff_ptr,
                                  intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
                                  tran_low_t *qcoeff_ptr,
                                  const MACROBLOCKD_PLANE *pd,
                                  tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                                  const scan_order *sc);

typedef struct {
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

void vp10_regular_quantize_b_4x4(MACROBLOCK *x, int plane, int block,
                                 const int16_t *scan, const int16_t *iscan);

struct VP10_COMP;
struct VP10Common;

void vp10_frame_init_quantizer(struct VP10_COMP *cpi);

void vp10_init_plane_quantizers(struct VP10_COMP *cpi, MACROBLOCK *x);

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
                             const scan_order *sc);

void vp10_quantize_b_facade(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                            const MACROBLOCK_PLANE *p, tran_low_t *qcoeff_ptr,
                            const MACROBLOCKD_PLANE *pd,
                            tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                            const scan_order *sc);

void vp10_quantize_dc_facade(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                             const MACROBLOCK_PLANE *p, tran_low_t *qcoeff_ptr,
                             const MACROBLOCKD_PLANE *pd,
                             tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                             const scan_order *sc);
#if CONFIG_VP9_HIGHBITDEPTH
void vp10_highbd_quantize_fp_facade(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
    tran_low_t *qcoeff_ptr, const MACROBLOCKD_PLANE *pd,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr, const scan_order *sc);

void vp10_highbd_quantize_b_facade(const tran_low_t *coeff_ptr,
                                   intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
                                   tran_low_t *qcoeff_ptr,
                                   const MACROBLOCKD_PLANE *pd,
                                   tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                                   const scan_order *sc);

void vp10_highbd_quantize_dc_facade(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
    tran_low_t *qcoeff_ptr, const MACROBLOCKD_PLANE *pd,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr, const scan_order *sc);
#endif  // CONFIG_VP9_HIGHBITDEPTH

void vp10_quantize_fp_32x32_facade(const tran_low_t *coeff_ptr,
                                   intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
                                   tran_low_t *qcoeff_ptr,
                                   const MACROBLOCKD_PLANE *pd,
                                   tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                                   const scan_order *sc);

void vp10_quantize_b_32x32_facade(const tran_low_t *coeff_ptr,
                                  intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
                                  tran_low_t *qcoeff_ptr,
                                  const MACROBLOCKD_PLANE *pd,
                                  tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                                  const scan_order *sc);

void vp10_quantize_dc_32x32_facade(const tran_low_t *coeff_ptr,
                                   intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
                                   tran_low_t *qcoeff_ptr,
                                   const MACROBLOCKD_PLANE *pd,
                                   tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                                   const scan_order *sc);
#if CONFIG_VP9_HIGHBITDEPTH
void vp10_highbd_quantize_fp_32x32_facade(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
    tran_low_t *qcoeff_ptr, const MACROBLOCKD_PLANE *pd,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr, const scan_order *sc);

void vp10_highbd_quantize_b_32x32_facade(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
    tran_low_t *qcoeff_ptr, const MACROBLOCKD_PLANE *pd,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr, const scan_order *sc);

void vp10_highbd_quantize_dc_32x32_facade(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
    tran_low_t *qcoeff_ptr, const MACROBLOCKD_PLANE *pd,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr, const scan_order *sc);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_QUANTIZE_H_
