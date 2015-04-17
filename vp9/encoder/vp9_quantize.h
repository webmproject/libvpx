/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_QUANTIZE_H_
#define VP9_ENCODER_VP9_QUANTIZE_H_

#include "./vpx_config.h"
#include "vp9/encoder/vp9_block.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
#if CONFIG_NEW_QUANT
  DECLARE_ALIGNED(16, tran_low_t,
                  y_cumbins_nuq[QINDEX_RANGE][COEF_BANDS][NUQ_KNOTES]);
  DECLARE_ALIGNED(16, tran_low_t,
                  uv_cumbins_nuq[QINDEX_RANGE][COEF_BANDS][NUQ_KNOTES]);
#endif  // CONFIG_NEW_QUANT

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

#if CONFIG_TX_SKIP
  DECLARE_ALIGNED(16, int16_t, y_quant_pxd[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, y_quant_shift_pxd[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, y_zbin_pxd[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, y_round_pxd[QINDEX_RANGE][8]);

  // TODO(jingning): in progress of re-working the quantization. will decide
  // if we want to deprecate the current use of y_quant.
  DECLARE_ALIGNED(16, int16_t, y_quant_pxd_fp[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, uv_quant_pxd_fp[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, y_round_pxd_fp[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, uv_round_pxd_fp[QINDEX_RANGE][8]);

  DECLARE_ALIGNED(16, int16_t, uv_quant_pxd[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, uv_quant_shift_pxd[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, uv_zbin_pxd[QINDEX_RANGE][8]);
  DECLARE_ALIGNED(16, int16_t, uv_round_pxd[QINDEX_RANGE][8]);
#if CONFIG_NEW_QUANT
  DECLARE_ALIGNED(16, tran_low_t,
                  y_cumbins_nuq_pxd[QINDEX_RANGE][COEF_BANDS][NUQ_KNOTES]);
  DECLARE_ALIGNED(16, tran_low_t,
                  uv_cumbins_nuq_pxd[QINDEX_RANGE][COEF_BANDS][NUQ_KNOTES]);
#endif  // CONFIG_NEW_QUANT
#endif  // CONFIG_TX_SKIP
} QUANTS;

void vp9_quantize_dc(const tran_low_t *coeff_ptr, int skip_block,
                     const int16_t *round_ptr, const int16_t quant_ptr,
                     tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                     const int16_t dequant_ptr, uint16_t *eob_ptr);
void vp9_quantize_dc_32x32(const tran_low_t *coeff_ptr, int skip_block,
                           const int16_t *round_ptr, const int16_t quant_ptr,
                           tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                           const int16_t dequant_ptr, uint16_t *eob_ptr);
#if CONFIG_NEW_QUANT
void vp9_quantize_dc_nuq(const tran_low_t *coeff_ptr,
                         int skip_block,
                         const int16_t quant,
                         const int16_t quant_shift,
                         const int16_t dequant,
                         const tran_low_t *cumbins_ptr,
                         const tran_low_t *dequant_val,
                         tran_low_t *qcoeff_ptr,
                         tran_low_t *dqcoeff_ptr,
                         uint16_t *eob_ptr);
void vp9_quantize_dc_32x32_nuq(const tran_low_t *coeff_ptr,
                               int skip_block,
                               const int16_t quant,
                               const int16_t quant_shift,
                               const int16_t dequant,
                               const tran_low_t *cumbins_ptr,
                               const tran_low_t *dequant_val,
                               tran_low_t *qcoeff_ptr,
                               tran_low_t *dqcoeff_ptr,
                               uint16_t *eob_ptr);
void vp9_quantize_dc_fp_nuq(const tran_low_t *coeff_ptr,
                            int skip_block,
                            const int16_t quant,
                            const int16_t dequant,
                            const tran_low_t *cumbins_ptr,
                            const tran_low_t *dequant_val,
                            tran_low_t *qcoeff_ptr,
                            tran_low_t *dqcoeff_ptr,
                            uint16_t *eob_ptr);
void vp9_quantize_dc_32x32_fp_nuq(const tran_low_t *coeff_ptr,
                                  int skip_block,
                                  const int16_t quant,
                                  const int16_t dequant,
                                  const tran_low_t *cumbins_ptr,
                                  const tran_low_t *dequant_val,
                                  tran_low_t *qcoeff_ptr,
                                  tran_low_t *dqcoeff_ptr,
                                  uint16_t *eob_ptr);
#endif  // CONFIG_NEW_QUANT

#if CONFIG_TX64X64
void vp9_quantize_dc_64x64(const tran_low_t *coeff_ptr, int skip_block,
                           const int16_t *round_ptr, const int16_t quant_ptr,
                           tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                           const int16_t dequant_ptr, uint16_t *eob_ptr);
#if CONFIG_NEW_QUANT
void vp9_quantize_dc_64x64_nuq(const tran_low_t *coeff_ptr,
                               int skip_block,
                               const int16_t quant,
                               const int16_t quant_shift,
                               const int16_t dequant,
                               const tran_low_t *cumbins_ptr,
                               const tran_low_t *dequant_val,
                               tran_low_t *qcoeff_ptr,
                               tran_low_t *dqcoeff_ptr,
                               uint16_t *eob_ptr);
void vp9_quantize_dc_64x64_fp_nuq(const tran_low_t *coeff_ptr,
                                  int skip_block,
                                  const int16_t quant,
                                  const int16_t dequant,
                                  const tran_low_t *cumbins_ptr,
                                  const tran_low_t *dequant_val,
                                  tran_low_t *qcoeff_ptr,
                                  tran_low_t *dqcoeff_ptr,
                                  uint16_t *eob_ptr);
#endif  // CONFIG_NEW_QUANT
#endif  // CONFIG_TX64X64

void vp9_regular_quantize_b_4x4(MACROBLOCK *x, int plane, int block,
                                const int16_t *scan, const int16_t *iscan);

#if CONFIG_VP9_HIGHBITDEPTH
void vp9_highbd_quantize_dc(const tran_low_t *coeff_ptr, int skip_block,
                            const int16_t *round_ptr, const int16_t quant_ptr,
                            tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                            const int16_t dequant_ptr, uint16_t *eob_ptr);
void vp9_highbd_quantize_dc_32x32(const tran_low_t *coeff_ptr,
                                  int skip_block,
                                  const int16_t *round_ptr,
                                  const int16_t quant_ptr,
                                  tran_low_t *qcoeff_ptr,
                                  tran_low_t *dqcoeff_ptr,
                                  const int16_t dequant_ptr,
                                  uint16_t *eob_ptr);
#if CONFIG_NEW_QUANT
void vp9_highbd_quantize_dc_nuq(const tran_low_t *coeff_ptr,
                                int skip_block,
                                const int16_t quant,
                                const int16_t quant_shift,
                                const int16_t dequant,
                                const tran_low_t *cumbins_ptr,
                                const tran_low_t *dequant_val,
                                tran_low_t *qcoeff_ptr,
                                tran_low_t *dqcoeff_ptr,
                                uint16_t *eob_ptr);
void vp9_highbd_quantize_dc_32x32_nuq(const tran_low_t *coeff_ptr,
                                      int skip_block,
                                      const int16_t quant,
                                      const int16_t quant_shift,
                                      const int16_t dequant,
                                      const tran_low_t *cumbins_ptr,
                                      const tran_low_t *dequant_val,
                                      tran_low_t *qcoeff_ptr,
                                      tran_low_t *dqcoeff_ptr,
                                      uint16_t *eob_ptr);
void vp9_highbd_quantize_dc_fp_nuq(const tran_low_t *coeff_ptr,
                                   int skip_block,
                                   const int16_t quant,
                                   const int16_t dequant,
                                   const tran_low_t *cumbins_ptr,
                                   const tran_low_t *dequant_val,
                                   tran_low_t *qcoeff_ptr,
                                   tran_low_t *dqcoeff_ptr,
                                   uint16_t *eob_ptr);
void vp9_highbd_quantize_dc_32x32_fp_nuq(const tran_low_t *coeff_ptr,
                                         int skip_block,
                                         const int16_t quant,
                                         const int16_t dequant,
                                         const tran_low_t *cumbins_ptr,
                                         const tran_low_t *dequant_val,
                                         tran_low_t *qcoeff_ptr,
                                         tran_low_t *dqcoeff_ptr,
                                         uint16_t *eob_ptr);
#endif  // CONFIG_NEW_QUANT
#if CONFIG_TX64X64
void vp9_highbd_quantize_dc_64x64(const tran_low_t *coeff_ptr,
                                  int skip_block,
                                  const int16_t *round_ptr,
                                  const int16_t quant_ptr,
                                  tran_low_t *qcoeff_ptr,
                                  tran_low_t *dqcoeff_ptr,
                                  const int16_t dequant_ptr,
                                  uint16_t *eob_ptr);
#if CONFIG_NEW_QUANT
void vp9_highbd_quantize_dc_64x64_nuq(const tran_low_t *coeff_ptr,
                                      int skip_block,
                                      const int16_t quant,
                                      const int16_t quant_shift,
                                      const int16_t dequant,
                                      const tran_low_t *cumbins_ptr,
                                      const tran_low_t *dequant_val,
                                      tran_low_t *qcoeff_ptr,
                                      tran_low_t *dqcoeff_ptr,
                                      uint16_t *eob_ptr);
void vp9_highbd_quantize_dc_64x64_fp_nuq(const tran_low_t *coeff_ptr,
                                         int skip_block,
                                         const int16_t quant,
                                         const int16_t dequant,
                                         const tran_low_t *cumbins_ptr,
                                         const tran_low_t *dequant_val,
                                         tran_low_t *qcoeff_ptr,
                                         tran_low_t *dqcoeff_ptr,
                                         uint16_t *eob_ptr);
#endif  // CONFIG_NEW_QUANT
#endif  // CONFIG_TX64X64
#endif  // CONFIG_VP9_HIGHBITDEPTH

struct VP9_COMP;
struct VP9Common;

void vp9_frame_init_quantizer(struct VP9_COMP *cpi);

void vp9_init_plane_quantizers(struct VP9_COMP *cpi, MACROBLOCK *x);

void vp9_init_quantizer(struct VP9_COMP *cpi);

void vp9_set_quantizer(struct VP9Common *cm, int q);

int vp9_quantizer_to_qindex(int quantizer);

int vp9_qindex_to_quantizer(int qindex);

#if CONFIG_TX_SKIP
void vp9_quantize_rect(const tran_low_t *coeff_ptr, int row, int col,
                       const int16_t *zbin_ptr, const int16_t *round_ptr,
                       const int16_t *quant_ptr, const int16_t *quant_shift_ptr,
                       tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                       const int16_t *dequant_ptr,
                       int logsizeby32, int stride, int has_dc, int hbd);
#if CONFIG_NEW_QUANT
void vp9_quantize_rect_nuq(const tran_low_t *coeff_ptr,
                           int row,
                           int col,
                           int stride,
                           const int16_t *quant_ptr,
                           const int16_t *quant_shift_ptr,
                           const int16_t *dequant_ptr,
                           const cumbins_type_nuq *cumbins_ptr,
                           const dequant_val_type_nuq *dequant_val,
                           tran_low_t *qcoeff_ptr,
                           tran_low_t *dqcoeff_ptr,
                           uint16_t *eob_ptr,
                           int logsizeby32,
                           const int16_t *scan,
                           const uint8_t *band);
#endif  // CONFIG_NEW_QUANT
int get_eob(tran_low_t *qcoeff_ptr, intptr_t n_coeffs, const int16_t *scan);
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_QUANTIZE_H_
