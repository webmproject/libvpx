/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>
#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/quantize.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem.h"

#include "vp10/common/quant_common.h"
#include "vp10/common/scan.h"
#include "vp10/common/seg_common.h"

#include "vp10/encoder/encoder.h"
#include "vp10/encoder/quantize.h"
#include "vp10/encoder/rd.h"

#if CONFIG_NEW_QUANT
static INLINE int quantize_coeff_nuq(const tran_low_t coeffv,
                                     const int16_t quant,
                                     const int16_t quant_shift,
                                     const int16_t dequant,
                                     const tran_low_t *cuml_bins_ptr,
                                     const tran_low_t *dequant_val,
                                     tran_low_t *qcoeff_ptr,
                                     tran_low_t *dqcoeff_ptr) {
  const int coeff = coeffv;
  const int coeff_sign = (coeff >> 31);
  const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
  int i, q;
  int tmp = clamp(abs_coeff, INT16_MIN, INT16_MAX);
  for (i = 0; i < NUQ_KNOTS; i++) {
    if (tmp < cuml_bins_ptr[i]) {
      q = i;
      break;
    }
  }
  if (i == NUQ_KNOTS) {
    tmp -= cuml_bins_ptr[NUQ_KNOTS - 1];
    q = NUQ_KNOTS + (((((tmp * quant) >> 16) + tmp) * quant_shift) >> 16);
  }
  if (q) {
    *dqcoeff_ptr =
        vp10_dequant_abscoeff_nuq(q, dequant, dequant_val);
    *qcoeff_ptr  = (q ^ coeff_sign) - coeff_sign;
    *dqcoeff_ptr = *qcoeff_ptr < 0 ? -*dqcoeff_ptr : *dqcoeff_ptr;
  } else {
    *qcoeff_ptr = 0;
    *dqcoeff_ptr = 0;
  }
  return (q != 0);
}

static INLINE int quantize_coeff_bigtx_nuq(const tran_low_t coeffv,
                                           const int16_t quant,
                                           const int16_t quant_shift,
                                           const int16_t dequant,
                                           const tran_low_t *cuml_bins_ptr,
                                           const tran_low_t *dequant_val,
                                           tran_low_t *qcoeff_ptr,
                                           tran_low_t *dqcoeff_ptr,
                                           int logsizeby32) {
  const int coeff = coeffv;
  const int coeff_sign = (coeff >> 31);
  const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
  int i, q;
  int tmp = clamp(abs_coeff, INT16_MIN, INT16_MAX);
  for (i = 0; i < NUQ_KNOTS; i++) {
    if (tmp < ROUND_POWER_OF_TWO(cuml_bins_ptr[i], 1 + logsizeby32)) {
      q = i;
      break;
    }
  }
  if (i == NUQ_KNOTS) {
    tmp -= ROUND_POWER_OF_TWO(cuml_bins_ptr[NUQ_KNOTS - 1], 1 + logsizeby32);
    q = NUQ_KNOTS +
        (((((tmp * quant) >> 16) + tmp) * quant_shift) >> (15 - logsizeby32));
  }
  if (q) {
    *dqcoeff_ptr =
         ROUND_POWER_OF_TWO(vp10_dequant_abscoeff_nuq(q, dequant, dequant_val),
                            1 + logsizeby32);
    // *dqcoeff_ptr = vp10_dequant_abscoeff_nuq(q, dequant, dequant_val) >>
    // (1 + logsizeby32);
    *qcoeff_ptr  = (q ^ coeff_sign) - coeff_sign;
    *dqcoeff_ptr = *qcoeff_ptr < 0 ? -*dqcoeff_ptr : *dqcoeff_ptr;
  } else {
    *qcoeff_ptr = 0;
    *dqcoeff_ptr = 0;
  }
  return (q != 0);
}

static INLINE int quantize_coeff_fp_nuq(const tran_low_t coeffv,
                                        const int16_t quant,
                                        const int16_t dequant,
                                        const tran_low_t *cuml_bins_ptr,
                                        const tran_low_t *dequant_val,
                                        tran_low_t *qcoeff_ptr,
                                        tran_low_t *dqcoeff_ptr) {
  const int coeff = coeffv;
  const int coeff_sign = (coeff >> 31);
  const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
  int i, q;
  int tmp = clamp(abs_coeff, INT16_MIN, INT16_MAX);
  for (i = 0; i < NUQ_KNOTS; i++) {
    if (tmp < cuml_bins_ptr[i]) {
      q = i;
      break;
    }
  }
  if (i == NUQ_KNOTS) {
    q = NUQ_KNOTS +
        ((((int64_t)tmp - cuml_bins_ptr[NUQ_KNOTS - 1]) * quant) >> 16);
  }
  if (q) {
    *dqcoeff_ptr =
        vp10_dequant_abscoeff_nuq(q, dequant, dequant_val);
    *qcoeff_ptr  = (q ^ coeff_sign) - coeff_sign;
    *dqcoeff_ptr = *qcoeff_ptr < 0 ? -*dqcoeff_ptr : *dqcoeff_ptr;
  } else {
    *qcoeff_ptr = 0;
    *dqcoeff_ptr = 0;
  }
  return (q != 0);
}

static INLINE int quantize_coeff_bigtx_fp_nuq(const tran_low_t coeffv,
                                              const int16_t quant,
                                              const int16_t dequant,
                                              const tran_low_t *cuml_bins_ptr,
                                              const tran_low_t *dequant_val,
                                              tran_low_t *qcoeff_ptr,
                                              tran_low_t *dqcoeff_ptr,
                                              int logsizeby32) {
  const int coeff = coeffv;
  const int coeff_sign = (coeff >> 31);
  const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
  int i, q;
  int tmp = clamp(abs_coeff, INT16_MIN, INT16_MAX);
  for (i = 0; i < NUQ_KNOTS; i++) {
    if (tmp < ROUND_POWER_OF_TWO(cuml_bins_ptr[i], 1 + logsizeby32)) {
      q = i;
      break;
    }
  }
  if (i == NUQ_KNOTS) {
    q = NUQ_KNOTS +
        ((((int64_t)tmp - ROUND_POWER_OF_TWO(cuml_bins_ptr[NUQ_KNOTS - 1],
                                             1 + logsizeby32)) * quant) >>
         (15 - logsizeby32));
  }
  if (q) {
    *dqcoeff_ptr =
        ROUND_POWER_OF_TWO(vp10_dequant_abscoeff_nuq(q, dequant, dequant_val),
                           1 + logsizeby32);
    // *dqcoeff_ptr = vp10_dequant_abscoeff_nuq(q, dequant, dequant_val) >>
    // (1 + logsizeby32);
    *qcoeff_ptr  = (q ^ coeff_sign) - coeff_sign;
    *dqcoeff_ptr = *qcoeff_ptr < 0 ? -*dqcoeff_ptr : *dqcoeff_ptr;
  } else {
    *qcoeff_ptr = 0;
    *dqcoeff_ptr = 0;
  }
  return (q != 0);
}

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
                     uint16_t *eob_ptr) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    const int rc = 0;
    if (quantize_coeff_nuq(coeff_ptr[rc],
                           quant,
                           quant_shift,
                           dequant,
                           cuml_bins_ptr,
                           dequant_val,
                           qcoeff_ptr,
                           dqcoeff_ptr))
      eob = 0;
  }
  *eob_ptr = eob + 1;
}

void quantize_dc_fp_nuq(const tran_low_t *coeff_ptr,
                        intptr_t n_coeffs,
                        int skip_block,
                        const int16_t quant,
                        const int16_t dequant,
                        const tran_low_t *cuml_bins_ptr,
                        const tran_low_t *dequant_val,
                        tran_low_t *qcoeff_ptr,
                        tran_low_t *dqcoeff_ptr,
                        uint16_t *eob_ptr) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    const int rc = 0;
    if (quantize_coeff_fp_nuq(coeff_ptr[rc],
                              quant,
                              dequant,
                              cuml_bins_ptr,
                              dequant_val,
                              qcoeff_ptr,
                              dqcoeff_ptr))
      eob = 0;
  }
  *eob_ptr = eob + 1;
}

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
                           uint16_t *eob_ptr) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    const int rc = 0;
    if (quantize_coeff_bigtx_nuq(coeff_ptr[rc],
                                 quant,
                                 quant_shift,
                                 dequant,
                                 cuml_bins_ptr,
                                 dequant_val,
                                 qcoeff_ptr,
                                 dqcoeff_ptr,
                                 0))
      eob = 0;
  }
  *eob_ptr = eob + 1;
}

void quantize_dc_32x32_fp_nuq(const tran_low_t *coeff_ptr,
                              intptr_t n_coeffs,
                              int skip_block,
                              const int16_t quant,
                              const int16_t dequant,
                              const tran_low_t *cuml_bins_ptr,
                              const tran_low_t *dequant_val,
                              tran_low_t *qcoeff_ptr,
                              tran_low_t *dqcoeff_ptr,
                              uint16_t *eob_ptr) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    const int rc = 0;
    if (quantize_coeff_bigtx_fp_nuq(coeff_ptr[rc],
                                    quant,
                                    dequant,
                                    cuml_bins_ptr,
                                    dequant_val,
                                    qcoeff_ptr,
                                    dqcoeff_ptr,
                                    0))
      eob = 0;
  }
  *eob_ptr = eob + 1;
}

void quantize_nuq_c(const tran_low_t *coeff_ptr,
                    intptr_t n_coeffs,
                    int skip_block,
                    const int16_t *quant_ptr,
                    const int16_t *quant_shift_ptr,
                    const int16_t *dequant_ptr,
                    const cuml_bins_type_nuq *cuml_bins_ptr,
                    const dequant_val_type_nuq *dequant_val,
                    tran_low_t *qcoeff_ptr,
                    tran_low_t *dqcoeff_ptr,
                    uint16_t *eob_ptr,
                    const int16_t *scan,
                    const uint8_t *band) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    int i;
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      if (quantize_coeff_nuq(coeff_ptr[rc],
                             quant_ptr[rc != 0],
                             quant_shift_ptr[rc != 0],
                             dequant_ptr[rc != 0],
                             cuml_bins_ptr[band[i]],
                             dequant_val[band[i]],
                             &qcoeff_ptr[rc],
                             &dqcoeff_ptr[rc]))
        eob = i;
    }
  }
  *eob_ptr = eob + 1;
}

void quantize_fp_nuq_c(const tran_low_t *coeff_ptr,
                       intptr_t n_coeffs,
                       int skip_block,
                       const int16_t *quant_ptr,
                       const int16_t *dequant_ptr,
                       const cuml_bins_type_nuq *cuml_bins_ptr,
                       const dequant_val_type_nuq *dequant_val,
                       tran_low_t *qcoeff_ptr,
                       tran_low_t *dqcoeff_ptr,
                       uint16_t *eob_ptr,
                       const int16_t *scan,
                       const uint8_t *band) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    int i;
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      if (quantize_coeff_fp_nuq(coeff_ptr[rc],
                                quant_ptr[rc != 0],
                                dequant_ptr[rc != 0],
                                cuml_bins_ptr[band[i]],
                                dequant_val[band[i]],
                                &qcoeff_ptr[rc],
                                &dqcoeff_ptr[rc]))
        eob = i;
    }
  }
  *eob_ptr = eob + 1;
}

void quantize_32x32_nuq_c(const tran_low_t *coeff_ptr,
                          intptr_t n_coeffs,
                          int skip_block,
                          const int16_t *quant_ptr,
                          const int16_t *quant_shift_ptr,
                          const int16_t *dequant_ptr,
                          const cuml_bins_type_nuq *cuml_bins_ptr,
                          const dequant_val_type_nuq *dequant_val,
                          tran_low_t *qcoeff_ptr,
                          tran_low_t *dqcoeff_ptr,
                          uint16_t *eob_ptr,
                          const int16_t *scan,
                          const uint8_t *band) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    int i;
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      if (quantize_coeff_bigtx_nuq(coeff_ptr[rc],
                                   quant_ptr[rc != 0],
                                   quant_shift_ptr[rc != 0],
                                   dequant_ptr[rc != 0],
                                   cuml_bins_ptr[band[i]],
                                   dequant_val[band[i]],
                                   &qcoeff_ptr[rc],
                                   &dqcoeff_ptr[rc],
                                   0))
        eob = i;
    }
  }
  *eob_ptr = eob + 1;
}

void quantize_32x32_fp_nuq_c(const tran_low_t *coeff_ptr,
                             intptr_t n_coeffs,
                             int skip_block,
                             const int16_t *quant_ptr,
                             const int16_t *dequant_ptr,
                             const cuml_bins_type_nuq *cuml_bins_ptr,
                             const dequant_val_type_nuq *dequant_val,
                             tran_low_t *qcoeff_ptr,
                             tran_low_t *dqcoeff_ptr,
                             uint16_t *eob_ptr,
                             const int16_t *scan,
                             const uint8_t *band) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    int i;
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      if (quantize_coeff_bigtx_fp_nuq(coeff_ptr[rc],
                                      quant_ptr[rc != 0],
                                      dequant_ptr[rc != 0],
                                      cuml_bins_ptr[band[i]],
                                      dequant_val[band[i]],
                                      &qcoeff_ptr[rc],
                                      &dqcoeff_ptr[rc],
                                      0))
        eob = i;
    }
  }
  *eob_ptr = eob + 1;
}
#endif  // CONFIG_NEW_QUANT

void vp10_quantize_skip(intptr_t n_coeffs, tran_low_t *qcoeff_ptr,
                        tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr) {
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  *eob_ptr = 0;
}

void vp10_quantize_fp_facade(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                             const MACROBLOCK_PLANE *p, tran_low_t *qcoeff_ptr,
                             const MACROBLOCKD_PLANE *pd,
                             tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                             const scan_order *sc, const QUANT_PARAM *qparam) {
  // obsolete skip_block
  const int skip_block = 0;

  if (qparam->log_scale == 0) {
    vp10_quantize_fp(coeff_ptr, n_coeffs, skip_block, p->zbin, p->round_fp,
                     p->quant_fp, p->quant_shift, qcoeff_ptr, dqcoeff_ptr,
                     pd->dequant, eob_ptr, sc->scan, sc->iscan);
  } else {
    vp10_quantize_fp_32x32(coeff_ptr, n_coeffs, skip_block, p->zbin,
                           p->round_fp, p->quant_fp, p->quant_shift, qcoeff_ptr,
                           dqcoeff_ptr, pd->dequant, eob_ptr, sc->scan,
                           sc->iscan);
  }
}

void vp10_quantize_b_facade(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                            const MACROBLOCK_PLANE *p, tran_low_t *qcoeff_ptr,
                            const MACROBLOCKD_PLANE *pd,
                            tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                            const scan_order *sc, const QUANT_PARAM *qparam) {
  // obsolete skip_block
  const int skip_block = 0;

  if (qparam->log_scale == 0) {
    vpx_quantize_b(coeff_ptr, n_coeffs, skip_block, p->zbin, p->round, p->quant,
                   p->quant_shift, qcoeff_ptr, dqcoeff_ptr, pd->dequant,
                   eob_ptr, sc->scan, sc->iscan);
  } else {
    vpx_quantize_b_32x32(coeff_ptr, n_coeffs, skip_block, p->zbin, p->round,
                         p->quant, p->quant_shift, qcoeff_ptr, dqcoeff_ptr,
                         pd->dequant, eob_ptr, sc->scan, sc->iscan);
  }
}

void vp10_quantize_dc_facade(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                             const MACROBLOCK_PLANE *p, tran_low_t *qcoeff_ptr,
                             const MACROBLOCKD_PLANE *pd,
                             tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                             const scan_order *sc, const QUANT_PARAM *qparam) {
  // obsolete skip_block
  const int skip_block = 0;
  (void)sc;
  if (qparam->log_scale == 0) {
    vpx_quantize_dc(coeff_ptr, (int)n_coeffs, skip_block, p->round,
                    p->quant_fp[0], qcoeff_ptr, dqcoeff_ptr, pd->dequant[0],
                    eob_ptr);
  } else {
    vpx_quantize_dc_32x32(coeff_ptr, skip_block, p->round, p->quant_fp[0],
                          qcoeff_ptr, dqcoeff_ptr, pd->dequant[0], eob_ptr);
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
void vp10_highbd_quantize_fp_facade(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
    tran_low_t *qcoeff_ptr, const MACROBLOCKD_PLANE *pd,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr, const scan_order *sc,
    const QUANT_PARAM *qparam) {
  // obsolete skip_block
  const int skip_block = 0;

  vp10_highbd_quantize_fp(coeff_ptr, n_coeffs, skip_block, p->zbin,
                          p->round_fp, p->quant_fp, p->quant_shift,
                          qcoeff_ptr, dqcoeff_ptr, pd->dequant, eob_ptr,
                          sc->scan, sc->iscan, qparam->log_scale);
}

void vp10_highbd_quantize_b_facade(const tran_low_t *coeff_ptr,
                                   intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
                                   tran_low_t *qcoeff_ptr,
                                   const MACROBLOCKD_PLANE *pd,
                                   tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
                                   const scan_order *sc,
                                   const QUANT_PARAM *qparam) {
  // obsolete skip_block
  const int skip_block = 0;

  vp10_highbd_quantize_b(coeff_ptr, n_coeffs, skip_block, p->zbin, p->round,
                         p->quant, p->quant_shift, qcoeff_ptr, dqcoeff_ptr,
                         pd->dequant, eob_ptr, sc->scan, sc->iscan,
                         qparam->log_scale);
}

void vp10_highbd_quantize_dc_facade(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, const MACROBLOCK_PLANE *p,
    tran_low_t *qcoeff_ptr, const MACROBLOCKD_PLANE *pd,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr, const scan_order *sc,
    const QUANT_PARAM *qparam) {
  // obsolete skip_block
  const int skip_block = 0;

  (void)sc;

  vp10_highbd_quantize_dc(coeff_ptr, (int)n_coeffs, skip_block, p->round,
                         p->quant_fp[0], qcoeff_ptr, dqcoeff_ptr,
                         pd->dequant[0], eob_ptr, qparam->log_scale);
}

#if CONFIG_NEW_QUANT
static INLINE int highbd_quantize_coeff_nuq(const tran_low_t coeffv,
                                            const int16_t quant,
                                            const int16_t quant_shift,
                                            const int16_t dequant,
                                            const tran_low_t *cuml_bins_ptr,
                                            const tran_low_t *dequant_val,
                                            tran_low_t *qcoeff_ptr,
                                            tran_low_t *dqcoeff_ptr) {
  const int coeff = coeffv;
  const int coeff_sign = (coeff >> 31);
  const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
  int i, q;
  int64_t tmp = clamp(abs_coeff, INT32_MIN, INT32_MAX);
  for (i = 0; i < NUQ_KNOTS; i++) {
    if (tmp < cuml_bins_ptr[i]) {
      q = i;
      break;
    }
  }
  if (i == NUQ_KNOTS) {
    tmp -= cuml_bins_ptr[NUQ_KNOTS - 1];
    q = NUQ_KNOTS + (((((tmp * quant) >> 16) + tmp) * quant_shift) >> 16);
  }
  if (q) {
    *dqcoeff_ptr =
        vp10_dequant_abscoeff_nuq(q, dequant, dequant_val);
    *qcoeff_ptr  = (q ^ coeff_sign) - coeff_sign;
    *dqcoeff_ptr = *qcoeff_ptr < 0 ? -*dqcoeff_ptr : *dqcoeff_ptr;
  } else {
    *qcoeff_ptr = 0;
    *dqcoeff_ptr = 0;
  }
  return (q != 0);
}

static INLINE int highbd_quantize_coeff_fp_nuq(const tran_low_t coeffv,
                                               const int16_t quant,
                                               const int16_t dequant,
                                               const tran_low_t *cuml_bins_ptr,
                                               const tran_low_t *dequant_val,
                                               tran_low_t *qcoeff_ptr,
                                               tran_low_t *dqcoeff_ptr) {
  const int coeff = coeffv;
  const int coeff_sign = (coeff >> 31);
  const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
  int i, q;
  int64_t tmp = clamp(abs_coeff, INT32_MIN, INT32_MAX);
  for (i = 0; i < NUQ_KNOTS; i++) {
    if (tmp < cuml_bins_ptr[i]) {
      q = i;
      break;
    }
  }
  if (i == NUQ_KNOTS) {
    q = NUQ_KNOTS +
        (((tmp - cuml_bins_ptr[NUQ_KNOTS - 1]) * quant) >> 16);
  }
  if (q) {
    *dqcoeff_ptr =
        vp10_dequant_abscoeff_nuq(q, dequant, dequant_val);
    *qcoeff_ptr  = (q ^ coeff_sign) - coeff_sign;
    *dqcoeff_ptr = *qcoeff_ptr < 0 ? -*dqcoeff_ptr : *dqcoeff_ptr;
  } else {
    *qcoeff_ptr = 0;
    *dqcoeff_ptr = 0;
  }
  return (q != 0);
}

static INLINE int highbd_quantize_coeff_bigtx_fp_nuq(
    const tran_low_t coeffv,
    const int16_t quant,
    const int16_t dequant,
    const tran_low_t *cuml_bins_ptr,
    const tran_low_t *dequant_val,
    tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr,
    int logsizeby32) {
  const int coeff = coeffv;
  const int coeff_sign = (coeff >> 31);
  const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
  int i, q;
  int64_t tmp = clamp(abs_coeff, INT32_MIN, INT32_MAX);
  for (i = 0; i < NUQ_KNOTS; i++) {
    if (tmp < ROUND_POWER_OF_TWO(cuml_bins_ptr[i], 1 + logsizeby32)) {
      q = i;
      break;
    }
  }
  if (i == NUQ_KNOTS) {
    q = NUQ_KNOTS +
        (((tmp - ROUND_POWER_OF_TWO(cuml_bins_ptr[NUQ_KNOTS - 1],
                                    1 + logsizeby32)) * quant) >>
         (15 - logsizeby32));
  }
  if (q) {
    *dqcoeff_ptr =
        ROUND_POWER_OF_TWO(vp10_dequant_abscoeff_nuq(q, dequant, dequant_val),
                           1 + logsizeby32);
    *qcoeff_ptr  = (q ^ coeff_sign) - coeff_sign;
    *dqcoeff_ptr = *qcoeff_ptr < 0 ? -*dqcoeff_ptr : *dqcoeff_ptr;
  } else {
    *qcoeff_ptr = 0;
    *dqcoeff_ptr = 0;
  }
  return (q != 0);
}

static INLINE int highbd_quantize_coeff_bigtx_nuq(const tran_low_t coeffv,
                                                  const int16_t quant,
                                                  const int16_t quant_shift,
                                                  const int16_t dequant,
                                                  const tran_low_t
                                                        *cuml_bins_ptr,
                                                  const tran_low_t *dequant_val,
                                                  tran_low_t *qcoeff_ptr,
                                                  tran_low_t *dqcoeff_ptr,
                                                  int logsizeby32) {
  const int coeff = coeffv;
  const int coeff_sign = (coeff >> 31);
  const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
  int i, q;
  int64_t tmp = clamp(abs_coeff, INT32_MIN, INT32_MAX);
  for (i = 0; i < NUQ_KNOTS; i++) {
    if (tmp < ROUND_POWER_OF_TWO(cuml_bins_ptr[i], 1 + logsizeby32)) {
      q = i;
      break;
    }
  }
  if (i == NUQ_KNOTS) {
    tmp -= ROUND_POWER_OF_TWO(cuml_bins_ptr[NUQ_KNOTS - 1], 1 + logsizeby32);
    q = NUQ_KNOTS +
        (((((tmp * quant) >> 16) + tmp) * quant_shift) >> (15 - logsizeby32));
  }
  if (q) {
    *dqcoeff_ptr =
        ROUND_POWER_OF_TWO(vp10_dequant_abscoeff_nuq(q, dequant, dequant_val),
                           1 + logsizeby32);
    *qcoeff_ptr  = (q ^ coeff_sign) - coeff_sign;
    *dqcoeff_ptr = *qcoeff_ptr < 0 ? -*dqcoeff_ptr : *dqcoeff_ptr;
  } else {
    *qcoeff_ptr = 0;
    *dqcoeff_ptr = 0;
  }
  return (q != 0);
}

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
                            uint16_t *eob_ptr) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    const int rc = 0;
    if (highbd_quantize_coeff_nuq(coeff_ptr[rc],
                                  quant,
                                  quant_shift,
                                  dequant,
                                  cuml_bins_ptr,
                                  dequant_val,
                                  qcoeff_ptr,
                                  dqcoeff_ptr))
      eob = 0;
  }
  *eob_ptr = eob + 1;
}

void highbd_quantize_dc_fp_nuq(const tran_low_t *coeff_ptr,
                               intptr_t n_coeffs,
                               int skip_block,
                               const int16_t quant,
                               const int16_t dequant,
                               const tran_low_t *cuml_bins_ptr,
                               const tran_low_t *dequant_val,
                               tran_low_t *qcoeff_ptr,
                               tran_low_t *dqcoeff_ptr,
                               uint16_t *eob_ptr) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    const int rc = 0;
    if (highbd_quantize_coeff_fp_nuq(coeff_ptr[rc],
                                     quant,
                                     dequant,
                                     cuml_bins_ptr,
                                     dequant_val,
                                     qcoeff_ptr,
                                     dqcoeff_ptr))
      eob = 0;
  }
  *eob_ptr = eob + 1;
}

void highbd_quantize_nuq_c(const tran_low_t *coeff_ptr,
                           intptr_t n_coeffs,
                           int skip_block,
                           const int16_t *quant_ptr,
                           const int16_t *quant_shift_ptr,
                           const int16_t *dequant_ptr,
                           const cuml_bins_type_nuq *cuml_bins_ptr,
                           const dequant_val_type_nuq *dequant_val,
                           tran_low_t *qcoeff_ptr,
                           tran_low_t *dqcoeff_ptr,
                           uint16_t *eob_ptr,
                           const int16_t *scan,
                           const uint8_t *band) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    int i;
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      if (highbd_quantize_coeff_nuq(coeff_ptr[rc],
                                    quant_ptr[rc != 0],
                                    quant_shift_ptr[rc != 0],
                                    dequant_ptr[rc != 0],
                                    cuml_bins_ptr[band[i]],
                                    dequant_val[band[i]],
                                    &qcoeff_ptr[rc],
                                    &dqcoeff_ptr[rc]))
        eob = i;
    }
  }
  *eob_ptr = eob + 1;
}

void highbd_quantize_32x32_nuq_c(const tran_low_t *coeff_ptr,
                                 intptr_t n_coeffs,
                                 int skip_block,
                                 const int16_t *quant_ptr,
                                 const int16_t *quant_shift_ptr,
                                 const int16_t *dequant_ptr,
                                 const cuml_bins_type_nuq *cuml_bins_ptr,
                                 const dequant_val_type_nuq *dequant_val,
                                 tran_low_t *qcoeff_ptr,
                                 tran_low_t *dqcoeff_ptr,
                                 uint16_t *eob_ptr,
                                 const int16_t *scan,
                                 const uint8_t *band) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    int i;
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      if (highbd_quantize_coeff_bigtx_nuq(coeff_ptr[rc],
                                          quant_ptr[rc != 0],
                                          quant_shift_ptr[rc != 0],
                                          dequant_ptr[rc != 0],
                                          cuml_bins_ptr[band[i]],
                                          dequant_val[band[i]],
                                          &qcoeff_ptr[rc],
                                          &dqcoeff_ptr[rc],
                                          0))
        eob = i;
    }
  }
  *eob_ptr = eob + 1;
}

void highbd_quantize_32x32_fp_nuq_c(const tran_low_t *coeff_ptr,
                                    intptr_t n_coeffs,
                                    int skip_block,
                                    const int16_t *quant_ptr,
                                    const int16_t *dequant_ptr,
                                    const cuml_bins_type_nuq *cuml_bins_ptr,
                                    const dequant_val_type_nuq *dequant_val,
                                    tran_low_t *qcoeff_ptr,
                                    tran_low_t *dqcoeff_ptr,
                                    uint16_t *eob_ptr,
                                    const int16_t *scan,
                                    const uint8_t *band) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    int i;
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      if (highbd_quantize_coeff_bigtx_fp_nuq(coeff_ptr[rc],
                                             quant_ptr[rc != 0],
                                             dequant_ptr[rc != 0],
                                             cuml_bins_ptr[band[i]],
                                             dequant_val[band[i]],
                                             &qcoeff_ptr[rc],
                                             &dqcoeff_ptr[rc],
                                             0))
        eob = i;
    }
  }
  *eob_ptr = eob + 1;
}

void highbd_quantize_fp_nuq_c(const tran_low_t *coeff_ptr,
                              intptr_t n_coeffs,
                              int skip_block,
                              const int16_t *quant_ptr,
                              const int16_t *dequant_ptr,
                              const cuml_bins_type_nuq *cuml_bins_ptr,
                              const dequant_val_type_nuq *dequant_val,
                              tran_low_t *qcoeff_ptr,
                              tran_low_t *dqcoeff_ptr,
                              uint16_t *eob_ptr,
                              const int16_t *scan,
                              const uint8_t *band) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    int i;
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      if (highbd_quantize_coeff_fp_nuq(coeff_ptr[rc],
                                       quant_ptr[rc != 0],
                                       dequant_ptr[rc != 0],
                                       cuml_bins_ptr[band[i]],
                                       dequant_val[band[i]],
                                       &qcoeff_ptr[rc],
                                       &dqcoeff_ptr[rc]))
        eob = i;
    }
  }
  *eob_ptr = eob + 1;
}

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
                                  uint16_t *eob_ptr) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    const int rc = 0;
    if (highbd_quantize_coeff_bigtx_nuq(coeff_ptr[rc],
                                        quant,
                                        quant_shift,
                                        dequant,
                                        cuml_bins_ptr,
                                        dequant_val,
                                        qcoeff_ptr,
                                        dqcoeff_ptr,
                                        0))
      eob = 0;
  }
  *eob_ptr = eob + 1;
}

void highbd_quantize_dc_32x32_fp_nuq(const tran_low_t *coeff_ptr,
                                     intptr_t n_coeffs,
                                     int skip_block,
                                     const int16_t quant,
                                     const int16_t dequant,
                                     const tran_low_t *cuml_bins_ptr,
                                     const tran_low_t *dequant_val,
                                     tran_low_t *qcoeff_ptr,
                                     tran_low_t *dqcoeff_ptr,
                                     uint16_t *eob_ptr) {
  int eob = -1;
  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));
  if (!skip_block) {
    const int rc = 0;
    if (highbd_quantize_coeff_bigtx_fp_nuq(coeff_ptr[rc],
                                           quant,
                                           dequant,
                                           cuml_bins_ptr,
                                           dequant_val,
                                           qcoeff_ptr,
                                           dqcoeff_ptr,
                                           0))
      eob = 0;
  }
  *eob_ptr = eob + 1;
}
#endif  // CONFIG_NEW_QUANT
#endif  // CONFIG_VP9_HIGHBITDEPTH

void vp10_quantize_fp_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                       int skip_block,
                       const int16_t *zbin_ptr, const int16_t *round_ptr,
                       const int16_t *quant_ptr, const int16_t *quant_shift_ptr,
                       tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                       const int16_t *dequant_ptr,
                       uint16_t *eob_ptr,
                       const int16_t *scan, const int16_t *iscan) {
  int i, eob = -1;
  // TODO(jingning) Decide the need of these arguments after the
  // quantization process is completed.
  (void)zbin_ptr;
  (void)quant_shift_ptr;
  (void)iscan;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    // Quantization pass: All coefficients with index >= zero_flag are
    // skippable. Note: zero_flag can be zero.
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      const int coeff = coeff_ptr[rc];
      const int coeff_sign = (coeff >> 31);
      const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;

      int tmp = clamp(abs_coeff + round_ptr[rc != 0], INT16_MIN, INT16_MAX);
      tmp = (tmp * quant_ptr[rc != 0]) >> 16;

      qcoeff_ptr[rc] = (tmp ^ coeff_sign) - coeff_sign;
      dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant_ptr[rc != 0];

      if (tmp)
        eob = i;
    }
  }
  *eob_ptr = eob + 1;
}

#if CONFIG_VP9_HIGHBITDEPTH
void vp10_highbd_quantize_fp_c(const tran_low_t *coeff_ptr,
                              intptr_t count,
                              int skip_block,
                              const int16_t *zbin_ptr,
                              const int16_t *round_ptr,
                              const int16_t *quant_ptr,
                              const int16_t *quant_shift_ptr,
                              tran_low_t *qcoeff_ptr,
                              tran_low_t *dqcoeff_ptr,
                              const int16_t *dequant_ptr,
                              uint16_t *eob_ptr,
                              const int16_t *scan,
                              const int16_t *iscan, const int log_scale) {
  int i;
  int eob = -1;
  const int scale = 1 << log_scale;
  const int shift = 16 - log_scale;
  // TODO(jingning) Decide the need of these arguments after the
  // quantization process is completed.
  (void)zbin_ptr;
  (void)quant_shift_ptr;
  (void)iscan;

  memset(qcoeff_ptr, 0, count * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, count * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    // Quantization pass: All coefficients with index >= zero_flag are
    // skippable. Note: zero_flag can be zero.
    for (i = 0; i < count; i++) {
      const int rc = scan[i];
      const int coeff = coeff_ptr[rc];
      const int coeff_sign = (coeff >> 31);
      const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
      const int64_t tmp = abs_coeff + round_ptr[rc != 0];
      const uint32_t abs_qcoeff =
          (uint32_t)((tmp * quant_ptr[rc != 0]) >> shift);
      qcoeff_ptr[rc] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
      dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant_ptr[rc != 0]  / scale;
      if (abs_qcoeff)
        eob = i;
    }
  }
  *eob_ptr = eob + 1;
}

#endif  // CONFIG_VP9_HIGHBITDEPTH

// TODO(jingning) Refactor this file and combine functions with similar
// operations.
void vp10_quantize_fp_32x32_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                             int skip_block,
                             const int16_t *zbin_ptr, const int16_t *round_ptr,
                             const int16_t *quant_ptr,
                             const int16_t *quant_shift_ptr,
                             tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                             const int16_t *dequant_ptr,
                             uint16_t *eob_ptr,
                             const int16_t *scan, const int16_t *iscan) {
  int i, eob = -1;
  (void)zbin_ptr;
  (void)quant_shift_ptr;
  (void)iscan;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      const int coeff = coeff_ptr[rc];
      const int coeff_sign = (coeff >> 31);
      int tmp = 0;
      int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;

      if (abs_coeff >= (dequant_ptr[rc != 0] >> 2)) {
        abs_coeff += ROUND_POWER_OF_TWO(round_ptr[rc != 0], 1);
        abs_coeff = clamp(abs_coeff, INT16_MIN, INT16_MAX);
        tmp = (abs_coeff * quant_ptr[rc != 0]) >> 15;
        qcoeff_ptr[rc] = (tmp ^ coeff_sign) - coeff_sign;
        dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant_ptr[rc != 0] / 2;
      }

      if (tmp)
        eob = i;
    }
  }
  *eob_ptr = eob + 1;
}

#if CONFIG_VP9_HIGHBITDEPTH
void vp10_highbd_quantize_b_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                              int skip_block, const int16_t *zbin_ptr,
                              const int16_t *round_ptr,
                              const int16_t *quant_ptr,
                              const int16_t *quant_shift_ptr,
                              tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                              const int16_t *dequant_ptr,
                              uint16_t *eob_ptr, const int16_t *scan,
                              const int16_t *iscan, const int log_scale) {
  int i, non_zero_count = (int)n_coeffs, eob = -1;
  int zbins[2] = {zbin_ptr[0], zbin_ptr[1]};
  int round[2] = {round_ptr[0], round_ptr[1]};
  int nzbins[2];
  int scale = 1;
  int shift = 16;
  (void)iscan;

  if (log_scale > 0) {
    zbins[0] = ROUND_POWER_OF_TWO(zbin_ptr[0], log_scale);
    zbins[1] = ROUND_POWER_OF_TWO(zbin_ptr[1], log_scale);
    round[0] = ROUND_POWER_OF_TWO(round_ptr[0], log_scale);
    round[1] = ROUND_POWER_OF_TWO(round_ptr[1], log_scale);
    scale = 1 << log_scale;
    shift = 16 - log_scale;
  }

  nzbins[0] = zbins[0] * -1;
  nzbins[1] = zbins[1] * -1;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    // Pre-scan pass
    for (i = (int)n_coeffs - 1; i >= 0; i--) {
      const int rc = scan[i];
      const int coeff = coeff_ptr[rc];

      if (coeff < zbins[rc != 0] && coeff > nzbins[rc != 0])
        non_zero_count--;
      else
        break;
    }

    // Quantization pass: All coefficients with index >= zero_flag are
    // skippable. Note: zero_flag can be zero.
    for (i = 0; i < non_zero_count; i++) {
      const int rc = scan[i];
      const int coeff = coeff_ptr[rc];
      const int coeff_sign = (coeff >> 31);
      const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;

      if (abs_coeff >= zbins[rc != 0]) {
        const int64_t tmp1 = abs_coeff + round[rc != 0];
        const int64_t tmp2 = ((tmp1 * quant_ptr[rc != 0]) >> 16) + tmp1;
        const uint32_t abs_qcoeff =
            (uint32_t)((tmp2 * quant_shift_ptr[rc != 0]) >> shift);
        qcoeff_ptr[rc] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
        dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant_ptr[rc != 0] / scale;
        if (abs_qcoeff)
          eob = i;
      }
    }
  }
  *eob_ptr = eob + 1;
}
#endif

#if CONFIG_VP9_HIGHBITDEPTH
void vp10_highbd_quantize_dc(const tran_low_t *coeff_ptr,
                            int n_coeffs, int skip_block,
                            const int16_t *round_ptr, const int16_t quant,
                            tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                            const int16_t dequant_ptr, uint16_t *eob_ptr,
                            const int log_scale) {
  int eob = -1;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    const int coeff = coeff_ptr[0];
    const int coeff_sign = (coeff >> 31);
    const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
    const int64_t tmp = abs_coeff + round_ptr[0];
    const uint32_t abs_qcoeff = (uint32_t)((tmp * quant) >> (16 - log_scale));
    qcoeff_ptr[0] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
    dqcoeff_ptr[0] = qcoeff_ptr[0] * dequant_ptr / (1 << log_scale);
    if (abs_qcoeff)
      eob = 0;
  }
  *eob_ptr = eob + 1;
}
#endif

static void invert_quant(int16_t *quant, int16_t *shift, int d) {
  unsigned t;
  int l, m;
  t = d;
  for (l = 0; t > 1; l++)
    t >>= 1;
  m = 1 + (1 << (16 + l)) / d;
  *quant = (int16_t)(m - (1 << 16));
  *shift = 1 << (16 - l);
}

static int get_qzbin_factor(int q, vpx_bit_depth_t bit_depth) {
  const int quant = vp10_dc_quant(q, 0, bit_depth);
#if CONFIG_VP9_HIGHBITDEPTH
  switch (bit_depth) {
    case VPX_BITS_8:
      return q == 0 ? 64 : (quant < 148 ? 84 : 80);
    case VPX_BITS_10:
      return q == 0 ? 64 : (quant < 592 ? 84 : 80);
    case VPX_BITS_12:
      return q == 0 ? 64 : (quant < 2368 ? 84 : 80);
    default:
      assert(0 && "bit_depth should be VPX_BITS_8, VPX_BITS_10 or VPX_BITS_12");
      return -1;
  }
#else
  (void) bit_depth;
  return q == 0 ? 64 : (quant < 148 ? 84 : 80);
#endif
}

void vp10_init_quantizer(VP10_COMP *cpi) {
  VP10_COMMON *const cm = &cpi->common;
  QUANTS *const quants = &cpi->quants;
  int i, q, quant;
#if CONFIG_NEW_QUANT
  int dq;
#endif

  for (q = 0; q < QINDEX_RANGE; q++) {
    const int qzbin_factor = get_qzbin_factor(q, cm->bit_depth);
    const int qrounding_factor = q == 0 ? 64 : 48;

    for (i = 0; i < 2; ++i) {
      int qrounding_factor_fp = 64;
      // y
      quant = i == 0 ? vp10_dc_quant(q, cm->y_dc_delta_q, cm->bit_depth)
                     : vp10_ac_quant(q, 0, cm->bit_depth);
      invert_quant(&quants->y_quant[q][i], &quants->y_quant_shift[q][i], quant);
      quants->y_quant_fp[q][i] = (1 << 16) / quant;
      quants->y_round_fp[q][i] = (qrounding_factor_fp * quant) >> 7;
      quants->y_zbin[q][i] = ROUND_POWER_OF_TWO(qzbin_factor * quant, 7);
      quants->y_round[q][i] = (qrounding_factor * quant) >> 7;
      cpi->y_dequant[q][i] = quant;

      // uv
      quant = i == 0 ? vp10_dc_quant(q, cm->uv_dc_delta_q, cm->bit_depth)
                     : vp10_ac_quant(q, cm->uv_ac_delta_q, cm->bit_depth);
      invert_quant(&quants->uv_quant[q][i],
                   &quants->uv_quant_shift[q][i], quant);
      quants->uv_quant_fp[q][i] = (1 << 16) / quant;
      quants->uv_round_fp[q][i] = (qrounding_factor_fp * quant) >> 7;
      quants->uv_zbin[q][i] = ROUND_POWER_OF_TWO(qzbin_factor * quant, 7);
      quants->uv_round[q][i] = (qrounding_factor * quant) >> 7;
      cpi->uv_dequant[q][i] = quant;
    }

#if CONFIG_NEW_QUANT
    for (dq = 0; dq < QUANT_PROFILES; dq++) {
      for (i = 0; i < COEF_BANDS; i++) {
        const int quant = cpi->y_dequant[q][i != 0];
        const int uvquant = cpi->uv_dequant[q][i != 0];
        vp10_get_dequant_val_nuq(quant, q, i,
                                 cpi->y_dequant_val_nuq[dq][q][i],
                                 quants->y_cuml_bins_nuq[dq][q][i], dq);
        vp10_get_dequant_val_nuq(uvquant, q, i,
                                 cpi->uv_dequant_val_nuq[dq][q][i],
                                 quants->uv_cuml_bins_nuq[dq][q][i], dq);
      }
    }
#endif  // CONFIG_NEW_QUANT

    for (i = 2; i < 8; i++) {  // 8: SIMD width
      quants->y_quant[q][i] = quants->y_quant[q][1];
      quants->y_quant_fp[q][i] = quants->y_quant_fp[q][1];
      quants->y_round_fp[q][i] = quants->y_round_fp[q][1];
      quants->y_quant_shift[q][i] = quants->y_quant_shift[q][1];
      quants->y_zbin[q][i] = quants->y_zbin[q][1];
      quants->y_round[q][i] = quants->y_round[q][1];
      cpi->y_dequant[q][i] = cpi->y_dequant[q][1];

      quants->uv_quant[q][i] = quants->uv_quant[q][1];
      quants->uv_quant_fp[q][i] = quants->uv_quant_fp[q][1];
      quants->uv_round_fp[q][i] = quants->uv_round_fp[q][1];
      quants->uv_quant_shift[q][i] = quants->uv_quant_shift[q][1];
      quants->uv_zbin[q][i] = quants->uv_zbin[q][1];
      quants->uv_round[q][i] = quants->uv_round[q][1];
      cpi->uv_dequant[q][i] = cpi->uv_dequant[q][1];
    }
  }
}

void vp10_init_plane_quantizers(const VP10_COMP *cpi, MACROBLOCK *x,
                                const int segment_id) {
  const VP10_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  const QUANTS *const quants = &cpi->quants;
  const int qindex = vp10_get_qindex(&cm->seg, segment_id, cm->base_qindex);
  const int rdmult = vp10_compute_rd_mult(cpi, qindex + cm->y_dc_delta_q);
  int i;
#if CONFIG_NEW_QUANT
  int dq;
#endif

  // Y
  x->plane[0].quant = quants->y_quant[qindex];
  x->plane[0].quant_fp = quants->y_quant_fp[qindex];
  x->plane[0].round_fp = quants->y_round_fp[qindex];
  x->plane[0].quant_shift = quants->y_quant_shift[qindex];
  x->plane[0].zbin = quants->y_zbin[qindex];
  x->plane[0].round = quants->y_round[qindex];
  xd->plane[0].dequant = cpi->y_dequant[qindex];
#if CONFIG_NEW_QUANT
  for (dq = 0; dq < QUANT_PROFILES; dq ++) {
    x->plane[0].cuml_bins_nuq[dq] = (cuml_bins_type_nuq*)
                                    quants->y_cuml_bins_nuq[dq][qindex];
    xd->plane[0].dequant_val_nuq[dq] = (const dequant_val_type_nuq*)
                                   cpi->y_dequant_val_nuq[dq][qindex];
  }
#endif  // CONFIG_NEW_QUANT

  x->plane[0].quant_thred[0] = x->plane[0].zbin[0] * x->plane[0].zbin[0];
  x->plane[0].quant_thred[1] = x->plane[0].zbin[1] * x->plane[0].zbin[1];

  // UV
  for (i = 1; i < 3; i++) {
    x->plane[i].quant = quants->uv_quant[qindex];
    x->plane[i].quant_fp = quants->uv_quant_fp[qindex];
    x->plane[i].round_fp = quants->uv_round_fp[qindex];
    x->plane[i].quant_shift = quants->uv_quant_shift[qindex];
    x->plane[i].zbin = quants->uv_zbin[qindex];
    x->plane[i].round = quants->uv_round[qindex];
    xd->plane[i].dequant = cpi->uv_dequant[qindex];
#if CONFIG_NEW_QUANT
    for (dq = 0; dq < QUANT_PROFILES; dq ++) {
      x->plane[i].cuml_bins_nuq[dq] = (cuml_bins_type_nuq*)
                                      quants->uv_cuml_bins_nuq[dq][qindex];
      xd->plane[i].dequant_val_nuq[dq] = (const dequant_val_type_nuq*)
                                     cpi->uv_dequant_val_nuq[dq][qindex];
    }
#endif  // CONFIG_NEW_QUANT

    x->plane[i].quant_thred[0] = x->plane[i].zbin[0] * x->plane[i].zbin[0];
    x->plane[i].quant_thred[1] = x->plane[i].zbin[1] * x->plane[i].zbin[1];
  }

  x->skip_block = segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP);
  x->q_index = qindex;

  set_error_per_bit(x, rdmult);

  vp10_initialize_me_consts(cpi, x, x->q_index);
}

void vp10_frame_init_quantizer(VP10_COMP *cpi) {
  MACROBLOCK *const x = &cpi->td.mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  vp10_init_plane_quantizers(cpi, x, xd->mi[0]->mbmi.segment_id);
}

void vp10_set_quantizer(VP10_COMMON *cm, int q) {
  // quantizer has to be reinitialized with vp10_init_quantizer() if any
  // delta_q changes.
  cm->base_qindex = q;
  cm->y_dc_delta_q = 0;
  cm->uv_dc_delta_q = 0;
  cm->uv_ac_delta_q = 0;
}

// Table that converts 0-63 Q-range values passed in outside to the Qindex
// range used internally.
static const int quantizer_to_qindex[] = {
  0,    4,   8,  12,  16,  20,  24,  28,
  32,   36,  40,  44,  48,  52,  56,  60,
  64,   68,  72,  76,  80,  84,  88,  92,
  96,  100, 104, 108, 112, 116, 120, 124,
  128, 132, 136, 140, 144, 148, 152, 156,
  160, 164, 168, 172, 176, 180, 184, 188,
  192, 196, 200, 204, 208, 212, 216, 220,
  224, 228, 232, 236, 240, 244, 249, 255,
};

int vp10_quantizer_to_qindex(int quantizer) {
  return quantizer_to_qindex[quantizer];
}

int vp10_qindex_to_quantizer(int qindex) {
  int quantizer;

  for (quantizer = 0; quantizer < 64; ++quantizer)
    if (quantizer_to_qindex[quantizer] >= qindex)
      return quantizer;

  return 63;
}
