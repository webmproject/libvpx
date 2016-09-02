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

#include "aom_dsp/quantize.h"
#include "aom_mem/aom_mem.h"

#if CONFIG_AOM_QM
void aom_quantize_dc(const tran_low_t *coeff_ptr, int n_coeffs, int skip_block,
                     const int16_t *round_ptr, const int16_t quant,
                     tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                     const int16_t dequant_ptr, uint16_t *eob_ptr,
                     const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr) {
  const int rc = 0;
  const int coeff = coeff_ptr[rc];
  const int coeff_sign = (coeff >> 31);
  const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
  int64_t tmp, eob = -1;
  int32_t tmp32;
  int dequant =
      (dequant_ptr * iqm_ptr[rc] + (1 << (AOM_QM_BITS - 1))) >> AOM_QM_BITS;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    tmp = clamp(abs_coeff + round_ptr[rc != 0], INT16_MIN, INT16_MAX);
    tmp32 = (int32_t)((tmp * qm_ptr[rc] * quant) >> (16 + AOM_QM_BITS));
    qcoeff_ptr[rc] = (tmp32 ^ coeff_sign) - coeff_sign;
    dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant;
    if (tmp32) eob = 0;
  }
  *eob_ptr = eob + 1;
}

#if CONFIG_AOM_HIGHBITDEPTH
void aom_highbd_quantize_dc(const tran_low_t *coeff_ptr, int n_coeffs,
                            int skip_block, const int16_t *round_ptr,
                            const int16_t quant, tran_low_t *qcoeff_ptr,
                            tran_low_t *dqcoeff_ptr, const int16_t dequant_ptr,
                            uint16_t *eob_ptr, const qm_val_t *qm_ptr,
                            const qm_val_t *iqm_ptr) {
  int eob = -1;
  int dequant =
      (dequant_ptr * iqm_ptr[0] + (1 << (AOM_QM_BITS - 1))) >> AOM_QM_BITS;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    const int coeff = coeff_ptr[0];
    const int coeff_sign = (coeff >> 31);
    const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
    const int64_t tmp = abs_coeff + round_ptr[0];
    const uint32_t abs_qcoeff =
        (uint32_t)((tmp * qm_ptr[0] * quant) >> (16 + AOM_QM_BITS));
    qcoeff_ptr[0] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
    dqcoeff_ptr[0] = qcoeff_ptr[0] * dequant;
    if (abs_qcoeff) eob = 0;
  }
  *eob_ptr = eob + 1;
}
#endif

void aom_quantize_dc_32x32(const tran_low_t *coeff_ptr, int skip_block,
                           const int16_t *round_ptr, const int16_t quant,
                           tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                           const int16_t dequant_ptr, uint16_t *eob_ptr,
                           const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr) {
  const int n_coeffs = 1024;
  const int rc = 0;
  const int coeff = coeff_ptr[rc];
  const int coeff_sign = (coeff >> 31);
  const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
  int64_t tmp, eob = -1;
  int32_t tmp32;
  int dequant;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    tmp = clamp(abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], 1),
                INT16_MIN, INT16_MAX);
    tmp32 = (int32_t)((tmp * qm_ptr[rc] * quant) >> (15 + AOM_QM_BITS));
    qcoeff_ptr[rc] = (tmp32 ^ coeff_sign) - coeff_sign;
    dequant =
        (dequant_ptr * iqm_ptr[rc] + (1 << (AOM_QM_BITS - 1))) >> AOM_QM_BITS;
    dqcoeff_ptr[rc] = (qcoeff_ptr[rc] * dequant) / 2;
    if (tmp32) eob = 0;
  }
  *eob_ptr = eob + 1;
}

#if CONFIG_AOM_HIGHBITDEPTH
void aom_highbd_quantize_dc_32x32(const tran_low_t *coeff_ptr, int skip_block,
                                  const int16_t *round_ptr, const int16_t quant,
                                  tran_low_t *qcoeff_ptr,
                                  tran_low_t *dqcoeff_ptr,
                                  const int16_t dequant_ptr, uint16_t *eob_ptr,
                                  const qm_val_t *qm_ptr,
                                  const qm_val_t *iqm_ptr) {
  const int n_coeffs = 1024;
  int eob = -1;
  int dequant;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    const int coeff = coeff_ptr[0];
    const int coeff_sign = (coeff >> 31);
    const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
    const int64_t tmp = abs_coeff + ROUND_POWER_OF_TWO(round_ptr[0], 1);
    const uint32_t abs_qcoeff =
        (uint32_t)((tmp * qm_ptr[0] * quant) >> (15 + AOM_QM_BITS));
    qcoeff_ptr[0] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
    dequant =
        (dequant_ptr * iqm_ptr[0] + (1 << (AOM_QM_BITS - 1))) >> AOM_QM_BITS;
    dqcoeff_ptr[0] = (qcoeff_ptr[0] * dequant) / 2;
    if (abs_qcoeff) eob = 0;
  }
  *eob_ptr = eob + 1;
}
#endif

void aom_quantize_b_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                      int skip_block, const int16_t *zbin_ptr,
                      const int16_t *round_ptr, const int16_t *quant_ptr,
                      const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
                      tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr,
                      uint16_t *eob_ptr, const int16_t *scan,
                      const int16_t *iscan, const qm_val_t *qm_ptr,
                      const qm_val_t *iqm_ptr) {
  int i, non_zero_count = (int)n_coeffs, eob = -1;
  const int zbins[2] = { zbin_ptr[0], zbin_ptr[1] };
  const int nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };
  (void)iscan;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    // Pre-scan pass
    for (i = (int)n_coeffs - 1; i >= 0; i--) {
      const int rc = scan[i];
      const qm_val_t wt = qm_ptr[rc];
      const int coeff = coeff_ptr[rc] * wt;

      if (coeff < (zbins[rc != 0] << AOM_QM_BITS) &&
          coeff > (nzbins[rc != 0] << AOM_QM_BITS))
        non_zero_count--;
      else
        break;
    }

    // Quantization pass: All coefficients with index >= zero_flag are
    // skippable. Note: zero_flag can be zero.
    for (i = 0; i < non_zero_count; i++) {
      const int rc = scan[i];
      const qm_val_t wt = qm_ptr[rc];
      const int coeff = coeff_ptr[rc];
      const int coeff_sign = (coeff >> 31);
      const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
      int dequant;

      if (abs_coeff * wt >= (zbins[rc != 0] << AOM_QM_BITS)) {
        int32_t tmp32;
        int64_t tmp =
            clamp(abs_coeff + round_ptr[rc != 0], INT16_MIN, INT16_MAX);
        tmp = tmp * wt;
        tmp32 = ((((tmp * quant_ptr[rc != 0]) >> 16) + tmp) *
                 quant_shift_ptr[rc != 0]) >>
                (16 + AOM_QM_BITS);  // quantization
        dequant =
            (dequant_ptr[rc != 0] * iqm_ptr[rc] + (1 << (AOM_QM_BITS - 1))) >>
            AOM_QM_BITS;
        qcoeff_ptr[rc] = (tmp32 ^ coeff_sign) - coeff_sign;
        dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant;

        if (tmp32) eob = i;
      }
    }
  }
  *eob_ptr = eob + 1;
}

#if CONFIG_AOM_HIGHBITDEPTH
void aom_highbd_quantize_b_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                             int skip_block, const int16_t *zbin_ptr,
                             const int16_t *round_ptr, const int16_t *quant_ptr,
                             const int16_t *quant_shift_ptr,
                             tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                             const int16_t *dequant_ptr, uint16_t *eob_ptr,
                             const int16_t *scan, const int16_t *iscan,
                             const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr) {
  int i, non_zero_count = (int)n_coeffs, eob = -1;
  const int zbins[2] = { zbin_ptr[0], zbin_ptr[1] };
  const int nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };
  int dequant;
  (void)iscan;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    // Pre-scan pass
    for (i = (int)n_coeffs - 1; i >= 0; i--) {
      const int rc = scan[i];
      const qm_val_t wt = qm_ptr[rc];
      const int coeff = coeff_ptr[rc] * wt;

      if (coeff < (zbins[rc != 0] << AOM_QM_BITS) &&
          coeff > (nzbins[rc != 0] << AOM_QM_BITS))
        non_zero_count--;
      else
        break;
    }

    // Quantization pass: All coefficients with index >= zero_flag are
    // skippable. Note: zero_flag can be zero.
    for (i = 0; i < non_zero_count; i++) {
      const int rc = scan[i];
      const int coeff = coeff_ptr[rc];
      const qm_val_t wt = qm_ptr[rc];
      const int coeff_sign = (coeff >> 31);
      const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;

      if (abs_coeff * wt >= (zbins[rc != 0] << AOM_QM_BITS)) {
        const int64_t tmp1 = abs_coeff + round_ptr[rc != 0];
        const int64_t tmpw = tmp1 * wt;
        const int64_t tmp2 = ((tmpw * quant_ptr[rc != 0]) >> 16) + tmpw;
        const uint32_t abs_qcoeff =
            (uint32_t)((tmp2 * quant_shift_ptr[rc != 0]) >> (16 + AOM_QM_BITS));
        qcoeff_ptr[rc] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
        dequant =
            (dequant_ptr[rc != 0] * iqm_ptr[rc] + (1 << (AOM_QM_BITS - 1))) >>
            AOM_QM_BITS;
        dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant;
        if (abs_qcoeff) eob = i;
      }
    }
  }
  *eob_ptr = eob + 1;
}
#endif

void aom_quantize_b_32x32_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                            int skip_block, const int16_t *zbin_ptr,
                            const int16_t *round_ptr, const int16_t *quant_ptr,
                            const int16_t *quant_shift_ptr,
                            tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                            const int16_t *dequant_ptr, uint16_t *eob_ptr,
                            const int16_t *scan, const int16_t *iscan,
                            const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr) {
  const int zbins[2] = { ROUND_POWER_OF_TWO(zbin_ptr[0], 1),
                         ROUND_POWER_OF_TWO(zbin_ptr[1], 1) };
  const int nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };

  int idx = 0;
  int idx_arr[1024];
  int i, eob = -1;
  int dequant;
  (void)iscan;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    // Pre-scan pass
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      const qm_val_t wt = qm_ptr[rc];
      const int coeff = coeff_ptr[rc] * wt;

      // If the coefficient is out of the base ZBIN range, keep it for
      // quantization.
      if (coeff >= (zbins[rc != 0] << AOM_QM_BITS) ||
          coeff <= (nzbins[rc != 0] << AOM_QM_BITS))
        idx_arr[idx++] = i;
    }

    // Quantization pass: only process the coefficients selected in
    // pre-scan pass. Note: idx can be zero.
    for (i = 0; i < idx; i++) {
      const int rc = scan[idx_arr[i]];
      const int coeff = coeff_ptr[rc];
      const int coeff_sign = (coeff >> 31);
      const qm_val_t wt = qm_ptr[rc];
      int64_t tmp;
      int tmp32;
      int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
      abs_coeff += ROUND_POWER_OF_TWO(round_ptr[rc != 0], 1);
      tmp = clamp(abs_coeff, INT16_MIN, INT16_MAX);
      tmp = tmp * wt;
      tmp32 = ((((tmp * quant_ptr[rc != 0]) >> 16) + tmp) *
               quant_shift_ptr[rc != 0]) >>
              (15 + AOM_QM_BITS);

      qcoeff_ptr[rc] = (tmp32 ^ coeff_sign) - coeff_sign;
      dequant =
          (dequant_ptr[rc != 0] * iqm_ptr[rc] + (1 << (AOM_QM_BITS - 1))) >>
          AOM_QM_BITS;
      dqcoeff_ptr[rc] = (qcoeff_ptr[rc] * dequant) / 2;

      if (tmp32) eob = idx_arr[i];
    }
  }
  *eob_ptr = eob + 1;
}

#if CONFIG_AOM_HIGHBITDEPTH
void aom_highbd_quantize_b_32x32_c(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block,
    const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan, const qm_val_t *qm_ptr,
    const qm_val_t *iqm_ptr) {
  const int zbins[2] = { ROUND_POWER_OF_TWO(zbin_ptr[0], 1),
                         ROUND_POWER_OF_TWO(zbin_ptr[1], 1) };
  const int nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };

  int idx = 0;
  int idx_arr[1024];
  int i, eob = -1;
  int dequant;
  (void)iscan;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    // Pre-scan pass
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      const qm_val_t wt = qm_ptr[rc];
      const int coeff = coeff_ptr[rc] * wt;

      // If the coefficient is out of the base ZBIN range, keep it for
      // quantization.
      if (coeff >= (zbins[rc != 0] << AOM_QM_BITS) ||
          coeff <= (nzbins[rc != 0] << AOM_QM_BITS))
        idx_arr[idx++] = i;
    }

    // Quantization pass: only process the coefficients selected in
    // pre-scan pass. Note: idx can be zero.
    for (i = 0; i < idx; i++) {
      const int rc = scan[idx_arr[i]];
      const int coeff = coeff_ptr[rc];
      const int coeff_sign = (coeff >> 31);
      const qm_val_t wt = qm_ptr[rc];
      const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
      const int64_t tmp1 =
          abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], 1);
      const int64_t tmpw = tmp1 * wt;
      const int64_t tmp2 = ((tmpw * quant_ptr[rc != 0]) >> 16) + tmpw;
      const uint32_t abs_qcoeff =
          (uint32_t)((tmp2 * quant_shift_ptr[rc != 0]) >> (15 + AOM_QM_BITS));
      qcoeff_ptr[rc] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
      dequant =
          (dequant_ptr[rc != 0] * iqm_ptr[rc] + (1 << (AOM_QM_BITS - 1))) >>
          AOM_QM_BITS;
      dqcoeff_ptr[rc] = (qcoeff_ptr[rc] * dequant) / 2;
      if (abs_qcoeff) eob = idx_arr[i];
    }
  }
  *eob_ptr = eob + 1;
}
#endif
#else
void aom_quantize_dc(const tran_low_t *coeff_ptr, int n_coeffs, int skip_block,
                     const int16_t *round_ptr, const int16_t quant,
                     tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                     const int16_t dequant_ptr, uint16_t *eob_ptr) {
  const int rc = 0;
  const int coeff = coeff_ptr[rc];
  const int coeff_sign = (coeff >> 31);
  const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
  int tmp, eob = -1;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    tmp = clamp(abs_coeff + round_ptr[rc != 0], INT16_MIN, INT16_MAX);
    tmp = (tmp * quant) >> 16;
    qcoeff_ptr[rc] = (tmp ^ coeff_sign) - coeff_sign;
    dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant_ptr;
    if (tmp) eob = 0;
  }
  *eob_ptr = eob + 1;
}

#if CONFIG_AOM_HIGHBITDEPTH
void aom_highbd_quantize_dc(const tran_low_t *coeff_ptr, int n_coeffs,
                            int skip_block, const int16_t *round_ptr,
                            const int16_t quant, tran_low_t *qcoeff_ptr,
                            tran_low_t *dqcoeff_ptr, const int16_t dequant_ptr,
                            uint16_t *eob_ptr) {
  int eob = -1;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    const int coeff = coeff_ptr[0];
    const int coeff_sign = (coeff >> 31);
    const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
    const int64_t tmp = abs_coeff + round_ptr[0];
    const uint32_t abs_qcoeff = (uint32_t)((tmp * quant) >> 16);
    qcoeff_ptr[0] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
    dqcoeff_ptr[0] = qcoeff_ptr[0] * dequant_ptr;
    if (abs_qcoeff) eob = 0;
  }
  *eob_ptr = eob + 1;
}
#endif

void aom_quantize_dc_32x32(const tran_low_t *coeff_ptr, int skip_block,
                           const int16_t *round_ptr, const int16_t quant,
                           tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                           const int16_t dequant_ptr, uint16_t *eob_ptr) {
  const int n_coeffs = 1024;
  const int rc = 0;
  const int coeff = coeff_ptr[rc];
  const int coeff_sign = (coeff >> 31);
  const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
  int tmp, eob = -1;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    tmp = clamp(abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], 1),
                INT16_MIN, INT16_MAX);
    tmp = (tmp * quant) >> 15;
    qcoeff_ptr[rc] = (tmp ^ coeff_sign) - coeff_sign;
    dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant_ptr / 2;
    if (tmp) eob = 0;
  }
  *eob_ptr = eob + 1;
}

#if CONFIG_AOM_HIGHBITDEPTH
void aom_highbd_quantize_dc_32x32(const tran_low_t *coeff_ptr, int skip_block,
                                  const int16_t *round_ptr, const int16_t quant,
                                  tran_low_t *qcoeff_ptr,
                                  tran_low_t *dqcoeff_ptr,
                                  const int16_t dequant_ptr,
                                  uint16_t *eob_ptr) {
  const int n_coeffs = 1024;
  int eob = -1;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    const int coeff = coeff_ptr[0];
    const int coeff_sign = (coeff >> 31);
    const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
    const int64_t tmp = abs_coeff + ROUND_POWER_OF_TWO(round_ptr[0], 1);
    const uint32_t abs_qcoeff = (uint32_t)((tmp * quant) >> 15);
    qcoeff_ptr[0] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
    dqcoeff_ptr[0] = qcoeff_ptr[0] * dequant_ptr / 2;
    if (abs_qcoeff) eob = 0;
  }
  *eob_ptr = eob + 1;
}
#endif

void aom_quantize_b_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                      int skip_block, const int16_t *zbin_ptr,
                      const int16_t *round_ptr, const int16_t *quant_ptr,
                      const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
                      tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr,
                      uint16_t *eob_ptr, const int16_t *scan,
                      const int16_t *iscan) {
  int i, non_zero_count = (int)n_coeffs, eob = -1;
  const int zbins[2] = { zbin_ptr[0], zbin_ptr[1] };
  const int nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };
  (void)iscan;

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
        int tmp = clamp(abs_coeff + round_ptr[rc != 0], INT16_MIN, INT16_MAX);
        tmp = ((((tmp * quant_ptr[rc != 0]) >> 16) + tmp) *
               quant_shift_ptr[rc != 0]) >>
              16;  // quantization
        qcoeff_ptr[rc] = (tmp ^ coeff_sign) - coeff_sign;
        dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant_ptr[rc != 0];

        if (tmp) eob = i;
      }
    }
  }
  *eob_ptr = eob + 1;
}

#if CONFIG_AOM_HIGHBITDEPTH
void aom_highbd_quantize_b_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                             int skip_block, const int16_t *zbin_ptr,
                             const int16_t *round_ptr, const int16_t *quant_ptr,
                             const int16_t *quant_shift_ptr,
                             tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                             const int16_t *dequant_ptr, uint16_t *eob_ptr,
                             const int16_t *scan, const int16_t *iscan) {
  int i, non_zero_count = (int)n_coeffs, eob = -1;
  const int zbins[2] = { zbin_ptr[0], zbin_ptr[1] };
  const int nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };
  (void)iscan;

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
        const int64_t tmp1 = abs_coeff + round_ptr[rc != 0];
        const int64_t tmp2 = ((tmp1 * quant_ptr[rc != 0]) >> 16) + tmp1;
        const uint32_t abs_qcoeff =
            (uint32_t)((tmp2 * quant_shift_ptr[rc != 0]) >> 16);
        qcoeff_ptr[rc] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
        dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant_ptr[rc != 0];
        if (abs_qcoeff) eob = i;
      }
    }
  }
  *eob_ptr = eob + 1;
}
#endif

void aom_quantize_b_32x32_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                            int skip_block, const int16_t *zbin_ptr,
                            const int16_t *round_ptr, const int16_t *quant_ptr,
                            const int16_t *quant_shift_ptr,
                            tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                            const int16_t *dequant_ptr, uint16_t *eob_ptr,
                            const int16_t *scan, const int16_t *iscan) {
  const int zbins[2] = { ROUND_POWER_OF_TWO(zbin_ptr[0], 1),
                         ROUND_POWER_OF_TWO(zbin_ptr[1], 1) };
  const int nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };

  int idx = 0;
  int idx_arr[1024];
  int i, eob = -1;
  (void)iscan;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    // Pre-scan pass
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      const int coeff = coeff_ptr[rc];

      // If the coefficient is out of the base ZBIN range, keep it for
      // quantization.
      if (coeff >= zbins[rc != 0] || coeff <= nzbins[rc != 0])
        idx_arr[idx++] = i;
    }

    // Quantization pass: only process the coefficients selected in
    // pre-scan pass. Note: idx can be zero.
    for (i = 0; i < idx; i++) {
      const int rc = scan[idx_arr[i]];
      const int coeff = coeff_ptr[rc];
      const int coeff_sign = (coeff >> 31);
      int tmp;
      int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
      abs_coeff += ROUND_POWER_OF_TWO(round_ptr[rc != 0], 1);
      abs_coeff = clamp(abs_coeff, INT16_MIN, INT16_MAX);
      tmp = ((((abs_coeff * quant_ptr[rc != 0]) >> 16) + abs_coeff) *
             quant_shift_ptr[rc != 0]) >>
            15;

      qcoeff_ptr[rc] = (tmp ^ coeff_sign) - coeff_sign;
      dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant_ptr[rc != 0] / 2;

      if (tmp) eob = idx_arr[i];
    }
  }
  *eob_ptr = eob + 1;
}

#if CONFIG_AOM_HIGHBITDEPTH
void aom_highbd_quantize_b_32x32_c(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block,
    const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan) {
  const int zbins[2] = { ROUND_POWER_OF_TWO(zbin_ptr[0], 1),
                         ROUND_POWER_OF_TWO(zbin_ptr[1], 1) };
  const int nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };

  int idx = 0;
  int idx_arr[1024];
  int i, eob = -1;
  (void)iscan;

  memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
  memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

  if (!skip_block) {
    // Pre-scan pass
    for (i = 0; i < n_coeffs; i++) {
      const int rc = scan[i];
      const int coeff = coeff_ptr[rc];

      // If the coefficient is out of the base ZBIN range, keep it for
      // quantization.
      if (coeff >= zbins[rc != 0] || coeff <= nzbins[rc != 0])
        idx_arr[idx++] = i;
    }

    // Quantization pass: only process the coefficients selected in
    // pre-scan pass. Note: idx can be zero.
    for (i = 0; i < idx; i++) {
      const int rc = scan[idx_arr[i]];
      const int coeff = coeff_ptr[rc];
      const int coeff_sign = (coeff >> 31);
      const int abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
      const int64_t tmp1 =
          abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], 1);
      const int64_t tmp2 = ((tmp1 * quant_ptr[rc != 0]) >> 16) + tmp1;
      const uint32_t abs_qcoeff =
          (uint32_t)((tmp2 * quant_shift_ptr[rc != 0]) >> 15);
      qcoeff_ptr[rc] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
      dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant_ptr[rc != 0] / 2;
      if (abs_qcoeff) eob = idx_arr[i];
    }
  }
  *eob_ptr = eob + 1;
}
#endif
#endif
