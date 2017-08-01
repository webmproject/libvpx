/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <arm_neon.h>

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/arm/mem_neon.h"

void vpx_quantize_b_neon(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                         int skip_block, const int16_t *zbin_ptr,
                         const int16_t *round_ptr, const int16_t *quant_ptr,
                         const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
                         tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr,
                         uint16_t *eob_ptr, const int16_t *scan_ptr,
                         const int16_t *iscan_ptr) {
  const int16x8_t zero = vdupq_n_s16(0);
  const int16x8_t one = vdupq_n_s16(1);
  const int16x8_t neg_one = vdupq_n_s16(-1);
  uint16x8_t eob_max;
  (void)scan_ptr;

  if (skip_block) {
    do {
      store_s16q_to_tran_low(qcoeff_ptr, zero);
      store_s16q_to_tran_low(dqcoeff_ptr, zero);
      qcoeff_ptr += 8;
      dqcoeff_ptr += 8;
      n_coeffs -= 8;
    } while (n_coeffs > 0);

    *eob_ptr = 0;
    return;
  }

  // Process first 8 values which include a dc component.
  {
    // Only the first element of each vector is DC.
    const int16x8_t zbin = vld1q_s16(zbin_ptr);
    const int16x8_t round = vld1q_s16(round_ptr);
    const int16x8_t quant = vld1q_s16(quant_ptr);
    const int16x8_t quant_shift = vld1q_s16(quant_shift_ptr);
    const int16x8_t dequant = vld1q_s16(dequant_ptr);
    // Add one because the eob does not index from 0.
    const uint16x8_t iscan =
        vreinterpretq_u16_s16(vaddq_s16(vld1q_s16(iscan_ptr), one));

    const int16x8_t coeff = load_tran_low_to_s16q(coeff_ptr);
    const int16x8_t coeff_sign = vshrq_n_s16(coeff, 15);
    const int16x8_t coeff_abs = vabsq_s16(coeff);

    const int16x8_t zbin_mask =
        vreinterpretq_s16_u16(vcgeq_s16(coeff_abs, zbin));

    const int16x8_t rounded = vqaddq_s16(coeff_abs, round);

    // (round * quant * 2) >> 16 >> 1 == (round * quant) >> 16
    int16x8_t qcoeff = vshrq_n_s16(vqdmulhq_s16(rounded, quant), 1);

    qcoeff = vaddq_s16(qcoeff, rounded);

    // (qcoeff * quant_shift * 2) >> 16 >> 1 == (qcoeff * quant_shift) >> 16
    qcoeff = vshrq_n_s16(vqdmulhq_s16(qcoeff, quant_shift), 1);

    // Restore the sign bit.
    qcoeff = veorq_s16(qcoeff, coeff_sign);
    qcoeff = vsubq_s16(qcoeff, coeff_sign);

    qcoeff = vandq_s16(qcoeff, zbin_mask);

    // Set non-zero elements to -1 and use that to extract values for eob.
    eob_max = vandq_u16(vtstq_s16(qcoeff, neg_one), iscan);

    coeff_ptr += 8;
    iscan_ptr += 8;

    store_s16q_to_tran_low(qcoeff_ptr, qcoeff);
    qcoeff_ptr += 8;

    qcoeff = vmulq_s16(qcoeff, dequant);

    store_s16q_to_tran_low(dqcoeff_ptr, qcoeff);
    dqcoeff_ptr += 8;
  }

  n_coeffs -= 8;

  {
    const int16x8_t zbin = vdupq_n_s16(zbin_ptr[1]);
    const int16x8_t round = vdupq_n_s16(round_ptr[1]);
    const int16x8_t quant = vdupq_n_s16(quant_ptr[1]);
    const int16x8_t quant_shift = vdupq_n_s16(quant_shift_ptr[1]);
    const int16x8_t dequant = vdupq_n_s16(dequant_ptr[1]);

    do {
      // Add one because the eob is not it's index.
      const uint16x8_t iscan =
          vreinterpretq_u16_s16(vaddq_s16(vld1q_s16(iscan_ptr), one));

      const int16x8_t coeff = load_tran_low_to_s16q(coeff_ptr);
      const int16x8_t coeff_sign = vshrq_n_s16(coeff, 15);
      const int16x8_t coeff_abs = vabsq_s16(coeff);

      const int16x8_t zbin_mask =
          vreinterpretq_s16_u16(vcgeq_s16(coeff_abs, zbin));

      const int16x8_t rounded = vqaddq_s16(coeff_abs, round);

      // (round * quant * 2) >> 16 >> 1 == (round * quant) >> 16
      int16x8_t qcoeff = vshrq_n_s16(vqdmulhq_s16(rounded, quant), 1);

      qcoeff = vaddq_s16(qcoeff, rounded);

      // (qcoeff * quant_shift * 2) >> 16 >> 1 == (qcoeff * quant_shift) >> 16
      qcoeff = vshrq_n_s16(vqdmulhq_s16(qcoeff, quant_shift), 1);

      // Restore the sign bit.
      qcoeff = veorq_s16(qcoeff, coeff_sign);
      qcoeff = vsubq_s16(qcoeff, coeff_sign);

      qcoeff = vandq_s16(qcoeff, zbin_mask);

      // Set non-zero elements to -1 and use that to extract values for eob.
      eob_max =
          vmaxq_u16(eob_max, vandq_u16(vtstq_s16(qcoeff, neg_one), iscan));

      coeff_ptr += 8;
      iscan_ptr += 8;

      store_s16q_to_tran_low(qcoeff_ptr, qcoeff);
      qcoeff_ptr += 8;

      qcoeff = vmulq_s16(qcoeff, dequant);

      store_s16q_to_tran_low(dqcoeff_ptr, qcoeff);
      dqcoeff_ptr += 8;

      n_coeffs -= 8;
    } while (n_coeffs > 0);
  }

  {
    const uint16x4_t eob_max_0 =
        vmax_u16(vget_low_u16(eob_max), vget_high_u16(eob_max));
    const uint16x4_t eob_max_1 = vpmax_u16(eob_max_0, eob_max_0);
    const uint16x4_t eob_max_2 = vpmax_u16(eob_max_1, eob_max_1);
    vst1_lane_u16(eob_ptr, eob_max_2, 0);
  }
}
