/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <tmmintrin.h>

#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/x86/bitdepth_conversion_sse2.h"

void vpx_quantize_b_ssse3(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                          int skip_block, const int16_t *zbin_ptr,
                          const int16_t *round_ptr, const int16_t *quant_ptr,
                          const int16_t *quant_shift_ptr,
                          tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                          const int16_t *dequant_ptr, uint16_t *eob_ptr,
                          const int16_t *scan_ptr, const int16_t *iscan_ptr) {
  const __m128i zero = _mm_setzero_si128();
  __m128i coeff0, coeff1;
  __m128i eob;
  __m128i zbin;
  __m128i round, quant, dequant, shift;
  intptr_t index = 0;
  (void)scan_ptr;

  if (skip_block) {
    do {
      store_tran_low(zero, dqcoeff_ptr + index);
      store_tran_low(zero, dqcoeff_ptr + index + 8);
      store_tran_low(zero, qcoeff_ptr + index);
      store_tran_low(zero, qcoeff_ptr + index + 8);
      index += 16;
    } while (index < n_coeffs);
    *eob_ptr = 0;
    return;
  }

  // Setup global values
  {
    const __m128i one = _mm_set1_epi16(1);
    zbin = _mm_load_si128((const __m128i *)zbin_ptr);
    // x86 has no "greater *or equal* comparison. Subtract 1 from zbin so
    // it is a strict "greater" comparison.
    zbin = _mm_sub_epi16(zbin, one);
    round = _mm_load_si128((const __m128i *)round_ptr);
    quant = _mm_load_si128((const __m128i *)quant_ptr);
    dequant = _mm_load_si128((const __m128i *)dequant_ptr);
    shift = _mm_load_si128((const __m128i *)quant_shift_ptr);
  }

  {
    __m128i qcoeff0, qcoeff1;
    __m128i qtmp0, qtmp1;
    __m128i cmp_mask0, cmp_mask1;
    __m128i zero_coeff0, zero_coeff1;
    __m128i iscan0, iscan1;
    __m128i eob1;

    // Do DC and first 15 AC
    coeff0 = load_tran_low(coeff_ptr + index);
    coeff1 = load_tran_low(coeff_ptr + index + 8);

    qcoeff0 = _mm_abs_epi16(coeff0);
    qcoeff1 = _mm_abs_epi16(coeff1);

    cmp_mask0 = _mm_cmpgt_epi16(qcoeff0, zbin);
    // Overwrite DC component.
    zbin = _mm_unpackhi_epi64(zbin, zbin);
    cmp_mask1 = _mm_cmpgt_epi16(qcoeff1, zbin);

    qcoeff0 = _mm_adds_epi16(qcoeff0, round);
    round = _mm_unpackhi_epi64(round, round);
    qcoeff1 = _mm_adds_epi16(qcoeff1, round);

    qtmp0 = _mm_mulhi_epi16(qcoeff0, quant);
    quant = _mm_unpackhi_epi64(quant, quant);
    qtmp1 = _mm_mulhi_epi16(qcoeff1, quant);

    qtmp0 = _mm_add_epi16(qtmp0, qcoeff0);
    qtmp1 = _mm_add_epi16(qtmp1, qcoeff1);

    qcoeff0 = _mm_mulhi_epi16(qtmp0, shift);
    shift = _mm_unpackhi_epi64(shift, shift);
    qcoeff1 = _mm_mulhi_epi16(qtmp1, shift);

    // Reinsert signs
    qcoeff0 = _mm_sign_epi16(qcoeff0, coeff0);
    qcoeff1 = _mm_sign_epi16(qcoeff1, coeff1);

    // Mask out zbin threshold coeffs
    qcoeff0 = _mm_and_si128(qcoeff0, cmp_mask0);
    qcoeff1 = _mm_and_si128(qcoeff1, cmp_mask1);

    store_tran_low(qcoeff0, qcoeff_ptr + index);
    store_tran_low(qcoeff1, qcoeff_ptr + index + 8);

    coeff0 = _mm_mullo_epi16(qcoeff0, dequant);
    dequant = _mm_unpackhi_epi64(dequant, dequant);
    coeff1 = _mm_mullo_epi16(qcoeff1, dequant);

    store_tran_low(coeff0, dqcoeff_ptr + index);
    store_tran_low(coeff1, dqcoeff_ptr + index + 8);

    // Scan for eob
    zero_coeff0 = _mm_cmpeq_epi16(coeff0, zero);
    zero_coeff1 = _mm_cmpeq_epi16(coeff1, zero);
    iscan0 = _mm_load_si128((const __m128i *)(iscan_ptr + index));
    iscan1 = _mm_load_si128((const __m128i *)(iscan_ptr + index + 8));
    // Add one to convert from indices to counts
    iscan0 = _mm_sub_epi16(iscan0, cmp_mask0);
    iscan1 = _mm_sub_epi16(iscan1, cmp_mask1);
    eob = _mm_andnot_si128(zero_coeff0, iscan0);
    eob1 = _mm_andnot_si128(zero_coeff1, iscan1);
    eob = _mm_max_epi16(eob, eob1);
  }
  index += 16;

  // AC only loop
  while (index < n_coeffs) {
    __m128i qcoeff0, qcoeff1;
    __m128i qtmp0, qtmp1;
    __m128i cmp_mask0, cmp_mask1;
    __m128i zero_coeff0, zero_coeff1;
    __m128i iscan0, iscan1;
    __m128i eob0, eob1;

    coeff0 = load_tran_low(coeff_ptr + index);
    coeff1 = load_tran_low(coeff_ptr + index + 8);

    qcoeff0 = _mm_abs_epi16(coeff0);
    qcoeff1 = _mm_abs_epi16(coeff1);

    cmp_mask0 = _mm_cmpgt_epi16(qcoeff0, zbin);
    cmp_mask1 = _mm_cmpgt_epi16(qcoeff1, zbin);

    qcoeff0 = _mm_adds_epi16(qcoeff0, round);
    qcoeff1 = _mm_adds_epi16(qcoeff1, round);

    qtmp0 = _mm_mulhi_epi16(qcoeff0, quant);
    qtmp1 = _mm_mulhi_epi16(qcoeff1, quant);

    qtmp0 = _mm_add_epi16(qtmp0, qcoeff0);
    qtmp1 = _mm_add_epi16(qtmp1, qcoeff1);

    qcoeff0 = _mm_mulhi_epi16(qtmp0, shift);
    qcoeff1 = _mm_mulhi_epi16(qtmp1, shift);

    // Reinsert signs
    qcoeff0 = _mm_sign_epi16(qcoeff0, coeff0);
    qcoeff1 = _mm_sign_epi16(qcoeff1, coeff1);

    // Mask out zbin threshold coeffs
    qcoeff0 = _mm_and_si128(qcoeff0, cmp_mask0);
    qcoeff1 = _mm_and_si128(qcoeff1, cmp_mask1);

    store_tran_low(qcoeff0, qcoeff_ptr + index);
    store_tran_low(qcoeff1, qcoeff_ptr + index + 8);

    coeff0 = _mm_mullo_epi16(qcoeff0, dequant);
    coeff1 = _mm_mullo_epi16(qcoeff1, dequant);

    store_tran_low(coeff0, dqcoeff_ptr + index);
    store_tran_low(coeff1, dqcoeff_ptr + index + 8);

    // Scan for eob
    zero_coeff0 = _mm_cmpeq_epi16(coeff0, zero);
    zero_coeff1 = _mm_cmpeq_epi16(coeff1, zero);
    iscan0 = _mm_load_si128((const __m128i *)(iscan_ptr + index));
    iscan1 = _mm_load_si128((const __m128i *)(iscan_ptr + index + 8));
    // Add one to convert from indices to counts
    iscan0 = _mm_sub_epi16(iscan0, cmp_mask0);
    iscan1 = _mm_sub_epi16(iscan1, cmp_mask1);
    eob0 = _mm_andnot_si128(zero_coeff0, iscan0);
    eob1 = _mm_andnot_si128(zero_coeff1, iscan1);
    eob0 = _mm_max_epi16(eob0, eob1);
    eob = _mm_max_epi16(eob, eob0);

    index += 16;
  }

  // Accumulate EOB
  {
    __m128i eob_shuffled;
    eob_shuffled = _mm_shuffle_epi32(eob, 0xe);
    eob = _mm_max_epi16(eob, eob_shuffled);
    eob_shuffled = _mm_shufflelo_epi16(eob, 0xe);
    eob = _mm_max_epi16(eob, eob_shuffled);
    eob_shuffled = _mm_shufflelo_epi16(eob, 0x1);
    eob = _mm_max_epi16(eob, eob_shuffled);
    *eob_ptr = _mm_extract_epi16(eob, 1);
  }
}
