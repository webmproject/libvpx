/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <emmintrin.h>
#include <xmmintrin.h>

#include "./vp9_rtcd.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_dsp/x86/bitdepth_conversion_sse2.h"
#include "vpx_dsp/x86/quantize_sse2.h"

void vp9_quantize_fp_sse2(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                          const int16_t *round_ptr, const int16_t *quant_ptr,
                          tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                          const int16_t *dequant_ptr, uint16_t *eob_ptr,
                          const int16_t *scan, const int16_t *iscan) {
  const __m128i zero = _mm_setzero_si128();
  __m128i thr;
  int nzflag;
  __m128i eob;
  __m128i round, quant, dequant;

  (void)scan;

  coeff_ptr += n_coeffs;
  iscan += n_coeffs;
  qcoeff_ptr += n_coeffs;
  dqcoeff_ptr += n_coeffs;
  n_coeffs = -n_coeffs;

  // Setup global values.
  load_fp_values(round_ptr, &round, quant_ptr, &quant, dequant_ptr, &dequant);

  {
    __m128i coeff0, coeff1;
    __m128i coeff0_sign, coeff1_sign;
    __m128i qcoeff0, qcoeff1;
    // Do DC and first 15 AC.
    coeff0 = load_tran_low(coeff_ptr + n_coeffs);
    coeff1 = load_tran_low(coeff_ptr + n_coeffs + 8);

    // Poor man's abs().
    coeff0_sign = _mm_srai_epi16(coeff0, 15);
    coeff1_sign = _mm_srai_epi16(coeff1, 15);
    qcoeff0 = invert_sign_sse2(coeff0, coeff0_sign);
    qcoeff1 = invert_sign_sse2(coeff1, coeff1_sign);

    qcoeff0 = _mm_adds_epi16(qcoeff0, round);
    qcoeff0 = _mm_mulhi_epi16(qcoeff0, quant);

    round = _mm_unpackhi_epi64(round, round);
    quant = _mm_unpackhi_epi64(quant, quant);

    qcoeff1 = _mm_adds_epi16(qcoeff1, round);
    qcoeff1 = _mm_mulhi_epi16(qcoeff1, quant);

    // Reinsert signs.
    qcoeff0 = invert_sign_sse2(qcoeff0, coeff0_sign);
    qcoeff1 = invert_sign_sse2(qcoeff1, coeff1_sign);

    store_tran_low(qcoeff0, qcoeff_ptr + n_coeffs);
    store_tran_low(qcoeff1, qcoeff_ptr + n_coeffs + 8);

    qcoeff0 = _mm_mullo_epi16(qcoeff0, dequant);
    dequant = _mm_unpackhi_epi64(dequant, dequant);
    qcoeff1 = _mm_mullo_epi16(qcoeff1, dequant);

    store_tran_low(qcoeff0, dqcoeff_ptr + n_coeffs);
    store_tran_low(qcoeff1, dqcoeff_ptr + n_coeffs + 8);

    eob = scan_for_eob(&qcoeff0, &qcoeff1, iscan + n_coeffs, 0, zero);

    n_coeffs += 8 * 2;
  }

  thr = _mm_srai_epi16(dequant, 1);

  // AC only loop.
  while (n_coeffs < 0) {
    __m128i coeff0, coeff1;
    __m128i coeff0_sign, coeff1_sign;
    __m128i qcoeff0, qcoeff1;

    coeff0 = load_tran_low(coeff_ptr + n_coeffs);
    coeff1 = load_tran_low(coeff_ptr + n_coeffs + 8);

    // Poor man's abs().
    coeff0_sign = _mm_srai_epi16(coeff0, 15);
    coeff1_sign = _mm_srai_epi16(coeff1, 15);
    qcoeff0 = invert_sign_sse2(coeff0, coeff0_sign);
    qcoeff1 = invert_sign_sse2(coeff1, coeff1_sign);

    nzflag = _mm_movemask_epi8(_mm_cmpgt_epi16(qcoeff0, thr)) |
             _mm_movemask_epi8(_mm_cmpgt_epi16(qcoeff1, thr));

    if (nzflag) {
      qcoeff0 = _mm_adds_epi16(qcoeff0, round);
      qcoeff1 = _mm_adds_epi16(qcoeff1, round);
      qcoeff0 = _mm_mulhi_epi16(qcoeff0, quant);
      qcoeff1 = _mm_mulhi_epi16(qcoeff1, quant);

      // Reinsert signs.
      qcoeff0 = invert_sign_sse2(qcoeff0, coeff0_sign);
      qcoeff1 = invert_sign_sse2(qcoeff1, coeff1_sign);

      store_tran_low(qcoeff0, qcoeff_ptr + n_coeffs);
      store_tran_low(qcoeff1, qcoeff_ptr + n_coeffs + 8);

      coeff0 = _mm_mullo_epi16(qcoeff0, dequant);
      coeff1 = _mm_mullo_epi16(qcoeff1, dequant);

      store_tran_low(coeff0, dqcoeff_ptr + n_coeffs);
      store_tran_low(coeff1, dqcoeff_ptr + n_coeffs + 8);
    } else {
      store_zero_tran_low(qcoeff_ptr + n_coeffs);
      store_zero_tran_low(qcoeff_ptr + n_coeffs + 8);

      store_zero_tran_low(dqcoeff_ptr + n_coeffs);
      store_zero_tran_low(dqcoeff_ptr + n_coeffs + 8);
    }

    if (nzflag) {
      const __m128i eob0 =
          scan_for_eob(&coeff0, &coeff1, iscan + n_coeffs, 0, zero);
      eob = _mm_max_epi16(eob, eob0);
    }
    n_coeffs += 8 * 2;
  }

  *eob_ptr = accumulate_eob(eob);
}
