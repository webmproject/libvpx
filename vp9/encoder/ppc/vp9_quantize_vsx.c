/*
 *  Copyright (c) 2018 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_config.h"

#include "./vp9_rtcd.h"
#include "vpx_dsp/ppc/types_vsx.h"

// Multiply the packed 16-bit integers in a and b, producing intermediate 32-bit
// integers, and return the high 16 bits of the intermediate integers.
// (a * b) >> 16
// Note: Because this is done in 2 operations, a and b cannot both be UINT16_MIN
static INLINE int16x8_t vec_mulhi(int16x8_t a, int16x8_t b) {
  // madds does ((A * B) >> 15) + C, we need >> 16, so we perform an extra right
  // shift.
  return vec_sra(vec_madds(a, b, vec_zeros_s16), vec_ones_u16);
}

// Negate 16-bit integers in a when the corresponding signed 16-bit
// integer in b is negative.
static INLINE int16x8_t vec_sign(int16x8_t a, int16x8_t b) {
  const int16x8_t mask = vec_sra(b, vec_shift_sign_s16);
  return vec_xor(vec_add(a, mask), mask);
}

// Compare packed 16-bit integers across a, and return the maximum value in
// every element. Returns a vector containing the biggest value across vector a.
static INLINE int16x8_t vec_max_across(int16x8_t a) {
  a = vec_max(a, vec_perm(a, a, vec_perm64));
  a = vec_max(a, vec_perm(a, a, vec_perm32));
  return vec_max(a, vec_perm(a, a, vec_perm16));
}

void vp9_quantize_fp_vsx(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
                         int skip_block, const int16_t *round_ptr,
                         const int16_t *quant_ptr, tran_low_t *qcoeff_ptr,
                         tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr,
                         uint16_t *eob_ptr, const int16_t *scan_ptr,
                         const int16_t *iscan_ptr) {
  int16x8_t qcoeff0, qcoeff1, dqcoeff0, dqcoeff1, eob;
  bool16x8_t zero_coeff0, zero_coeff1;

  int16x8_t round = vec_vsx_ld(0, round_ptr);
  int16x8_t quant = vec_vsx_ld(0, quant_ptr);
  int16x8_t dequant = vec_vsx_ld(0, dequant_ptr);
  int16x8_t coeff0 = vec_vsx_ld(0, coeff_ptr);
  int16x8_t coeff1 = vec_vsx_ld(16, coeff_ptr);
  int16x8_t scan0 = vec_vsx_ld(0, iscan_ptr);
  int16x8_t scan1 = vec_vsx_ld(16, iscan_ptr);

  (void)scan_ptr;
  (void)skip_block;
  assert(!skip_block);

  // First set of 8 coeff starts with DC + 7 AC
  qcoeff0 = vec_mulhi(vec_vaddshs(vec_abs(coeff0), round), quant);
  zero_coeff0 = vec_cmpeq(qcoeff0, vec_zeros_s16);
  qcoeff0 = vec_sign(qcoeff0, coeff0);
  vec_vsx_st(qcoeff0, 0, qcoeff_ptr);

  dqcoeff0 = vec_mladd(qcoeff0, dequant, vec_zeros_s16);
  vec_vsx_st(dqcoeff0, 0, dqcoeff_ptr);

  // Remove DC value from round and quant
  round = vec_splat(round, 1);
  quant = vec_splat(quant, 1);

  // Remove DC value from dequant
  dequant = vec_splat(dequant, 1);

  // Second set of 8 coeff starts with (all AC)
  qcoeff1 = vec_mulhi(vec_vaddshs(vec_abs(coeff1), round), quant);
  zero_coeff1 = vec_cmpeq(qcoeff1, vec_zeros_s16);
  qcoeff1 = vec_sign(qcoeff1, coeff1);
  vec_vsx_st(qcoeff1, 16, qcoeff_ptr);

  dqcoeff1 = vec_mladd(qcoeff1, dequant, vec_zeros_s16);
  vec_vsx_st(dqcoeff1, 16, dqcoeff_ptr);

  eob = vec_max(vec_or(scan0, zero_coeff0), vec_or(scan1, zero_coeff1));

  // We quantize 16 coeff up front (enough for a 4x4) and process 24 coeff per
  // loop iteration.
  // for 8x8: 16 + 2 x 24 = 64
  // for 16x16: 16 + 10 x 24 = 256
  if (n_coeffs > 16) {
    int16x8_t coeff2, qcoeff2, dqcoeff2, eob2, scan2;
    bool16x8_t zero_coeff2;

    int index = 16;
    int off0 = 32;
    int off1 = 48;
    int off2 = 64;

    do {
      coeff0 = vec_vsx_ld(off0, coeff_ptr);
      coeff1 = vec_vsx_ld(off1, coeff_ptr);
      coeff2 = vec_vsx_ld(off2, coeff_ptr);
      scan0 = vec_vsx_ld(off0, iscan_ptr);
      scan1 = vec_vsx_ld(off1, iscan_ptr);
      scan2 = vec_vsx_ld(off2, iscan_ptr);

      qcoeff0 = vec_mulhi(vec_vaddshs(vec_abs(coeff0), round), quant);
      zero_coeff0 = vec_cmpeq(qcoeff0, vec_zeros_s16);
      qcoeff0 = vec_sign(qcoeff0, coeff0);
      vec_vsx_st(qcoeff0, off0, qcoeff_ptr);
      dqcoeff0 = vec_mladd(qcoeff0, dequant, vec_zeros_s16);
      vec_vsx_st(dqcoeff0, off0, dqcoeff_ptr);

      qcoeff1 = vec_mulhi(vec_vaddshs(vec_abs(coeff1), round), quant);
      zero_coeff1 = vec_cmpeq(qcoeff1, vec_zeros_s16);
      qcoeff1 = vec_sign(qcoeff1, coeff1);
      vec_vsx_st(qcoeff1, off1, qcoeff_ptr);
      dqcoeff1 = vec_mladd(qcoeff1, dequant, vec_zeros_s16);
      vec_vsx_st(dqcoeff1, off1, dqcoeff_ptr);

      qcoeff2 = vec_mulhi(vec_vaddshs(vec_abs(coeff2), round), quant);
      zero_coeff2 = vec_cmpeq(qcoeff2, vec_zeros_s16);
      qcoeff2 = vec_sign(qcoeff2, coeff2);
      vec_vsx_st(qcoeff2, off2, qcoeff_ptr);
      dqcoeff2 = vec_mladd(qcoeff2, dequant, vec_zeros_s16);
      vec_vsx_st(dqcoeff2, off2, dqcoeff_ptr);

      eob = vec_max(eob, vec_or(scan0, zero_coeff0));
      eob2 = vec_max(vec_or(scan1, zero_coeff1), vec_or(scan2, zero_coeff2));
      eob = vec_max(eob, eob2);

      index += 24;
      off0 += 48;
      off1 += 48;
      off2 += 48;
    } while (index < n_coeffs);
  }

  eob = vec_max_across(eob);
  *eob_ptr = eob[0] + 1;
}
