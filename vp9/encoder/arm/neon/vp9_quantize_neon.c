/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <arm_neon.h>
#include <assert.h>
#include <math.h>

#include "./vpx_config.h"
#include "vpx_mem/vpx_mem.h"

#include "vp9/common/vp9_quant_common.h"
#include "vp9/common/vp9_seg_common.h"

#include "vp9/encoder/vp9_encoder.h"
#include "vp9/encoder/vp9_quantize.h"
#include "vp9/encoder/vp9_rd.h"

#include "vpx_dsp/arm/idct_neon.h"
#include "vpx_dsp/arm/mem_neon.h"
#include "vpx_dsp/vpx_dsp_common.h"

static VPX_FORCE_INLINE void calculate_dqcoeff_and_store(
    const int16x8_t qcoeff, const int16x8_t dequant, tran_low_t *dqcoeff) {
  const int32x4_t dqcoeff_0 =
      vmull_s16(vget_low_s16(qcoeff), vget_low_s16(dequant));
  const int32x4_t dqcoeff_1 =
      vmull_s16(vget_high_s16(qcoeff), vget_high_s16(dequant));

#if CONFIG_VP9_HIGHBITDEPTH
  vst1q_s32(dqcoeff, dqcoeff_0);
  vst1q_s32(dqcoeff + 4, dqcoeff_1);
#else
  vst1q_s16(dqcoeff, vcombine_s16(vmovn_s32(dqcoeff_0), vmovn_s32(dqcoeff_1)));
#endif  // CONFIG_VP9_HIGHBITDEPTH
}

static VPX_FORCE_INLINE int16x8_t get_max_lane_eob(const int16_t *iscan_ptr,
                                                   int16x8_t v_eobmax,
                                                   uint16x8_t v_nz_mask) {
  const int16x8_t v_iscan = vld1q_s16(&iscan_ptr[0]);
  const int16x8_t v_iscan_plus1 = vaddq_s16(v_iscan, vdupq_n_s16(1));
  const int16x8_t v_nz_iscan =
      vbslq_s16(v_nz_mask, vdupq_n_s16(0), v_iscan_plus1);
  return vmaxq_s16(v_eobmax, v_nz_iscan);
}

static VPX_FORCE_INLINE uint16_t get_max_eob(int16x8_t v_eobmax) {
#ifdef __aarch64__
  return (uint16_t)vmaxvq_s16(v_eobmax);
#else
  const int16x4_t v_eobmax_3210 =
      vmax_s16(vget_low_s16(v_eobmax), vget_high_s16(v_eobmax));
  const int64x1_t v_eobmax_xx32 =
      vshr_n_s64(vreinterpret_s64_s16(v_eobmax_3210), 32);
  const int16x4_t v_eobmax_tmp =
      vmax_s16(v_eobmax_3210, vreinterpret_s16_s64(v_eobmax_xx32));
  const int64x1_t v_eobmax_xxx3 =
      vshr_n_s64(vreinterpret_s64_s16(v_eobmax_tmp), 16);
  const int16x4_t v_eobmax_final =
      vmax_s16(v_eobmax_tmp, vreinterpret_s16_s64(v_eobmax_xxx3));

  return (uint16_t)vget_lane_s16(v_eobmax_final, 0);
#endif  // __aarch64__
}

static VPX_FORCE_INLINE void load_fp_values(const int16_t *round_ptr,
                                            const int16_t *quant_ptr,
                                            const int16_t *dequant_ptr,
                                            int16x8_t *round, int16x8_t *quant,
                                            int16x8_t *dequant) {
  *round = vld1q_s16(round_ptr);
  *quant = vld1q_s16(quant_ptr);
  *dequant = vld1q_s16(dequant_ptr);
}

static VPX_FORCE_INLINE void update_fp_values(int16x8_t *v_round,
                                              int16x8_t *v_quant,
                                              int16x8_t *v_dequant) {
#ifdef __aarch64__
  *v_round = vdupq_laneq_s16(*v_round, 1);
  *v_quant = vdupq_laneq_s16(*v_quant, 1);
  *v_dequant = vdupq_laneq_s16(*v_dequant, 1);
#else
  *v_round = vdupq_lane_s16(vget_low_s16(*v_round), 1);
  *v_quant = vdupq_lane_s16(vget_low_s16(*v_quant), 1);
  *v_dequant = vdupq_lane_s16(vget_low_s16(*v_dequant), 1);
#endif
}

static VPX_FORCE_INLINE void quantize_fp_8(
    const int16x8_t *v_round, const int16x8_t *v_quant,
    const int16x8_t *v_dequant, const tran_low_t *coeff_ptr,
    const int16_t *iscan_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
    int16x8_t *v_eobmax) {
  const int16x8_t v_zero = vdupq_n_s16(0);
  const int16x8_t v_coeff = load_tran_low_to_s16q(coeff_ptr);
  const int16x8_t v_coeff_sign = vshrq_n_s16(v_coeff, 15);
  const int16x8_t v_abs = vabsq_s16(v_coeff);
  const int16x8_t v_tmp = vqaddq_s16(v_abs, *v_round);
  const int32x4_t v_tmp_lo =
      vmull_s16(vget_low_s16(v_tmp), vget_low_s16(*v_quant));
  const int32x4_t v_tmp_hi =
      vmull_s16(vget_high_s16(v_tmp), vget_high_s16(*v_quant));
  const int16x8_t v_tmp2 =
      vcombine_s16(vshrn_n_s32(v_tmp_lo, 16), vshrn_n_s32(v_tmp_hi, 16));
  const uint16x8_t v_nz_mask = vceqq_s16(v_tmp2, v_zero);
  const int16x8_t v_qcoeff_a = veorq_s16(v_tmp2, v_coeff_sign);
  const int16x8_t v_qcoeff = vsubq_s16(v_qcoeff_a, v_coeff_sign);
  calculate_dqcoeff_and_store(v_qcoeff, *v_dequant, dqcoeff_ptr);
  store_s16q_to_tran_low(qcoeff_ptr, v_qcoeff);

  *v_eobmax = get_max_lane_eob(iscan_ptr, *v_eobmax, v_nz_mask);
}

void vp9_quantize_fp_neon(const tran_low_t *coeff_ptr, intptr_t count,
                          const int16_t *round_ptr, const int16_t *quant_ptr,
                          tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                          const int16_t *dequant_ptr, uint16_t *eob_ptr,
                          const int16_t *scan, const int16_t *iscan) {
  // Quantization pass: All coefficients with index >= zero_flag are
  // skippable. Note: zero_flag can be zero.
  int i;
  int16x8_t v_eobmax = vdupq_n_s16(-1);
  int16x8_t v_round, v_quant, v_dequant;
  (void)scan;

  load_fp_values(round_ptr, quant_ptr, dequant_ptr, &v_round, &v_quant,
                 &v_dequant);
  // process dc and the first seven ac coeffs
  quantize_fp_8(&v_round, &v_quant, &v_dequant, coeff_ptr, iscan, qcoeff_ptr,
                dqcoeff_ptr, &v_eobmax);

  // now process the rest of the ac coeffs
  update_fp_values(&v_round, &v_quant, &v_dequant);
  for (i = 8; i < count; i += 8) {
    quantize_fp_8(&v_round, &v_quant, &v_dequant, coeff_ptr + i, iscan + i,
                  qcoeff_ptr + i, dqcoeff_ptr + i, &v_eobmax);
  }

  *eob_ptr = get_max_eob(v_eobmax);
}

static INLINE int32x4_t extract_sign_bit(int32x4_t a) {
  return vreinterpretq_s32_u32(vshrq_n_u32(vreinterpretq_u32_s32(a), 31));
}

static VPX_FORCE_INLINE void quantize_fp_32x32_8(
    const int16x8_t *v_round, const int16x8_t *v_quant,
    const int16x8_t *v_dequant, const int16x8_t *dequant_thresh,
    const tran_low_t *coeff_ptr, const int16_t *iscan_ptr,
    tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, int16x8_t *v_eobmax) {
  const int16x8_t v_coeff = load_tran_low_to_s16q(coeff_ptr);
  const int16x8_t v_coeff_sign = vshrq_n_s16(v_coeff, 15);
  const int16x8_t v_coeff_abs = vabsq_s16(v_coeff);
  const int16x8_t v_thr_mask =
      vreinterpretq_s16_u16(vcgeq_s16(v_coeff_abs, *dequant_thresh));
  const int16x8_t v_tmp_rnd =
      vandq_s16(vqaddq_s16(v_coeff_abs, *v_round), v_thr_mask);
  const int16x8_t v_abs_qcoeff = vqdmulhq_s16(v_tmp_rnd, *v_quant);
  const int16x8_t v_qcoeff =
      vsubq_s16(veorq_s16(v_abs_qcoeff, v_coeff_sign), v_coeff_sign);
  const uint16x8_t v_nz_mask = vceqq_s16(v_abs_qcoeff, vdupq_n_s16(0));

  int32x4_t dqcoeff_0, dqcoeff_1;
  dqcoeff_0 = vmull_s16(vget_low_s16(v_qcoeff), vget_low_s16(*v_dequant));
  dqcoeff_1 = vmull_s16(vget_high_s16(v_qcoeff), vget_high_s16(*v_dequant));
  // Add 1 if negative to round towards zero because the C uses division.
  dqcoeff_0 = vaddq_s32(dqcoeff_0, extract_sign_bit(dqcoeff_0));
  dqcoeff_1 = vaddq_s32(dqcoeff_1, extract_sign_bit(dqcoeff_1));

#if CONFIG_VP9_HIGHBITDEPTH
  vst1q_s32(dqcoeff_ptr, vshrq_n_s32(dqcoeff_0, 1));
  vst1q_s32(dqcoeff_ptr + 4, vshrq_n_s32(dqcoeff_1, 1));
#else
  store_s16q_to_tran_low(dqcoeff_ptr, vcombine_s16(vshrn_n_s32(dqcoeff_0, 1),
                                                   vshrn_n_s32(dqcoeff_1, 1)));
#endif

  store_s16q_to_tran_low(qcoeff_ptr, v_qcoeff);

  *v_eobmax = get_max_lane_eob(iscan_ptr, *v_eobmax, v_nz_mask);
}

void vp9_quantize_fp_32x32_neon(const tran_low_t *coeff_ptr, intptr_t count,
                                const int16_t *round_ptr,
                                const int16_t *quant_ptr,
                                tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                                const int16_t *dequant_ptr, uint16_t *eob_ptr,
                                const int16_t *scan, const int16_t *iscan) {
  int16x8_t eob_max = vdupq_n_s16(-1);
  // ROUND_POWER_OF_TWO(round_ptr[], 1)
  int16x8_t round = vrshrq_n_s16(vld1q_s16(round_ptr), 1);
  int16x8_t quant = vld1q_s16(quant_ptr);
  int16x8_t dequant = vld1q_s16(dequant_ptr);
  // dequant >> 2 is used similar to zbin as a threshold.
  int16x8_t dequant_thresh = vshrq_n_s16(vld1q_s16(dequant_ptr), 2);
  int i;

  (void)scan;
  (void)count;

  // Process dc and the first seven ac coeffs.
  quantize_fp_32x32_8(&round, &quant, &dequant, &dequant_thresh, coeff_ptr,
                      iscan, qcoeff_ptr, dqcoeff_ptr, &eob_max);

  update_fp_values(&round, &quant, &dequant);
  dequant_thresh = vdupq_lane_s16(vget_low_s16(dequant_thresh), 1);

  iscan += 8;
  coeff_ptr += 8;
  qcoeff_ptr += 8;
  dqcoeff_ptr += 8;

  // Process the rest of the ac coeffs.
  for (i = 8; i < 32 * 32; i += 8) {
    quantize_fp_32x32_8(&round, &quant, &dequant, &dequant_thresh, coeff_ptr,
                        iscan, qcoeff_ptr, dqcoeff_ptr, &eob_max);

    iscan += 8;
    coeff_ptr += 8;
    qcoeff_ptr += 8;
    dqcoeff_ptr += 8;
  }

  *eob_ptr = get_max_eob(eob_max);
}
