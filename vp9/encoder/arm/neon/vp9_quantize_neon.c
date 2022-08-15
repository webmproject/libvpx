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

void vp9_quantize_fp_32x32_neon(const tran_low_t *coeff_ptr, intptr_t count,
                                const int16_t *round_ptr,
                                const int16_t *quant_ptr,
                                tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                                const int16_t *dequant_ptr, uint16_t *eob_ptr,
                                const int16_t *scan, const int16_t *iscan) {
  const int16x8_t one = vdupq_n_s16(1);
  const int16x8_t neg_one = vdupq_n_s16(-1);

  // ROUND_POWER_OF_TWO(round_ptr[], 1)
  const int16x8_t round = vrshrq_n_s16(vld1q_s16(round_ptr), 1);
  const int16x8_t quant = vld1q_s16(quant_ptr);
  const int16x4_t dequant = vld1_s16(dequant_ptr);
  // dequant >> 2 is used similar to zbin as a threshold.
  const int16x8_t dequant_thresh = vshrq_n_s16(vld1q_s16(dequant_ptr), 2);

  // Process dc and the first seven ac coeffs.
  const uint16x8_t v_iscan =
      vreinterpretq_u16_s16(vaddq_s16(vld1q_s16(iscan), one));
  const int16x8_t coeff = load_tran_low_to_s16q(coeff_ptr);
  const int16x8_t coeff_sign = vshrq_n_s16(coeff, 15);
  const int16x8_t coeff_abs = vabsq_s16(coeff);
  const int16x8_t dequant_mask =
      vreinterpretq_s16_u16(vcgeq_s16(coeff_abs, dequant_thresh));

  int16x8_t qcoeff = vqaddq_s16(coeff_abs, round);
  int32x4_t dqcoeff_0, dqcoeff_1;
  uint16x8_t eob_max;
  (void)scan;
  (void)count;

  // coeff * quant_ptr[]) >> 15
  qcoeff = vqdmulhq_s16(qcoeff, quant);

  // Restore sign.
  qcoeff = veorq_s16(qcoeff, coeff_sign);
  qcoeff = vsubq_s16(qcoeff, coeff_sign);
  qcoeff = vandq_s16(qcoeff, dequant_mask);

  // qcoeff * dequant[] / 2
  dqcoeff_0 = vmull_s16(vget_low_s16(qcoeff), dequant);
  dqcoeff_1 = vmull_n_s16(vget_high_s16(qcoeff), dequant_ptr[1]);

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

  eob_max = vandq_u16(vtstq_s16(qcoeff, neg_one), v_iscan);

  store_s16q_to_tran_low(qcoeff_ptr, qcoeff);

  iscan += 8;
  coeff_ptr += 8;
  qcoeff_ptr += 8;
  dqcoeff_ptr += 8;

  {
    int i;
    const int16x8_t round = vrshrq_n_s16(vmovq_n_s16(round_ptr[1]), 1);
    const int16x8_t quant = vmovq_n_s16(quant_ptr[1]);
    const int16x8_t dequant_thresh =
        vshrq_n_s16(vmovq_n_s16(dequant_ptr[1]), 2);

    // Process the rest of the ac coeffs.
    for (i = 8; i < 32 * 32; i += 8) {
      const uint16x8_t v_iscan =
          vreinterpretq_u16_s16(vaddq_s16(vld1q_s16(iscan), one));
      const int16x8_t coeff = load_tran_low_to_s16q(coeff_ptr);
      const int16x8_t coeff_sign = vshrq_n_s16(coeff, 15);
      const int16x8_t coeff_abs = vabsq_s16(coeff);
      const int16x8_t dequant_mask =
          vreinterpretq_s16_u16(vcgeq_s16(coeff_abs, dequant_thresh));

      int16x8_t qcoeff = vqaddq_s16(coeff_abs, round);
      int32x4_t dqcoeff_0, dqcoeff_1;

      qcoeff = vqdmulhq_s16(qcoeff, quant);
      qcoeff = veorq_s16(qcoeff, coeff_sign);
      qcoeff = vsubq_s16(qcoeff, coeff_sign);
      qcoeff = vandq_s16(qcoeff, dequant_mask);

      dqcoeff_0 = vmull_n_s16(vget_low_s16(qcoeff), dequant_ptr[1]);
      dqcoeff_1 = vmull_n_s16(vget_high_s16(qcoeff), dequant_ptr[1]);

      dqcoeff_0 = vaddq_s32(dqcoeff_0, extract_sign_bit(dqcoeff_0));
      dqcoeff_1 = vaddq_s32(dqcoeff_1, extract_sign_bit(dqcoeff_1));

#if CONFIG_VP9_HIGHBITDEPTH
      vst1q_s32(dqcoeff_ptr, vshrq_n_s32(dqcoeff_0, 1));
      vst1q_s32(dqcoeff_ptr + 4, vshrq_n_s32(dqcoeff_1, 1));
#else
      store_s16q_to_tran_low(
          dqcoeff_ptr,
          vcombine_s16(vshrn_n_s32(dqcoeff_0, 1), vshrn_n_s32(dqcoeff_1, 1)));
#endif

      eob_max =
          vmaxq_u16(eob_max, vandq_u16(vtstq_s16(qcoeff, neg_one), v_iscan));

      store_s16q_to_tran_low(qcoeff_ptr, qcoeff);

      iscan += 8;
      coeff_ptr += 8;
      qcoeff_ptr += 8;
      dqcoeff_ptr += 8;
    }

#ifdef __aarch64__
    *eob_ptr = vmaxvq_u16(eob_max);
#else
    {
      const uint16x4_t eob_max_0 =
          vmax_u16(vget_low_u16(eob_max), vget_high_u16(eob_max));
      const uint16x4_t eob_max_1 = vpmax_u16(eob_max_0, eob_max_0);
      const uint16x4_t eob_max_2 = vpmax_u16(eob_max_1, eob_max_1);
      vst1_lane_u16(eob_ptr, eob_max_2, 0);
    }
#endif  // __aarch64__
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
static VPX_FORCE_INLINE uint16x4_t
highbd_quantize_fp_4(const tran_low_t *coeff_ptr, tran_low_t *qcoeff_ptr,
                     tran_low_t *dqcoeff_ptr, int32x4_t v_quant_s32,
                     int32x4_t v_dequant_s32, int32x4_t v_round_s32) {
  const int32x4_t v_coeff = vld1q_s32(coeff_ptr);
  const int32x4_t v_coeff_sign =
      vreinterpretq_s32_u32(vcltq_s32(v_coeff, vdupq_n_s32(0)));
  const int32x4_t v_abs_coeff = vabsq_s32(v_coeff);
  const int32x4_t v_tmp = vaddq_s32(v_abs_coeff, v_round_s32);
  //  const int abs_qcoeff = (int)((tmp * quant) >> 16);
  const int32x4_t v_abs_qcoeff = vqdmulhq_s32(v_tmp, v_quant_s32);
  //  qcoeff_ptr[rc] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
  const int32x4_t v_qcoeff =
      vsubq_s32(veorq_s32(v_abs_qcoeff, v_coeff_sign), v_coeff_sign);
  const int32x4_t v_abs_dqcoeff = vmulq_s32(v_abs_qcoeff, v_dequant_s32);
  //  dqcoeff_ptr[rc] = (tran_low_t)((abs_dqcoeff ^ coeff_sign) - coeff_sign);
  const int32x4_t v_dqcoeff =
      vsubq_s32(veorq_s32(v_abs_dqcoeff, v_coeff_sign), v_coeff_sign);

  vst1q_s32(qcoeff_ptr, v_qcoeff);
  vst1q_s32(dqcoeff_ptr, v_dqcoeff);

  // Packed nz_qcoeff_mask. Used to find eob.
  return vmovn_u32(vceqq_s32(v_abs_qcoeff, vdupq_n_s32(0)));
}

void vp9_highbd_quantize_fp_neon(const tran_low_t *coeff_ptr, intptr_t count,
                                 const int16_t *round_ptr,
                                 const int16_t *quant_ptr,
                                 tran_low_t *qcoeff_ptr,
                                 tran_low_t *dqcoeff_ptr,
                                 const int16_t *dequant_ptr, uint16_t *eob_ptr,
                                 const int16_t *scan, const int16_t *iscan) {
  const int16x4_t v_zero = vdup_n_s16(0);
  const int16x4_t v_quant = vld1_s16(quant_ptr);
  const int16x4_t v_dequant = vld1_s16(dequant_ptr);
  const int16x4_t v_round = vld1_s16(round_ptr);
  int32x4_t v_round_s32 = vaddl_s16(v_round, v_zero);
  int32x4_t v_quant_s32 = vshlq_n_s32(vaddl_s16(v_quant, v_zero), 15);
  int32x4_t v_dequant_s32 = vaddl_s16(v_dequant, v_zero);
  uint16x4_t v_mask_lo, v_mask_hi;
  int16x8_t v_eobmax = vdupq_n_s16(-1);

  (void)scan;

  // DC and first 3 AC
  v_mask_lo = highbd_quantize_fp_4(coeff_ptr, qcoeff_ptr, dqcoeff_ptr,
                                   v_quant_s32, v_dequant_s32, v_round_s32);

  // overwrite the DC constants with AC constants
  v_round_s32 = vdupq_lane_s32(vget_low_s32(v_round_s32), 1);
  v_quant_s32 = vdupq_lane_s32(vget_low_s32(v_quant_s32), 1);
  v_dequant_s32 = vdupq_lane_s32(vget_low_s32(v_dequant_s32), 1);

  // 4 more AC
  v_mask_hi =
      highbd_quantize_fp_4(coeff_ptr + 4, qcoeff_ptr + 4, dqcoeff_ptr + 4,
                           v_quant_s32, v_dequant_s32, v_round_s32);

  // Find the max lane eob for the first 8 coeffs.
  v_eobmax =
      get_max_lane_eob(iscan, v_eobmax, vcombine_u16(v_mask_lo, v_mask_hi));

  count -= 8;
  do {
    coeff_ptr += 8;
    qcoeff_ptr += 8;
    dqcoeff_ptr += 8;
    iscan += 8;
    v_mask_lo = highbd_quantize_fp_4(coeff_ptr, qcoeff_ptr, dqcoeff_ptr,
                                     v_quant_s32, v_dequant_s32, v_round_s32);
    v_mask_hi =
        highbd_quantize_fp_4(coeff_ptr + 4, qcoeff_ptr + 4, dqcoeff_ptr + 4,
                             v_quant_s32, v_dequant_s32, v_round_s32);
    // Find the max lane eob for 8 coeffs.
    v_eobmax =
        get_max_lane_eob(iscan, v_eobmax, vcombine_u16(v_mask_lo, v_mask_hi));
    count -= 8;
  } while (count);

  *eob_ptr = get_max_eob(v_eobmax);
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
