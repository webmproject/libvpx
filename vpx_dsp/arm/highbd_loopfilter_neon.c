/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <arm_neon.h>
#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/arm/transpose_neon.h"

static INLINE void load_thresh(const uint8_t *blimit, const uint8_t *limit,
                               const uint8_t *thresh, uint16x8_t *blimit_vec,
                               uint16x8_t *limit_vec, uint16x8_t *thresh_vec,
                               const int bd) {
  const int16x8_t shift = vdupq_n_s16(bd - 8);
  *blimit_vec = vmovl_u8(vld1_dup_u8(blimit));
  *limit_vec = vmovl_u8(vld1_dup_u8(limit));
  *thresh_vec = vmovl_u8(vld1_dup_u8(thresh));
  *blimit_vec = vshlq_u16(*blimit_vec, shift);
  *limit_vec = vshlq_u16(*limit_vec, shift);
  *thresh_vec = vshlq_u16(*thresh_vec, shift);
}

static INLINE uint16x8_t
filter_hev_mask4(const uint16x8_t limit, const uint16x8_t blimit,
                 const uint16x8_t thresh, const uint16x8_t p3,
                 const uint16x8_t p2, const uint16x8_t p1, const uint16x8_t p0,
                 const uint16x8_t q0, const uint16x8_t q1, const uint16x8_t q2,
                 const uint16x8_t q3, uint16x8_t *hev, uint16x8_t *mask) {
  uint16x8_t max, t0, t1;

  max = vabdq_u16(p1, p0);
  max = vmaxq_u16(max, vabdq_u16(q1, q0));
  *hev = vcgtq_u16(max, thresh);
  *mask = vmaxq_u16(max, vabdq_u16(p3, p2));
  *mask = vmaxq_u16(*mask, vabdq_u16(p2, p1));
  *mask = vmaxq_u16(*mask, vabdq_u16(q2, q1));
  *mask = vmaxq_u16(*mask, vabdq_u16(q3, q2));
  t0 = vabdq_u16(p0, q0);
  t1 = vabdq_u16(p1, q1);
  t0 = vaddq_u16(t0, t0);
  t1 = vshrq_n_u16(t1, 1);
  t0 = vaddq_u16(t0, t1);
  *mask = vcleq_u16(*mask, limit);
  t0 = vcleq_u16(t0, blimit);
  *mask = vandq_u16(*mask, t0);

  return max;
}

static INLINE int16x8_t flip_sign(const uint16x8_t v, const int bd) {
  const uint16x8_t offset = vdupq_n_u16(0x80 << (bd - 8));
  return vreinterpretq_s16_u16(vsubq_u16(v, offset));
}

static INLINE uint16x8_t flip_sign_back(const int16x8_t v, const int bd) {
  const int16x8_t offset = vdupq_n_s16(0x80 << (bd - 8));
  return vreinterpretq_u16_s16(vaddq_s16(v, offset));
}

static INLINE void filter4(const uint16x8_t mask, const uint16x8_t hev,
                           const uint16x8_t p1, const uint16x8_t p0,
                           const uint16x8_t q0, const uint16x8_t q1,
                           uint16x8_t *op1, uint16x8_t *op0, uint16x8_t *oq0,
                           uint16x8_t *oq1, const int bd) {
  const int16x8_t max = vdupq_n_s16((1 << (bd - 1)) - 1);
  const int16x8_t min = vdupq_n_s16((int16_t)(((uint32_t)-1) << (bd - 1)));
  int16x8_t filter, filter1, filter2, t;
  int16x8_t ps1 = flip_sign(p1, bd);
  int16x8_t ps0 = flip_sign(p0, bd);
  int16x8_t qs0 = flip_sign(q0, bd);
  int16x8_t qs1 = flip_sign(q1, bd);

  /* add outer taps if we have high edge variance */
  filter = vsubq_s16(ps1, qs1);
  filter = vmaxq_s16(filter, min);
  filter = vminq_s16(filter, max);
  filter = vandq_s16(filter, vreinterpretq_s16_u16(hev));
  t = vsubq_s16(qs0, ps0);

  /* inner taps */
  filter = vaddq_s16(filter, t);
  filter = vaddq_s16(filter, t);
  filter = vaddq_s16(filter, t);
  filter = vmaxq_s16(filter, min);
  filter = vminq_s16(filter, max);
  filter = vandq_s16(filter, vreinterpretq_s16_u16(mask));

  /* save bottom 3 bits so that we round one side +4 and the other +3 */
  /* if it equals 4 we'll set it to adjust by -1 to account for the fact */
  /* we'd round it by 3 the other way */
  t = vaddq_s16(filter, vdupq_n_s16(4));
  t = vminq_s16(t, max);
  filter1 = vshrq_n_s16(t, 3);
  t = vaddq_s16(filter, vdupq_n_s16(3));
  t = vminq_s16(t, max);
  filter2 = vshrq_n_s16(t, 3);

  qs0 = vsubq_s16(qs0, filter1);
  qs0 = vmaxq_s16(qs0, min);
  qs0 = vminq_s16(qs0, max);
  ps0 = vaddq_s16(ps0, filter2);
  ps0 = vmaxq_s16(ps0, min);
  ps0 = vminq_s16(ps0, max);
  *oq0 = flip_sign_back(qs0, bd);
  *op0 = flip_sign_back(ps0, bd);

  /* outer tap adjustments */
  filter = vrshrq_n_s16(filter1, 1);
  filter = vbicq_s16(filter, vreinterpretq_s16_u16(hev));

  qs1 = vsubq_s16(qs1, filter);
  qs1 = vmaxq_s16(qs1, min);
  qs1 = vminq_s16(qs1, max);
  ps1 = vaddq_s16(ps1, filter);
  ps1 = vmaxq_s16(ps1, min);
  ps1 = vminq_s16(ps1, max);
  *oq1 = flip_sign_back(qs1, bd);
  *op1 = flip_sign_back(ps1, bd);
}

static INLINE void load_8x8(const uint16_t *s, const int p, uint16x8_t *p3,
                            uint16x8_t *p2, uint16x8_t *p1, uint16x8_t *p0,
                            uint16x8_t *q0, uint16x8_t *q1, uint16x8_t *q2,
                            uint16x8_t *q3) {
  *p3 = vld1q_u16(s);
  s += p;
  *p2 = vld1q_u16(s);
  s += p;
  *p1 = vld1q_u16(s);
  s += p;
  *p0 = vld1q_u16(s);
  s += p;
  *q0 = vld1q_u16(s);
  s += p;
  *q1 = vld1q_u16(s);
  s += p;
  *q2 = vld1q_u16(s);
  s += p;
  *q3 = vld1q_u16(s);
}

static INLINE void store_8x4(uint16_t *s, const int p, const uint16x8_t s0,
                             const uint16x8_t s1, const uint16x8_t s2,
                             const uint16x8_t s3) {
  vst1q_u16(s, s0);
  s += p;
  vst1q_u16(s, s1);
  s += p;
  vst1q_u16(s, s2);
  s += p;
  vst1q_u16(s, s3);
}

static INLINE void store_4x8(uint16_t *s, const int p, const uint16x8_t p1,
                             const uint16x8_t p0, const uint16x8_t q0,
                             const uint16x8_t q1) {
  uint16x8x4_t o;

  o.val[0] = p1;
  o.val[1] = p0;
  o.val[2] = q0;
  o.val[3] = q1;
  vst4q_lane_u16(s, o, 0);
  s += p;
  vst4q_lane_u16(s, o, 1);
  s += p;
  vst4q_lane_u16(s, o, 2);
  s += p;
  vst4q_lane_u16(s, o, 3);
  s += p;
  vst4q_lane_u16(s, o, 4);
  s += p;
  vst4q_lane_u16(s, o, 5);
  s += p;
  vst4q_lane_u16(s, o, 6);
  s += p;
  vst4q_lane_u16(s, o, 7);
}

void vpx_highbd_lpf_horizontal_4_neon(uint16_t *s, int p, const uint8_t *blimit,
                                      const uint8_t *limit,
                                      const uint8_t *thresh, int bd) {
  uint16x8_t blimit_vec, limit_vec, thresh_vec, p3, p2, p1, p0, q0, q1, q2, q3,
      mask, hev;

  load_thresh(blimit, limit, thresh, &blimit_vec, &limit_vec, &thresh_vec, bd);
  load_8x8(s - 4 * p, p, &p3, &p2, &p1, &p0, &q0, &q1, &q2, &q3);
  filter_hev_mask4(limit_vec, blimit_vec, thresh_vec, p3, p2, p1, p0, q0, q1,
                   q2, q3, &hev, &mask);
  filter4(mask, hev, p1, p0, q0, q1, &p1, &p0, &q0, &q1, bd);
  store_8x4(s - 2 * p, p, p1, p0, q0, q1);
}

void vpx_highbd_lpf_horizontal_4_dual_neon(
    uint16_t *s, int p, const uint8_t *blimit0, const uint8_t *limit0,
    const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1,
    const uint8_t *thresh1, int bd) {
  vpx_highbd_lpf_horizontal_4_neon(s, p, blimit0, limit0, thresh0, bd);
  vpx_highbd_lpf_horizontal_4_neon(s + 8, p, blimit1, limit1, thresh1, bd);
}

void vpx_highbd_lpf_vertical_4_neon(uint16_t *s, int p, const uint8_t *blimit,
                                    const uint8_t *limit, const uint8_t *thresh,
                                    int bd) {
  uint16x8_t blimit_vec, limit_vec, thresh_vec, p3, p2, p1, p0, q0, q1, q2, q3,
      mask, hev;

  load_8x8(s - 4, p, &p3, &p2, &p1, &p0, &q0, &q1, &q2, &q3);
  transpose_s16_8x8((int16x8_t *)&p3, (int16x8_t *)&p2, (int16x8_t *)&p1,
                    (int16x8_t *)&p0, (int16x8_t *)&q0, (int16x8_t *)&q1,
                    (int16x8_t *)&q2, (int16x8_t *)&q3);
  load_thresh(blimit, limit, thresh, &blimit_vec, &limit_vec, &thresh_vec, bd);
  filter_hev_mask4(limit_vec, blimit_vec, thresh_vec, p3, p2, p1, p0, q0, q1,
                   q2, q3, &hev, &mask);
  filter4(mask, hev, p1, p0, q0, q1, &p1, &p0, &q0, &q1, bd);
  store_4x8(s - 2, p, p1, p0, q0, q1);
}

void vpx_highbd_lpf_vertical_4_dual_neon(
    uint16_t *s, int p, const uint8_t *blimit0, const uint8_t *limit0,
    const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1,
    const uint8_t *thresh1, int bd) {
  vpx_highbd_lpf_vertical_4_neon(s, p, blimit0, limit0, thresh0, bd);
  vpx_highbd_lpf_vertical_4_neon(s + 8 * p, p, blimit1, limit1, thresh1, bd);
}
