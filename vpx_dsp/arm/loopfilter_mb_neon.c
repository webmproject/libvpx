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

// Should we apply any filter at all: 11111111 yes, 00000000 no
static INLINE uint8x8_t filter_mask(
    const uint8x8_t limit, const uint8x8_t blimit, const uint8x8_t thresh,
    const uint8x8_t p3, const uint8x8_t p2, const uint8x8_t p1,
    const uint8x8_t p0, const uint8x8_t q0, const uint8x8_t q1,
    const uint8x8_t q2, const uint8x8_t q3, uint8x8_t *flat, uint8x8_t *hev) {
  uint8x8_t t0, t1;
  uint8x8_t max = vabd_u8(p1, p0);
  max = vmax_u8(max, vabd_u8(q1, q0));

  // Is there high edge variance internal edge: 11111111 yes, 00000000 no
  *hev = vcgt_u8(max, thresh);
  *flat = vmax_u8(max, vabd_u8(p2, p0));
  max = vmax_u8(max, vabd_u8(p3, p2));
  max = vmax_u8(max, vabd_u8(p2, p1));
  max = vmax_u8(max, vabd_u8(q2, q1));
  max = vmax_u8(max, vabd_u8(q3, q2));
  t0 = vabd_u8(p0, q0);
  t1 = vabd_u8(p1, q1);
  t0 = vqshl_n_u8(t0, 1);
  t1 = vshr_n_u8(t1, 1);
  t0 = vqadd_u8(t0, t1);
  max = vcle_u8(max, limit);
  t0 = vcle_u8(t0, blimit);
  max = vand_u8(max, t0);

  *flat = vmax_u8(*flat, vabd_u8(q2, q0));
  *flat = vmax_u8(*flat, vabd_u8(p3, p0));
  *flat = vmax_u8(*flat, vabd_u8(q3, q0));
  *flat = vcle_u8(*flat, vdup_n_u8(1));  // flat_mask4()

  return max;
}

static INLINE uint8x8_t flat_mask5(const uint8x8_t p4, const uint8x8_t p3,
                                   const uint8x8_t p2, const uint8x8_t p1,
                                   const uint8x8_t p0, const uint8x8_t q0,
                                   const uint8x8_t q1, const uint8x8_t q2,
                                   const uint8x8_t q3, const uint8x8_t q4) {
  uint8x8_t max = vabd_u8(p4, p0);
  max = vmax_u8(max, vabd_u8(p3, p0));
  max = vmax_u8(max, vabd_u8(p2, p0));
  max = vmax_u8(max, vabd_u8(p1, p0));
  max = vmax_u8(max, vabd_u8(q1, q0));
  max = vmax_u8(max, vabd_u8(q2, q0));
  max = vmax_u8(max, vabd_u8(q3, q0));
  max = vmax_u8(max, vabd_u8(q4, q0));
  max = vcle_u8(max, vdup_n_u8(1));

  return max;
}

static INLINE int8x8_t flip_sign(const uint8x8_t v) {
  const uint8x8_t sign_bit = vdup_n_u8(0x80);
  return vreinterpret_s8_u8(veor_u8(v, sign_bit));
}

static INLINE uint8x8_t flip_sign_back(const int8x8_t v) {
  const int8x8_t sign_bit = vdup_n_s8(0x80);
  return vreinterpret_u8_s8(veor_s8(v, sign_bit));
}

static INLINE uint8x8_t filter_tap7(const uint8x8_t flat, const uint8x8_t sub0,
                                    const uint8x8_t sub1, const uint8x8_t add0,
                                    const uint8x8_t add1, const uint8x8_t in,
                                    uint16x8_t *sum) {
  *sum = vsubw_u8(*sum, sub0);
  *sum = vsubw_u8(*sum, sub1);
  *sum = vaddw_u8(*sum, add0);
  *sum = vaddw_u8(*sum, add1);
  return vbsl_u8(flat, vrshrn_n_u16(*sum, 3), in);
}

static INLINE uint8x8_t filter_tap15(const uint8x8_t flat, const uint8x8_t sub0,
                                     const uint8x8_t sub1, const uint8x8_t add0,
                                     const uint8x8_t add1, const uint8x8_t in,
                                     uint16x8_t *sum) {
  *sum = vsubw_u8(*sum, sub0);
  *sum = vsubw_u8(*sum, sub1);
  *sum = vaddw_u8(*sum, add0);
  *sum = vaddw_u8(*sum, add1);
  return vbsl_u8(flat, vrshrn_n_u16(*sum, 4), in);
}

// 7-tap filter [1, 1, 1, 2, 1, 1, 1]
static INLINE void apply_7_tap_filter(const uint8x8_t flat, const uint8x8_t p3,
                                      const uint8x8_t p2, const uint8x8_t p1,
                                      const uint8x8_t p0, const uint8x8_t q0,
                                      const uint8x8_t q1, const uint8x8_t q2,
                                      const uint8x8_t q3, uint8x8_t *op2,
                                      uint8x8_t *op1, uint8x8_t *op0,
                                      uint8x8_t *oq0, uint8x8_t *oq1,
                                      uint8x8_t *oq2) {
  uint16x8_t sum;
  sum = vaddl_u8(p3, p3);   // 2*p3
  sum = vaddw_u8(sum, p3);  // 3*p3
  sum = vaddw_u8(sum, p2);  // 3*p3+p2
  sum = vaddw_u8(sum, p2);  // 3*p3+2*p2
  sum = vaddw_u8(sum, p1);  // 3*p3+2*p2+p1
  sum = vaddw_u8(sum, p0);  // 3*p3+2*p2+p1+p0
  sum = vaddw_u8(sum, q0);  // 3*p3+2*p2+p1+p0+q0
  *op2 = vbsl_u8(flat, vrshrn_n_u16(sum, 3), p2);
  *op1 = filter_tap7(flat, p3, p2, p1, q1, *op1, &sum);
  *op0 = filter_tap7(flat, p3, p1, p0, q2, *op0, &sum);
  *oq0 = filter_tap7(flat, p3, p0, q0, q3, *oq0, &sum);
  *oq1 = filter_tap7(flat, p2, q0, q1, q3, *oq1, &sum);
  *oq2 = filter_tap7(flat, p1, q1, q2, q3, q2, &sum);
}

// 15-tap filter [1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1]
static INLINE void apply_15_tap_filter(
    const uint8x8_t flat2, const uint8x8_t p7, const uint8x8_t p6,
    const uint8x8_t p5, const uint8x8_t p4, const uint8x8_t p3,
    const uint8x8_t p2, const uint8x8_t p1, const uint8x8_t p0,
    const uint8x8_t q0, const uint8x8_t q1, const uint8x8_t q2,
    const uint8x8_t q3, const uint8x8_t q4, const uint8x8_t q5,
    const uint8x8_t q6, const uint8x8_t q7, uint8x8_t *op6, uint8x8_t *op5,
    uint8x8_t *op4, uint8x8_t *op3, uint8x8_t *op2, uint8x8_t *op1,
    uint8x8_t *op0, uint8x8_t *oq0, uint8x8_t *oq1, uint8x8_t *oq2,
    uint8x8_t *oq3, uint8x8_t *oq4, uint8x8_t *oq5, uint8x8_t *oq6) {
  uint16x8_t sum;
  sum = vshll_n_u8(p7, 3);  // 8*p7
  sum = vsubw_u8(sum, p7);  // 7*p7
  sum = vaddw_u8(sum, p6);  // 7*p7+p6
  sum = vaddw_u8(sum, p6);  // 7*p7+2*p6
  sum = vaddw_u8(sum, p5);  // 7*p7+2*p6+p5
  sum = vaddw_u8(sum, p4);  // 7*p7+2*p6+p5+p4
  sum = vaddw_u8(sum, p3);  // 7*p7+2*p6+p5+p4+p3
  sum = vaddw_u8(sum, p2);  // 7*p7+2*p6+p5+p4+p3+p2
  sum = vaddw_u8(sum, p1);  // 7*p7+2*p6+p5+p4+p3+p2+p1
  sum = vaddw_u8(sum, p0);  // 7*p7+2*p6+p5+p4+p3+p2+p1+p0
  sum = vaddw_u8(sum, q0);  // 7*p7+2*p6+p5+p4+p3+p2+p1+p0+q0
  *op6 = vbsl_u8(flat2, vrshrn_n_u16(sum, 4), p6);
  *op5 = filter_tap15(flat2, p7, p6, p5, q1, p5, &sum);
  *op4 = filter_tap15(flat2, p7, p5, p4, q2, p4, &sum);
  *op3 = filter_tap15(flat2, p7, p4, p3, q3, p3, &sum);
  *op2 = filter_tap15(flat2, p7, p3, p2, q4, *op2, &sum);
  *op1 = filter_tap15(flat2, p7, p2, p1, q5, *op1, &sum);
  *op0 = filter_tap15(flat2, p7, p1, p0, q6, *op0, &sum);
  *oq0 = filter_tap15(flat2, p7, p0, q0, q7, *oq0, &sum);
  *oq1 = filter_tap15(flat2, p6, q0, q1, q7, *oq1, &sum);
  *oq2 = filter_tap15(flat2, p5, q1, q2, q7, *oq2, &sum);
  *oq3 = filter_tap15(flat2, p4, q2, q3, q7, q3, &sum);
  *oq4 = filter_tap15(flat2, p3, q3, q4, q7, q4, &sum);
  *oq5 = filter_tap15(flat2, p2, q4, q5, q7, q5, &sum);
  *oq6 = filter_tap15(flat2, p1, q5, q6, q7, q6, &sum);
}

static INLINE void filter16(
    const uint8x8_t mask, const uint8x8_t flat, const uint64_t flat_u64,
    const uint8x8_t flat2, const uint64_t flat2_u64, const uint8x8_t hev,
    const uint8x8_t p7, const uint8x8_t p6, const uint8x8_t p5,
    const uint8x8_t p4, const uint8x8_t p3, const uint8x8_t p2,
    const uint8x8_t p1, const uint8x8_t p0, const uint8x8_t q0,
    const uint8x8_t q1, const uint8x8_t q2, const uint8x8_t q3,
    const uint8x8_t q4, const uint8x8_t q5, const uint8x8_t q6,
    const uint8x8_t q7, uint8x8_t *op6, uint8x8_t *op5, uint8x8_t *op4,
    uint8x8_t *op3, uint8x8_t *op2, uint8x8_t *op1, uint8x8_t *op0,
    uint8x8_t *oq0, uint8x8_t *oq1, uint8x8_t *oq2, uint8x8_t *oq3,
    uint8x8_t *oq4, uint8x8_t *oq5, uint8x8_t *oq6) {
  // add outer taps if we have high edge variance
  if (flat_u64 != (uint64_t)-1) {
    int8x8_t filter, filter1, filter2, t;
    int8x8_t ps1 = flip_sign(p1);
    int8x8_t ps0 = flip_sign(p0);
    int8x8_t qs0 = flip_sign(q0);
    int8x8_t qs1 = flip_sign(q1);

    filter = vqsub_s8(ps1, qs1);
    filter = vand_s8(filter, vreinterpret_s8_u8(hev));
    t = vqsub_s8(qs0, ps0);

    // inner taps
    filter = vqadd_s8(filter, t);
    filter = vqadd_s8(filter, t);
    filter = vqadd_s8(filter, t);
    filter = vand_s8(filter, vreinterpret_s8_u8(mask));

    // save bottom 3 bits so that we round one side +4 and the other +3
    // if it equals 4 we'll set to adjust by -1 to account for the fact
    // we'd round 3 the other way
    filter1 = vshr_n_s8(vqadd_s8(filter, vdup_n_s8(4)), 3);
    filter2 = vshr_n_s8(vqadd_s8(filter, vdup_n_s8(3)), 3);

    qs0 = vqsub_s8(qs0, filter1);
    ps0 = vqadd_s8(ps0, filter2);
    *oq0 = flip_sign_back(qs0);
    *op0 = flip_sign_back(ps0);

    // outer tap adjustments
    filter = vrshr_n_s8(filter1, 1);
    filter = vbic_s8(filter, vreinterpret_s8_u8(hev));

    qs1 = vqsub_s8(qs1, filter);
    ps1 = vqadd_s8(ps1, filter);
    *oq1 = flip_sign_back(qs1);
    *op1 = flip_sign_back(ps1);
  }

  if (flat_u64) {
    *op2 = p2;
    *oq2 = q2;
    if (flat2_u64 != (uint64_t)-1) {
      apply_7_tap_filter(flat, p3, p2, p1, p0, q0, q1, q2, q3, op2, op1, op0,
                         oq0, oq1, oq2);
    }
    if (flat2_u64) {
      apply_15_tap_filter(flat2, p7, p6, p5, p4, p3, p2, p1, p0, q0, q1, q2, q3,
                          q4, q5, q6, q7, op6, op5, op4, op3, op2, op1, op0,
                          oq0, oq1, oq2, oq3, oq4, oq5, oq6);
    }
  }
}

static void mb_lpf_horizontal_edge_w(uint8_t *s, int p, const uint8_t *blimit,
                                     const uint8_t *limit,
                                     const uint8_t *thresh, int count) {
  const uint8x8_t blimit_u8x8 = vld1_dup_u8(blimit);
  const uint8x8_t limit_u8x8 = vld1_dup_u8(limit);
  const uint8x8_t thresh_u8x8 = vld1_dup_u8(thresh);

  do {
    const uint8x8_t p7 = vld1_u8(s - 8 * p);
    const uint8x8_t p6 = vld1_u8(s - 7 * p);
    const uint8x8_t p5 = vld1_u8(s - 6 * p);
    const uint8x8_t p4 = vld1_u8(s - 5 * p);
    const uint8x8_t p3 = vld1_u8(s - 4 * p);
    const uint8x8_t p2 = vld1_u8(s - 3 * p);
    const uint8x8_t p1 = vld1_u8(s - 2 * p);
    const uint8x8_t p0 = vld1_u8(s - 1 * p);
    const uint8x8_t q0 = vld1_u8(s + 0 * p);
    const uint8x8_t q1 = vld1_u8(s + 1 * p);
    const uint8x8_t q2 = vld1_u8(s + 2 * p);
    const uint8x8_t q3 = vld1_u8(s + 3 * p);
    const uint8x8_t q4 = vld1_u8(s + 4 * p);
    const uint8x8_t q5 = vld1_u8(s + 5 * p);
    const uint8x8_t q6 = vld1_u8(s + 6 * p);
    const uint8x8_t q7 = vld1_u8(s + 7 * p);
    uint8x8_t op6, op5, op4, op3, op2, op1, op0, oq0, oq1, oq2, oq3, oq4, oq5,
        oq6, flat, hev;
    const uint8x8_t mask = filter_mask(limit_u8x8, blimit_u8x8, thresh_u8x8, p3,
                                       p2, p1, p0, q0, q1, q2, q3, &flat, &hev);
    uint8x8_t flat2 = flat_mask5(p7, p6, p5, p4, p0, q0, q4, q5, q6, q7);
    uint64_t flat_u64, flat2_u64;

    flat = vand_u8(flat, mask);
    flat2 = vand_u8(flat2, flat);
    flat_u64 = vget_lane_u64(vreinterpret_u64_u8(flat), 0);
    flat2_u64 = vget_lane_u64(vreinterpret_u64_u8(flat2), 0);

    filter16(mask, flat, flat_u64, flat2, flat2_u64, hev, p7, p6, p5, p4, p3,
             p2, p1, p0, q0, q1, q2, q3, q4, q5, q6, q7, &op6, &op5, &op4, &op3,
             &op2, &op1, &op0, &oq0, &oq1, &oq2, &oq3, &oq4, &oq5, &oq6);

    if (flat_u64) {
      if (flat2_u64) {
        vst1_u8(s - 7 * p, op6);
        vst1_u8(s - 6 * p, op5);
        vst1_u8(s - 5 * p, op4);
        vst1_u8(s - 4 * p, op3);
        vst1_u8(s + 3 * p, oq3);
        vst1_u8(s + 4 * p, oq4);
        vst1_u8(s + 5 * p, oq5);
        vst1_u8(s + 6 * p, oq6);
      }
      vst1_u8(s - 3 * p, op2);
      vst1_u8(s + 2 * p, oq2);
    }
    vst1_u8(s - 2 * p, op1);
    vst1_u8(s - 1 * p, op0);
    vst1_u8(s + 0 * p, oq0);
    vst1_u8(s + 1 * p, oq1);
    s += 8;
  } while (--count);
}

void vpx_lpf_horizontal_edge_8_neon(uint8_t *s, int p, const uint8_t *blimit,
                                    const uint8_t *limit,
                                    const uint8_t *thresh) {
  mb_lpf_horizontal_edge_w(s, p, blimit, limit, thresh, 1);
}

void vpx_lpf_horizontal_edge_16_neon(uint8_t *s, int p, const uint8_t *blimit,
                                     const uint8_t *limit,
                                     const uint8_t *thresh) {
  mb_lpf_horizontal_edge_w(s, p, blimit, limit, thresh, 2);
}

static void mb_lpf_vertical_edge_w(uint8_t *s, int p, const uint8_t *blimit,
                                   const uint8_t *limit, const uint8_t *thresh,
                                   int count) {
  const uint8x8_t blimit_u8x8 = vld1_dup_u8(blimit);
  const uint8x8_t limit_u8x8 = vld1_dup_u8(limit);
  const uint8x8_t thresh_u8x8 = vld1_dup_u8(thresh);
  uint8_t *d;

  s -= 8;
  d = s;
  do {
    uint8x16_t t0, t1, t2, t3, t4, t5, t6, t7;
    uint8x8_t p7, p6, p5, p4, p3, p2, p1, p0, q0, q1, q2, q3, q4, q5, q6, q7,
        op6, op5, op4, op3, op2, op1, op0, oq0, oq1, oq2, oq3, oq4, oq5, oq6,
        flat, hev, mask, flat2;
    uint64_t flat_u64, flat2_u64;

    t0 = vld1q_u8(s);
    s += p;
    t1 = vld1q_u8(s);
    s += p;
    t2 = vld1q_u8(s);
    s += p;
    t3 = vld1q_u8(s);
    s += p;
    t4 = vld1q_u8(s);
    s += p;
    t5 = vld1q_u8(s);
    s += p;
    t6 = vld1q_u8(s);
    s += p;
    t7 = vld1q_u8(s);
    s += p;

    transpose_u8_16x8(t0, t1, t2, t3, t4, t5, t6, t7, &p7, &p6, &p5, &p4, &p3,
                      &p2, &p1, &p0, &q0, &q1, &q2, &q3, &q4, &q5, &q6, &q7);

    mask = filter_mask(limit_u8x8, blimit_u8x8, thresh_u8x8, p3, p2, p1, p0, q0,
                       q1, q2, q3, &flat, &hev);
    flat2 = flat_mask5(p7, p6, p5, p4, p0, q0, q4, q5, q6, q7);
    flat = vand_u8(flat, mask);
    flat2 = vand_u8(flat2, flat);
    flat_u64 = vget_lane_u64(vreinterpret_u64_u8(flat), 0);
    flat2_u64 = vget_lane_u64(vreinterpret_u64_u8(flat2), 0);

    filter16(mask, flat, flat_u64, flat2, flat2_u64, hev, p7, p6, p5, p4, p3,
             p2, p1, p0, q0, q1, q2, q3, q4, q5, q6, q7, &op6, &op5, &op4, &op3,
             &op2, &op1, &op0, &oq0, &oq1, &oq2, &oq3, &oq4, &oq5, &oq6);

    if (flat_u64) {
      if (flat2_u64) {
        uint8x16_t o0, o1, o2, o3, o4, o5, o6, o7;
        transpose_u8_8x16(p7, op6, op5, op4, op3, op2, op1, op0, oq0, oq1, oq2,
                          oq3, oq4, oq5, oq6, q7, &o0, &o1, &o2, &o3, &o4, &o5,
                          &o6, &o7);

        vst1q_u8(d, o0);
        d += p;
        vst1q_u8(d, o1);
        d += p;
        vst1q_u8(d, o2);
        d += p;
        vst1q_u8(d, o3);
        d += p;
        vst1q_u8(d, o4);
        d += p;
        vst1q_u8(d, o5);
        d += p;
        vst1q_u8(d, o6);
        d += p;
        vst1q_u8(d, o7);
        d += p;
      } else {
        uint8x8x3_t o0, o1;
        d += 8;
        o0.val[0] = op2;
        o0.val[1] = op1;
        o0.val[2] = op0;
        o1.val[0] = oq0;
        o1.val[1] = oq1;
        o1.val[2] = oq2;
        vst3_lane_u8(d - 3, o0, 0);
        vst3_lane_u8(d + 0, o1, 0);
        d += p;
        vst3_lane_u8(d - 3, o0, 1);
        vst3_lane_u8(d + 0, o1, 1);
        d += p;
        vst3_lane_u8(d - 3, o0, 2);
        vst3_lane_u8(d + 0, o1, 2);
        d += p;
        vst3_lane_u8(d - 3, o0, 3);
        vst3_lane_u8(d + 0, o1, 3);
        d += p;
        vst3_lane_u8(d - 3, o0, 4);
        vst3_lane_u8(d + 0, o1, 4);
        d += p;
        vst3_lane_u8(d - 3, o0, 5);
        vst3_lane_u8(d + 0, o1, 5);
        d += p;
        vst3_lane_u8(d - 3, o0, 6);
        vst3_lane_u8(d + 0, o1, 6);
        d += p;
        vst3_lane_u8(d - 3, o0, 7);
        vst3_lane_u8(d + 0, o1, 7);
        d += p - 8;
      }
    } else {
      uint8x8x4_t o;
      d += 6;
      o.val[0] = op1;
      o.val[1] = op0;
      o.val[2] = oq0;
      o.val[3] = oq1;
      vst4_lane_u8(d, o, 0);
      d += p;
      vst4_lane_u8(d, o, 1);
      d += p;
      vst4_lane_u8(d, o, 2);
      d += p;
      vst4_lane_u8(d, o, 3);
      d += p;
      vst4_lane_u8(d, o, 4);
      d += p;
      vst4_lane_u8(d, o, 5);
      d += p;
      vst4_lane_u8(d, o, 6);
      d += p;
      vst4_lane_u8(d, o, 7);
      d += p - 6;
    }
  } while (--count);
}

void vpx_lpf_vertical_16_neon(uint8_t *s, int p, const uint8_t *blimit,
                              const uint8_t *limit, const uint8_t *thresh) {
  mb_lpf_vertical_edge_w(s, p, blimit, limit, thresh, 1);
}

void vpx_lpf_vertical_16_dual_neon(uint8_t *s, int p, const uint8_t *blimit,
                                   const uint8_t *limit,
                                   const uint8_t *thresh) {
  mb_lpf_vertical_edge_w(s, p, blimit, limit, thresh, 2);
}
