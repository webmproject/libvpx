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
#include "./vpx_config.h"

static INLINE void vp8_loop_filter_simple_vertical_edge_neon(
        unsigned char *s,
        int p,
        const unsigned char *blimit) {
    int pitch;
    unsigned char *src1, *src2;
    uint8x16_t qblimit, q0u8;
    uint8x16_t q3u8, q4u8, q5u8, q6u8, q7u8, q11u8, q12u8, q14u8, q15u8;
    int16x8_t q2s16, q13s16, q11s16;
    int8x8_t d28s8, d29s8;
    int8x16_t q2s8, q3s8, q10s8, q11s8, q14s8;
    uint8x8x4_t d0u8x4;  // d6, d7, d8, d9
    uint8x8x4_t d1u8x4;  // d10, d11, d12, d13
    uint8x8x2_t d2u8x2;  // d12, d13
    uint8x8x2_t d3u8x2;  // d14, d15

    pitch = p << 1;
    qblimit = vdupq_n_u8(*blimit);

    src1 = s - 2;
    d0u8x4 = vld4_lane_u8(src1, d0u8x4, 0);
    src1 += pitch;
    d0u8x4 = vld4_lane_u8(src1, d0u8x4, 2);
    src1 += pitch;
    d0u8x4 = vld4_lane_u8(src1, d0u8x4, 4);
    src1 += pitch;
    d0u8x4 = vld4_lane_u8(src1, d0u8x4, 6);
    src1 += pitch;
    d1u8x4 = vld4_lane_u8(src1, d1u8x4, 0);
    src1 += pitch;
    d1u8x4 = vld4_lane_u8(src1, d1u8x4, 2);
    src1 += pitch;
    d1u8x4 = vld4_lane_u8(src1, d1u8x4, 4);
    src1 += pitch;
    d1u8x4 = vld4_lane_u8(src1, d1u8x4, 6);

    src2 = s - 2 + p;
    d0u8x4 = vld4_lane_u8(src2, d0u8x4, 1);
    src2 += pitch;
    d0u8x4 = vld4_lane_u8(src2, d0u8x4, 3);
    src2 += pitch;
    d0u8x4 = vld4_lane_u8(src2, d0u8x4, 5);
    src2 += pitch;
    d0u8x4 = vld4_lane_u8(src2, d0u8x4, 7);
    src2 += pitch;
    d1u8x4 = vld4_lane_u8(src2, d1u8x4, 1);
    src2 += pitch;
    d1u8x4 = vld4_lane_u8(src2, d1u8x4, 3);
    src2 += pitch;
    d1u8x4 = vld4_lane_u8(src2, d1u8x4, 5);
    src2 += pitch;
    d1u8x4 = vld4_lane_u8(src2, d1u8x4, 7);

    // vswp    d7, d10
    // vswp    d12, d9
    q3u8 = vcombine_u8(d0u8x4.val[0], d1u8x4.val[0]);  // d6 d10
    q4u8 = vcombine_u8(d0u8x4.val[2], d1u8x4.val[2]);  // d8 d12
    q5u8 = vcombine_u8(d0u8x4.val[1], d1u8x4.val[1]);  // d7 d11
    q6u8 = vcombine_u8(d0u8x4.val[3], d1u8x4.val[3]);  // d9 d13

    q15u8 = vabdq_u8(q5u8, q4u8);
    q14u8 = vabdq_u8(q3u8, q6u8);

    q15u8 = vqaddq_u8(q15u8, q15u8);
    q14u8 = vshrq_n_u8(q14u8, 1);
    q0u8 = vdupq_n_u8(0x80);
    q11s16 = vdupq_n_s16(3);
    q15u8 = vqaddq_u8(q15u8, q14u8);

    q3u8 = veorq_u8(q3u8, q0u8);
    q4u8 = veorq_u8(q4u8, q0u8);
    q5u8 = veorq_u8(q5u8, q0u8);
    q6u8 = veorq_u8(q6u8, q0u8);

    q15u8 = vcgeq_u8(qblimit, q15u8);

    q2s16 = vsubl_s8(vget_low_s8(vreinterpretq_s8_u8(q4u8)),
                     vget_low_s8(vreinterpretq_s8_u8(q5u8)));
    q13s16 = vsubl_s8(vget_high_s8(vreinterpretq_s8_u8(q4u8)),
                      vget_high_s8(vreinterpretq_s8_u8(q5u8)));

    q14s8 = vqsubq_s8(vreinterpretq_s8_u8(q3u8),
                      vreinterpretq_s8_u8(q6u8));

    q2s16 = vmulq_s16(q2s16, q11s16);
    q13s16 = vmulq_s16(q13s16, q11s16);

    q11u8 = vdupq_n_u8(3);
    q12u8 = vdupq_n_u8(4);

    q2s16 = vaddw_s8(q2s16, vget_low_s8(q14s8));
    q13s16 = vaddw_s8(q13s16, vget_high_s8(q14s8));

    d28s8 = vqmovn_s16(q2s16);
    d29s8 = vqmovn_s16(q13s16);
    q14s8 = vcombine_s8(d28s8, d29s8);

    q14s8 = vandq_s8(q14s8, vreinterpretq_s8_u8(q15u8));

    q2s8 = vqaddq_s8(q14s8, vreinterpretq_s8_u8(q11u8));
    q3s8 = vqaddq_s8(q14s8, vreinterpretq_s8_u8(q12u8));
    q2s8 = vshrq_n_s8(q2s8, 3);
    q14s8 = vshrq_n_s8(q3s8, 3);

    q11s8 = vqaddq_s8(vreinterpretq_s8_u8(q5u8), q2s8);
    q10s8 = vqsubq_s8(vreinterpretq_s8_u8(q4u8), q14s8);

    q6u8 = veorq_u8(vreinterpretq_u8_s8(q11s8), q0u8);
    q7u8 = veorq_u8(vreinterpretq_u8_s8(q10s8), q0u8);

    d2u8x2.val[0] = vget_low_u8(q6u8);   // d12
    d2u8x2.val[1] = vget_low_u8(q7u8);   // d14
    d3u8x2.val[0] = vget_high_u8(q6u8);  // d13
    d3u8x2.val[1] = vget_high_u8(q7u8);  // d15

    src1 = s - 1;
    vst2_lane_u8(src1, d2u8x2, 0);
    src1 += pitch;
    vst2_lane_u8(src1, d2u8x2, 2);
    src1 += pitch;
    vst2_lane_u8(src1, d2u8x2, 4);
    src1 += pitch;
    vst2_lane_u8(src1, d2u8x2, 6);
    src1 += pitch;
    vst2_lane_u8(src1, d3u8x2, 0);
    src1 += pitch;
    vst2_lane_u8(src1, d3u8x2, 2);
    src1 += pitch;
    vst2_lane_u8(src1, d3u8x2, 4);
    src1 += pitch;
    vst2_lane_u8(src1, d3u8x2, 6);

    src2 = s - 1 + p;
    vst2_lane_u8(src2, d2u8x2, 1);
    src2 += pitch;
    vst2_lane_u8(src2, d2u8x2, 3);
    src2 += pitch;
    vst2_lane_u8(src2, d2u8x2, 5);
    src2 += pitch;
    vst2_lane_u8(src2, d2u8x2, 7);
    src2 += pitch;
    vst2_lane_u8(src2, d3u8x2, 1);
    src2 += pitch;
    vst2_lane_u8(src2, d3u8x2, 3);
    src2 += pitch;
    vst2_lane_u8(src2, d3u8x2, 5);
    src2 += pitch;
    vst2_lane_u8(src2, d3u8x2, 7);
    return;
}

void vp8_loop_filter_bvs_neon(
        unsigned char *y_ptr,
        int y_stride,
        const unsigned char *blimit) {
    y_ptr += 4;
    vp8_loop_filter_simple_vertical_edge_neon(y_ptr, y_stride, blimit);
    y_ptr += 4;
    vp8_loop_filter_simple_vertical_edge_neon(y_ptr, y_stride, blimit);
    y_ptr += 4;
    vp8_loop_filter_simple_vertical_edge_neon(y_ptr, y_stride, blimit);
    return;
}

void vp8_loop_filter_mbvs_neon(
        unsigned char *y_ptr,
        int y_stride,
        const unsigned char *blimit) {
    vp8_loop_filter_simple_vertical_edge_neon(y_ptr, y_stride, blimit);
    return;
}
