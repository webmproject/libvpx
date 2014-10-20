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
#include "vp8/encoder/block.h"
#include "vpx_mem/vpx_mem.h"

static const uint16_t inv_zig_zag[16] = {
    0x0001, 0x0002, 0x0006, 0x0007,
    0x0003, 0x0005, 0x0008, 0x000d,
    0x0004, 0x0009, 0x000c, 0x000e,
    0x000a, 0x000b, 0x000f, 0x0010
};

void vp8_fast_quantize_b_neon(BLOCK *b, BLOCKD *d) {
    const int16x8_t one_q = vdupq_n_s16(0xff),
                    z0 = vld1q_s16(b->coeff),
                    z1 = vld1q_s16(b->coeff + 8),
                    round0 = vld1q_s16(b->round),
                    round1 = vld1q_s16(b->round + 8),
                    quant0 = vld1q_s16(b->quant_fast),
                    quant1 = vld1q_s16(b->quant_fast + 8),
                    dequant0 = vld1q_s16(d->dequant),
                    dequant1 = vld1q_s16(d->dequant + 8);
    const uint16x8_t zig_zag0 = vld1q_u16(inv_zig_zag),
                     zig_zag1 = vld1q_u16(inv_zig_zag + 8);
    int16x8_t x0, x1, sz0, sz1, y0, y1;
    uint16x8_t eob0, eob1;
    uint16x4_t eob_d16;
    uint32x2_t eob_d32;
    uint32x4_t eob_q32;

    /* sign of z: z >> 15 */
    sz0 = vshrq_n_s16(z0, 15);
    sz1 = vshrq_n_s16(z1, 15);

    /* x = abs(z) */
    x0 = vabsq_s16(z0);
    x1 = vabsq_s16(z1);

    /* x += round */
    x0 = vaddq_s16(x0, round0);
    x1 = vaddq_s16(x1, round1);

    /* y = 2 * (x * quant) >> 16 */
    y0 = vqdmulhq_s16(x0, quant0);
    y1 = vqdmulhq_s16(x1, quant1);

    /* Compensate for doubling in vqdmulhq */
    y0 = vshrq_n_s16(y0, 1);
    y1 = vshrq_n_s16(y1, 1);

    /* Restore sign bit */
    y0 = veorq_s16(y0, sz0);
    y1 = veorq_s16(y1, sz1);
    x0 = vsubq_s16(y0, sz0);
    x1 = vsubq_s16(y1, sz1);

    /* find non-zero elements */
    eob0 = vtstq_s16(x0, one_q);
    eob1 = vtstq_s16(x1, one_q);

    /* mask zig zag */
    eob0 = vandq_u16(eob0, zig_zag0);
    eob1 = vandq_u16(eob1, zig_zag1);

    /* select the largest value */
    eob0 = vmaxq_u16(eob0, eob1);
    eob_d16 = vmax_u16(vget_low_u16(eob0), vget_high_u16(eob0));
    eob_q32 = vmovl_u16(eob_d16);
    eob_d32 = vmax_u32(vget_low_u32(eob_q32), vget_high_u32(eob_q32));
    eob_d32 = vpmax_u32(eob_d32, eob_d32);

    /* qcoeff = x */
    vst1q_s16(d->qcoeff, x0);
    vst1q_s16(d->qcoeff + 8, x1);

    /* dqcoeff = x * dequant */
    vst1q_s16(d->dqcoeff, vmulq_s16(dequant0, x0));
    vst1q_s16(d->dqcoeff + 8, vmulq_s16(dequant1, x1));

    vst1_lane_s8((int8_t *)d->eob, vreinterpret_s8_u32(eob_d32), 0);
}

void vp8_fast_quantize_b_pair_neon(BLOCK *b0, BLOCK *b1,
                                   BLOCKD *d0, BLOCKD *d1) {
    const int16x8_t one_q = vdupq_n_s16(0xff),
                    b0_z0 = vld1q_s16(b0->coeff),
                    b0_z1 = vld1q_s16(b0->coeff + 8),
                    b0_round0 = vld1q_s16(b0->round),
                    b0_round1 = vld1q_s16(b0->round + 8),
                    b0_quant0 = vld1q_s16(b0->quant_fast),
                    b0_quant1 = vld1q_s16(b0->quant_fast + 8),
                    d0_dequant0 = vld1q_s16(d0->dequant),
                    d0_dequant1 = vld1q_s16(d0->dequant + 8),
                    b1_z0 = vld1q_s16(b1->coeff),
                    b1_z1 = vld1q_s16(b1->coeff + 8),
                    b1_round0 = vld1q_s16(b1->round),
                    b1_round1 = vld1q_s16(b1->round + 8),
                    b1_quant0 = vld1q_s16(b1->quant_fast),
                    b1_quant1 = vld1q_s16(b1->quant_fast + 8),
                    d1_dequant0 = vld1q_s16(d1->dequant),
                    d1_dequant1 = vld1q_s16(d1->dequant + 8);
    const uint16x8_t zig_zag0 = vld1q_u16(inv_zig_zag),
                     zig_zag1 = vld1q_u16(inv_zig_zag + 8);
    int16x8_t b0_x0, b0_x1, b0_sz0, b0_sz1, b0_y0, b0_y1,
              b1_x0, b1_x1, b1_sz0, b1_sz1, b1_y0, b1_y1;
    uint16x8_t b0_eob0, b0_eob1,
               b1_eob0, b1_eob1;
    uint16x4_t b0_eob_d16, b1_eob_d16;
    uint32x2_t b0_eob_d32, b1_eob_d32;
    uint32x4_t b0_eob_q32, b1_eob_q32;

    /* sign of z: z >> 15 */
    b0_sz0 = vshrq_n_s16(b0_z0, 15);
    b0_sz1 = vshrq_n_s16(b0_z1, 15);
    b1_sz0 = vshrq_n_s16(b1_z0, 15);
    b1_sz1 = vshrq_n_s16(b1_z1, 15);

    /* x = abs(z) */
    b0_x0 = vabsq_s16(b0_z0);
    b0_x1 = vabsq_s16(b0_z1);
    b1_x0 = vabsq_s16(b1_z0);
    b1_x1 = vabsq_s16(b1_z1);

    /* x += round */
    b0_x0 = vaddq_s16(b0_x0, b0_round0);
    b0_x1 = vaddq_s16(b0_x1, b0_round1);
    b1_x0 = vaddq_s16(b1_x0, b1_round0);
    b1_x1 = vaddq_s16(b1_x1, b1_round1);

    /* y = 2 * (x * quant) >> 16 */
    b0_y0 = vqdmulhq_s16(b0_x0, b0_quant0);
    b0_y1 = vqdmulhq_s16(b0_x1, b0_quant1);
    b1_y0 = vqdmulhq_s16(b1_x0, b1_quant0);
    b1_y1 = vqdmulhq_s16(b1_x1, b1_quant1);

    /* Compensate for doubling in vqdmulhq */
    b0_y0 = vshrq_n_s16(b0_y0, 1);
    b0_y1 = vshrq_n_s16(b0_y1, 1);
    b1_y0 = vshrq_n_s16(b1_y0, 1);
    b1_y1 = vshrq_n_s16(b1_y1, 1);

    /* Restore sign bit */
    b0_y0 = veorq_s16(b0_y0, b0_sz0);
    b0_y1 = veorq_s16(b0_y1, b0_sz1);
    b0_x0 = vsubq_s16(b0_y0, b0_sz0);
    b0_x1 = vsubq_s16(b0_y1, b0_sz1);
    b1_y0 = veorq_s16(b1_y0, b1_sz0);
    b1_y1 = veorq_s16(b1_y1, b1_sz1);
    b1_x0 = vsubq_s16(b1_y0, b1_sz0);
    b1_x1 = vsubq_s16(b1_y1, b1_sz1);

    /* find non-zero elements */
    b0_eob0 = vtstq_s16(b0_x0, one_q);
    b0_eob1 = vtstq_s16(b0_x1, one_q);
    b1_eob0 = vtstq_s16(b1_x0, one_q);
    b1_eob1 = vtstq_s16(b1_x1, one_q);

    /* mask zig zag */
    b0_eob0 = vandq_u16(b0_eob0, zig_zag0);
    b0_eob1 = vandq_u16(b0_eob1, zig_zag1);
    b1_eob0 = vandq_u16(b1_eob0, zig_zag0);
    b1_eob1 = vandq_u16(b1_eob1, zig_zag1);

    /* select the largest value */
    b0_eob0 = vmaxq_u16(b0_eob0, b0_eob1);
    b0_eob_d16 = vmax_u16(vget_low_u16(b0_eob0),
                          vget_high_u16(b0_eob0));
    b0_eob_q32 = vmovl_u16(b0_eob_d16);
    b0_eob_d32 = vmax_u32(vget_low_u32(b0_eob_q32),
                          vget_high_u32(b0_eob_q32));
    b0_eob_d32 = vpmax_u32(b0_eob_d32, b0_eob_d32);

    b1_eob0 = vmaxq_u16(b1_eob0, b1_eob1);
    b1_eob_d16 = vmax_u16(vget_low_u16(b1_eob0),
                          vget_high_u16(b1_eob0));
    b1_eob_q32 = vmovl_u16(b1_eob_d16);
    b1_eob_d32 = vmax_u32(vget_low_u32(b1_eob_q32),
                          vget_high_u32(b1_eob_q32));
    b1_eob_d32 = vpmax_u32(b1_eob_d32, b1_eob_d32);

    /* qcoeff = x */
    vst1q_s16(d0->qcoeff, b0_x0);
    vst1q_s16(d0->qcoeff + 8, b0_x1);
    vst1q_s16(d1->qcoeff, b1_x0);
    vst1q_s16(d1->qcoeff + 8, b1_x1);

    /* dqcoeff = x * dequant */
    vst1q_s16(d0->dqcoeff, vmulq_s16(d0_dequant0, b0_x0));
    vst1q_s16(d0->dqcoeff + 8, vmulq_s16(d0_dequant1, b0_x1));
    vst1q_s16(d1->dqcoeff, vmulq_s16(d1_dequant0, b1_x0));
    vst1q_s16(d1->dqcoeff + 8, vmulq_s16(d1_dequant1, b1_x1));

    vst1_lane_s8((int8_t *)d0->eob, vreinterpret_s8_u32(b0_eob_d32), 0);
    vst1_lane_s8((int8_t *)d1->eob, vreinterpret_s8_u32(b1_eob_d32), 0);
    return;
}
