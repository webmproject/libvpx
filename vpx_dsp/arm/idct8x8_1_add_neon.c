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

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/inv_txfm.h"
#include "vpx_ports/mem.h"

void vpx_idct8x8_1_add_neon(const tran_low_t *input, uint8_t *dest,
                            int stride) {
  int i;
  const int16_t out0 = dct_const_round_shift(input[0] * cospi_16_64);
  const int16_t out1 = dct_const_round_shift(out0 * cospi_16_64);
  const int16_t out2 = ROUND_POWER_OF_TWO(out1, 5);
  const int16x8_t dc = vdupq_n_s16(out2);
  const uint16x8_t dc_u16 = vreinterpretq_u16_s16(dc);
  const uint8_t *dst = dest;
  uint8x8_t d0, d1, d2, d3;
  uint16x8_t d0_u16, d1_u16, d2_u16, d3_u16;

  for (i = 0; i < 2; i++) {
    d0 = vld1_u8(dst);
    dst += stride;
    d1 = vld1_u8(dst);
    dst += stride;
    d2 = vld1_u8(dst);
    dst += stride;
    d3 = vld1_u8(dst);
    dst += stride;

    d0_u16 = vaddw_u8(dc_u16, d0);
    d1_u16 = vaddw_u8(dc_u16, d1);
    d2_u16 = vaddw_u8(dc_u16, d2);
    d3_u16 = vaddw_u8(dc_u16, d3);

    d0 = vqmovun_s16(vreinterpretq_s16_u16(d0_u16));
    d1 = vqmovun_s16(vreinterpretq_s16_u16(d1_u16));
    d2 = vqmovun_s16(vreinterpretq_s16_u16(d2_u16));
    d3 = vqmovun_s16(vreinterpretq_s16_u16(d3_u16));

    vst1_u8(dest, d0);
    dest += stride;
    vst1_u8(dest, d1);
    dest += stride;
    vst1_u8(dest, d2);
    dest += stride;
    vst1_u8(dest, d3);
    dest += stride;
  }
}
