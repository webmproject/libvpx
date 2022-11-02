/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
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
#include "vpx_dsp/txfm_common.h"
#include "vpx_dsp/arm/mem_neon.h"
#include "vpx_dsp/arm/transpose_neon.h"
#include "vpx_dsp/arm/fdct16x16_neon.h"

// Some builds of gcc 4.9.2 and .3 have trouble with some of the inline
// functions.
#if !defined(__clang__) && !defined(__ANDROID__) && defined(__GNUC__) && \
    __GNUC__ == 4 && __GNUC_MINOR__ == 9 && __GNUC_PATCHLEVEL__ < 4

void vpx_fdct16x16_neon(const int16_t *input, tran_low_t *output, int stride) {
  vpx_fdct16x16_c(input, output, stride);
}

#else

void vpx_fdct16x16_neon(const int16_t *input, tran_low_t *output, int stride) {
  int16x8_t temp0[16];
  int16x8_t temp1[16];
  int16x8_t temp2[16];
  int16x8_t temp3[16];

  // Left half.
  load_cross(input, stride, temp0);
  scale_input(temp0, temp1);
  vpx_fdct8x16_body(temp1, temp0);

  // Right half.
  load_cross(input + 8, stride, temp1);
  scale_input(temp1, temp2);
  vpx_fdct8x16_body(temp2, temp1);

  // Transpose top left and top right quarters into one contiguous location to
  // process to the top half.

  transpose_s16_8x8_new(&temp0[0], &temp2[0]);
  transpose_s16_8x8_new(&temp1[0], &temp2[8]);
  partial_round_shift(temp2);
  cross_input(temp2, temp3);
  vpx_fdct8x16_body(temp3, temp2);
  transpose_s16_8x8(&temp2[0], &temp2[1], &temp2[2], &temp2[3], &temp2[4],
                    &temp2[5], &temp2[6], &temp2[7]);
  transpose_s16_8x8(&temp2[8], &temp2[9], &temp2[10], &temp2[11], &temp2[12],
                    &temp2[13], &temp2[14], &temp2[15]);
  store(output, temp2);
  store(output + 8, temp2 + 8);
  output += 8 * 16;

  // Transpose bottom left and bottom right quarters into one contiguous
  // location to process to the bottom half.
  transpose_s16_8x8_new(&temp0[8], &temp1[0]);

  transpose_s16_8x8(&temp1[8], &temp1[9], &temp1[10], &temp1[11], &temp1[12],
                    &temp1[13], &temp1[14], &temp1[15]);
  partial_round_shift(temp1);
  cross_input(temp1, temp0);
  vpx_fdct8x16_body(temp0, temp1);
  transpose_s16_8x8(&temp1[0], &temp1[1], &temp1[2], &temp1[3], &temp1[4],
                    &temp1[5], &temp1[6], &temp1[7]);
  transpose_s16_8x8(&temp1[8], &temp1[9], &temp1[10], &temp1[11], &temp1[12],
                    &temp1[13], &temp1[14], &temp1[15]);
  store(output, temp1);
  store(output + 8, temp1 + 8);
}

#if CONFIG_VP9_HIGHBITDEPTH

void vpx_highbd_fdct16x16_neon(const int16_t *input, tran_low_t *output,
                               int stride) {
  int16x8_t temp0[16];
  int32x4_t left1[16], left2[16], left3[16], left4[16], right1[16], right2[16],
      right3[16], right4[16];

  // Left half.
  load_cross(input, stride, temp0);
  highbd_scale_input(temp0, left1, right1);
  vpx_highbd_fdct8x16_body(left1, right1);

  // right half.
  load_cross(input + 8, stride, temp0);
  highbd_scale_input(temp0, left2, right2);
  vpx_highbd_fdct8x16_body(left2, right2);

  // Transpose top left and top right quarters into one contiguous location to
  // process to the top half.

  transpose_s32_8x8_2(left1, right1, left3, right3);
  transpose_s32_8x8_2(left2, right2, left3 + 8, right3 + 8);
  transpose_s32_8x8_2(left1 + 8, right1 + 8, left4, right4);
  transpose_s32_8x8_2(left2 + 8, right2 + 8, left4 + 8, right4 + 8);

  highbd_partial_round_shift(left3, right3);
  highbd_cross_input(left3, right3, left1, right1);
  vpx_highbd_fdct8x16_body(left1, right1);

  // Transpose bottom left and bottom right quarters into one contiguous
  // location to process to the bottom half.

  highbd_partial_round_shift(left4, right4);
  highbd_cross_input(left4, right4, left2, right2);
  vpx_highbd_fdct8x16_body(left2, right2);

  transpose_s32_8x8_2(left1, right1, left3, right3);
  transpose_s32_8x8_2(left2, right2, left3 + 8, right3 + 8);
  transpose_s32_8x8_2(left1 + 8, right1 + 8, left4, right4);
  transpose_s32_8x8_2(left2 + 8, right2 + 8, left4 + 8, right4 + 8);
  store16_s32(output, left3);
  output += 4;
  store16_s32(output, right3);
  output += 4;

  store16_s32(output, left4);
  output += 4;
  store16_s32(output, right4);
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

#endif  // !defined(__clang__) && !defined(__ANDROID__) && defined(__GNUC__) &&
        // __GNUC__ == 4 && __GNUC_MINOR__ == 9 && __GNUC_PATCHLEVEL__ < 4
