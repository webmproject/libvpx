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
#include "vpx_dsp/arm/fdct_neon.h"
#include "vpx_dsp/arm/fdct32x32_neon.h"

// Most gcc 4.9 distributions outside of Android do not generate correct code
// for this function.
#if !defined(__clang__) && !defined(__ANDROID__) && defined(__GNUC__) && \
    __GNUC__ == 4 && __GNUC_MINOR__ <= 9

void vpx_fdct32x32_neon(const int16_t *input, tran_low_t *output, int stride) {
  vpx_fdct32x32_c(input, output, stride);
}

void vpx_fdct32x32_rd_neon(const int16_t *input, tran_low_t *output,
                           int stride) {
  vpx_fdct32x32_rd_c(input, output, stride);
}

#else

void vpx_fdct32x32_neon(const int16_t *input, tran_low_t *output, int stride) {
  int16x8_t temp0[32];
  int16x8_t temp1[32];
  int16x8_t temp2[32];
  int16x8_t temp3[32];
  int16x8_t temp4[32];
  int16x8_t temp5[32];

  // Process in 8x32 columns.
  load_cross(input, stride, temp0);
  scale_input(temp0, temp5);
  dct_body_first_pass(temp5, temp1);

  load_cross(input + 8, stride, temp0);
  scale_input(temp0, temp5);
  dct_body_first_pass(temp5, temp2);

  load_cross(input + 16, stride, temp0);
  scale_input(temp0, temp5);
  dct_body_first_pass(temp5, temp3);

  load_cross(input + 24, stride, temp0);
  scale_input(temp0, temp5);
  dct_body_first_pass(temp5, temp4);

  // Generate the top row by munging the first set of 8 from each one together.
  transpose_s16_8x8_new(&temp1[0], &temp0[0]);
  transpose_s16_8x8_new(&temp2[0], &temp0[8]);
  transpose_s16_8x8_new(&temp3[0], &temp0[16]);
  transpose_s16_8x8_new(&temp4[0], &temp0[24]);

  dct_body_second_pass(temp0, temp5);

  transpose_s16_8x8(&temp5[0], &temp5[1], &temp5[2], &temp5[3], &temp5[4],
                    &temp5[5], &temp5[6], &temp5[7]);
  transpose_s16_8x8(&temp5[8], &temp5[9], &temp5[10], &temp5[11], &temp5[12],
                    &temp5[13], &temp5[14], &temp5[15]);
  transpose_s16_8x8(&temp5[16], &temp5[17], &temp5[18], &temp5[19], &temp5[20],
                    &temp5[21], &temp5[22], &temp5[23]);
  transpose_s16_8x8(&temp5[24], &temp5[25], &temp5[26], &temp5[27], &temp5[28],
                    &temp5[29], &temp5[30], &temp5[31]);
  store(output, temp5);

  // Second row of 8x32.
  transpose_s16_8x8_new(&temp1[8], &temp0[0]);
  transpose_s16_8x8_new(&temp2[8], &temp0[8]);
  transpose_s16_8x8_new(&temp3[8], &temp0[16]);
  transpose_s16_8x8_new(&temp4[8], &temp0[24]);

  dct_body_second_pass(temp0, temp5);

  transpose_s16_8x8(&temp5[0], &temp5[1], &temp5[2], &temp5[3], &temp5[4],
                    &temp5[5], &temp5[6], &temp5[7]);
  transpose_s16_8x8(&temp5[8], &temp5[9], &temp5[10], &temp5[11], &temp5[12],
                    &temp5[13], &temp5[14], &temp5[15]);
  transpose_s16_8x8(&temp5[16], &temp5[17], &temp5[18], &temp5[19], &temp5[20],
                    &temp5[21], &temp5[22], &temp5[23]);
  transpose_s16_8x8(&temp5[24], &temp5[25], &temp5[26], &temp5[27], &temp5[28],
                    &temp5[29], &temp5[30], &temp5[31]);
  store(output + 8 * 32, temp5);

  // Third row of 8x32
  transpose_s16_8x8_new(&temp1[16], &temp0[0]);
  transpose_s16_8x8_new(&temp2[16], &temp0[8]);
  transpose_s16_8x8_new(&temp3[16], &temp0[16]);
  transpose_s16_8x8_new(&temp4[16], &temp0[24]);

  dct_body_second_pass(temp0, temp5);

  transpose_s16_8x8(&temp5[0], &temp5[1], &temp5[2], &temp5[3], &temp5[4],
                    &temp5[5], &temp5[6], &temp5[7]);
  transpose_s16_8x8(&temp5[8], &temp5[9], &temp5[10], &temp5[11], &temp5[12],
                    &temp5[13], &temp5[14], &temp5[15]);
  transpose_s16_8x8(&temp5[16], &temp5[17], &temp5[18], &temp5[19], &temp5[20],
                    &temp5[21], &temp5[22], &temp5[23]);
  transpose_s16_8x8(&temp5[24], &temp5[25], &temp5[26], &temp5[27], &temp5[28],
                    &temp5[29], &temp5[30], &temp5[31]);
  store(output + 16 * 32, temp5);

  // Final row of 8x32.
  transpose_s16_8x8_new(&temp1[24], &temp0[0]);
  transpose_s16_8x8_new(&temp2[24], &temp0[8]);
  transpose_s16_8x8_new(&temp3[24], &temp0[16]);
  transpose_s16_8x8_new(&temp4[24], &temp0[24]);

  dct_body_second_pass(temp0, temp5);

  transpose_s16_8x8(&temp5[0], &temp5[1], &temp5[2], &temp5[3], &temp5[4],
                    &temp5[5], &temp5[6], &temp5[7]);
  transpose_s16_8x8(&temp5[8], &temp5[9], &temp5[10], &temp5[11], &temp5[12],
                    &temp5[13], &temp5[14], &temp5[15]);
  transpose_s16_8x8(&temp5[16], &temp5[17], &temp5[18], &temp5[19], &temp5[20],
                    &temp5[21], &temp5[22], &temp5[23]);
  transpose_s16_8x8(&temp5[24], &temp5[25], &temp5[26], &temp5[27], &temp5[28],
                    &temp5[29], &temp5[30], &temp5[31]);
  store(output + 24 * 32, temp5);
}

void vpx_fdct32x32_rd_neon(const int16_t *input, tran_low_t *output,
                           int stride) {
  int16x8_t temp0[32];
  int16x8_t temp1[32];
  int16x8_t temp2[32];
  int16x8_t temp3[32];
  int16x8_t temp4[32];
  int16x8_t temp5[32];

  // Process in 8x32 columns.
  load_cross(input, stride, temp0);
  scale_input(temp0, temp5);
  dct_body_first_pass(temp5, temp1);

  load_cross(input + 8, stride, temp0);
  scale_input(temp0, temp5);
  dct_body_first_pass(temp5, temp2);

  load_cross(input + 16, stride, temp0);
  scale_input(temp0, temp5);
  dct_body_first_pass(temp5, temp3);

  load_cross(input + 24, stride, temp0);
  scale_input(temp0, temp5);
  dct_body_first_pass(temp5, temp4);

  // Generate the top row by munging the first set of 8 from each one together.
  transpose_s16_8x8_new(&temp1[0], &temp0[0]);
  transpose_s16_8x8_new(&temp2[0], &temp0[8]);
  transpose_s16_8x8_new(&temp3[0], &temp0[16]);
  transpose_s16_8x8_new(&temp4[0], &temp0[24]);

  dct_body_second_pass_rd(temp0, temp5);

  transpose_s16_8x8(&temp5[0], &temp5[1], &temp5[2], &temp5[3], &temp5[4],
                    &temp5[5], &temp5[6], &temp5[7]);
  transpose_s16_8x8(&temp5[8], &temp5[9], &temp5[10], &temp5[11], &temp5[12],
                    &temp5[13], &temp5[14], &temp5[15]);
  transpose_s16_8x8(&temp5[16], &temp5[17], &temp5[18], &temp5[19], &temp5[20],
                    &temp5[21], &temp5[22], &temp5[23]);
  transpose_s16_8x8(&temp5[24], &temp5[25], &temp5[26], &temp5[27], &temp5[28],
                    &temp5[29], &temp5[30], &temp5[31]);
  store(output, temp5);

  // Second row of 8x32.
  transpose_s16_8x8_new(&temp1[8], &temp0[0]);
  transpose_s16_8x8_new(&temp2[8], &temp0[8]);
  transpose_s16_8x8_new(&temp3[8], &temp0[16]);
  transpose_s16_8x8_new(&temp4[8], &temp0[24]);

  dct_body_second_pass_rd(temp0, temp5);

  transpose_s16_8x8(&temp5[0], &temp5[1], &temp5[2], &temp5[3], &temp5[4],
                    &temp5[5], &temp5[6], &temp5[7]);
  transpose_s16_8x8(&temp5[8], &temp5[9], &temp5[10], &temp5[11], &temp5[12],
                    &temp5[13], &temp5[14], &temp5[15]);
  transpose_s16_8x8(&temp5[16], &temp5[17], &temp5[18], &temp5[19], &temp5[20],
                    &temp5[21], &temp5[22], &temp5[23]);
  transpose_s16_8x8(&temp5[24], &temp5[25], &temp5[26], &temp5[27], &temp5[28],
                    &temp5[29], &temp5[30], &temp5[31]);
  store(output + 8 * 32, temp5);

  // Third row of 8x32
  transpose_s16_8x8_new(&temp1[16], &temp0[0]);
  transpose_s16_8x8_new(&temp2[16], &temp0[8]);
  transpose_s16_8x8_new(&temp3[16], &temp0[16]);
  transpose_s16_8x8_new(&temp4[16], &temp0[24]);

  dct_body_second_pass_rd(temp0, temp5);

  transpose_s16_8x8(&temp5[0], &temp5[1], &temp5[2], &temp5[3], &temp5[4],
                    &temp5[5], &temp5[6], &temp5[7]);
  transpose_s16_8x8(&temp5[8], &temp5[9], &temp5[10], &temp5[11], &temp5[12],
                    &temp5[13], &temp5[14], &temp5[15]);
  transpose_s16_8x8(&temp5[16], &temp5[17], &temp5[18], &temp5[19], &temp5[20],
                    &temp5[21], &temp5[22], &temp5[23]);
  transpose_s16_8x8(&temp5[24], &temp5[25], &temp5[26], &temp5[27], &temp5[28],
                    &temp5[29], &temp5[30], &temp5[31]);
  store(output + 16 * 32, temp5);

  // Final row of 8x32.
  transpose_s16_8x8_new(&temp1[24], &temp0[0]);
  transpose_s16_8x8_new(&temp2[24], &temp0[8]);
  transpose_s16_8x8_new(&temp3[24], &temp0[16]);
  transpose_s16_8x8_new(&temp4[24], &temp0[24]);

  dct_body_second_pass_rd(temp0, temp5);

  transpose_s16_8x8(&temp5[0], &temp5[1], &temp5[2], &temp5[3], &temp5[4],
                    &temp5[5], &temp5[6], &temp5[7]);
  transpose_s16_8x8(&temp5[8], &temp5[9], &temp5[10], &temp5[11], &temp5[12],
                    &temp5[13], &temp5[14], &temp5[15]);
  transpose_s16_8x8(&temp5[16], &temp5[17], &temp5[18], &temp5[19], &temp5[20],
                    &temp5[21], &temp5[22], &temp5[23]);
  transpose_s16_8x8(&temp5[24], &temp5[25], &temp5[26], &temp5[27], &temp5[28],
                    &temp5[29], &temp5[30], &temp5[31]);
  store(output + 24 * 32, temp5);
}
#endif  // !defined(__clang__) && !defined(__ANDROID__) && defined(__GNUC__) &&
        // __GNUC__ == 4 && __GNUC_MINOR__ <= 9
