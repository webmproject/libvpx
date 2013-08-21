/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vp9_rtcd.h"
#include "vp9/common/vp9_common.h"

extern void vp9_short_idct16x16_add_neon_pass1(int16_t *input,
                                               int16_t *output,
                                               int output_stride);
extern void vp9_short_idct16x16_add_neon_pass2(int16_t *src,
                                               int16_t *output,
                                               int16_t *pass1Output,
                                               int16_t skip_adding,
                                               uint8_t *dest,
                                               int dest_stride);
extern void vp9_short_idct10_16x16_add_neon_pass1(int16_t *input,
                                               int16_t *output,
                                               int output_stride);
extern void vp9_short_idct10_16x16_add_neon_pass2(int16_t *src,
                                               int16_t *output,
                                               int16_t *pass1Output,
                                               int16_t skip_adding,
                                               uint8_t *dest,
                                               int dest_stride);
extern void save_registers();
extern void restore_registers();


void vp9_short_idct16x16_add_neon(int16_t *input,
                                  uint8_t *dest, int dest_stride) {
  int16_t pass1_output[16*16] = {0};
  int16_t row_idct_output[16*16] = {0};

  // save d8-d15 register values.
  save_registers();

  /* Parallel idct on the upper 8 rows */
  // First pass processes even elements 0, 2, 4, 6, 8, 10, 12, 14 and save the
  // stage 6 result in pass1_output.
  vp9_short_idct16x16_add_neon_pass1(input, pass1_output, 8);

  // Second pass processes odd elements 1, 3, 5, 7, 9, 11, 13, 15 and combines
  // with result in pass1(pass1_output) to calculate final result in stage 7
  // which will be saved into row_idct_output.
  vp9_short_idct16x16_add_neon_pass2(input+1,
                                     row_idct_output,
                                     pass1_output,
                                     0,
                                     dest,
                                     dest_stride);

  /* Parallel idct on the lower 8 rows */
  // First pass processes even elements 0, 2, 4, 6, 8, 10, 12, 14 and save the
  // stage 6 result in pass1_output.
  vp9_short_idct16x16_add_neon_pass1(input+8*16, pass1_output, 8);

  // Second pass processes odd elements 1, 3, 5, 7, 9, 11, 13, 15 and combines
  // with result in pass1(pass1_output) to calculate final result in stage 7
  // which will be saved into row_idct_output.
  vp9_short_idct16x16_add_neon_pass2(input+8*16+1,
                                     row_idct_output+8,
                                     pass1_output,
                                     0,
                                     dest,
                                     dest_stride);

  /* Parallel idct on the left 8 columns */
  // First pass processes even elements 0, 2, 4, 6, 8, 10, 12, 14 and save the
  // stage 6 result in pass1_output.
  vp9_short_idct16x16_add_neon_pass1(row_idct_output, pass1_output, 8);

  // Second pass processes odd elements 1, 3, 5, 7, 9, 11, 13, 15 and combines
  // with result in pass1(pass1_output) to calculate final result in stage 7.
  // Then add the result to the destination data.
  vp9_short_idct16x16_add_neon_pass2(row_idct_output+1,
                                     row_idct_output,
                                     pass1_output,
                                     1,
                                     dest,
                                     dest_stride);

  /* Parallel idct on the right 8 columns */
  // First pass processes even elements 0, 2, 4, 6, 8, 10, 12, 14 and save the
  // stage 6 result in pass1_output.
  vp9_short_idct16x16_add_neon_pass1(row_idct_output+8*16, pass1_output, 8);

  // Second pass processes odd elements 1, 3, 5, 7, 9, 11, 13, 15 and combines
  // with result in pass1(pass1_output) to calculate final result in stage 7.
  // Then add the result to the destination data.
  vp9_short_idct16x16_add_neon_pass2(row_idct_output+8*16+1,
                                     row_idct_output+8,
                                     pass1_output,
                                     1,
                                     dest+8,
                                     dest_stride);

  // restore d8-d15 register values.
  restore_registers();

  return;
}

void vp9_short_idct10_16x16_add_neon(int16_t *input,
                                  uint8_t *dest, int dest_stride) {
  int16_t pass1_output[16*16] = {0};
  int16_t row_idct_output[16*16] = {0};

  // save d8-d15 register values.
  save_registers();

  /* Parallel idct on the upper 8 rows */
  // First pass processes even elements 0, 2, 4, 6, 8, 10, 12, 14 and save the
  // stage 6 result in pass1_output.
  vp9_short_idct10_16x16_add_neon_pass1(input, pass1_output, 8);

  // Second pass processes odd elements 1, 3, 5, 7, 9, 11, 13, 15 and combines
  // with result in pass1(pass1_output) to calculate final result in stage 7
  // which will be saved into row_idct_output.
  vp9_short_idct10_16x16_add_neon_pass2(input+1,
                                        row_idct_output,
                                        pass1_output,
                                        0,
                                        dest,
                                        dest_stride);

  /* Skip Parallel idct on the lower 8 rows as they are all 0s */

  /* Parallel idct on the left 8 columns */
  // First pass processes even elements 0, 2, 4, 6, 8, 10, 12, 14 and save the
  // stage 6 result in pass1_output.
  vp9_short_idct16x16_add_neon_pass1(row_idct_output, pass1_output, 8);

  // Second pass processes odd elements 1, 3, 5, 7, 9, 11, 13, 15 and combines
  // with result in pass1(pass1_output) to calculate final result in stage 7.
  // Then add the result to the destination data.
  vp9_short_idct16x16_add_neon_pass2(row_idct_output+1,
                                     row_idct_output,
                                     pass1_output,
                                     1,
                                     dest,
                                     dest_stride);

  /* Parallel idct on the right 8 columns */
  // First pass processes even elements 0, 2, 4, 6, 8, 10, 12, 14 and save the
  // stage 6 result in pass1_output.
  vp9_short_idct16x16_add_neon_pass1(row_idct_output+8*16, pass1_output, 8);

  // Second pass processes odd elements 1, 3, 5, 7, 9, 11, 13, 15 and combines
  // with result in pass1(pass1_output) to calculate final result in stage 7.
  // Then add the result to the destination data.
  vp9_short_idct16x16_add_neon_pass2(row_idct_output+8*16+1,
                                     row_idct_output+8,
                                     pass1_output,
                                     1,
                                     dest+8,
                                     dest_stride);

  // restore d8-d15 register values.
  restore_registers();

  return;
}
