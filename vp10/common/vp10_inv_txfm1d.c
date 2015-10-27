/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp10/common/vp10_inv_txfm1d.h"
#if CONFIG_COEFFICIENT_RANGE_CHECKING
#define range_check(stage, input, buf, size, bit)                         \
  {                                                                       \
    int i, j;                                                             \
    for (i = 0; i < size; ++i) {                                          \
      int buf_bit = get_max_bit(abs(buf[i])) + 1;                         \
      if (buf_bit > bit) {                                                \
        printf("======== %s overflow ========\n", __func__);              \
        printf("stage: %d node: %d\n", stage, i);                         \
        printf("bit: %d buf_bit: %d buf[i]: %d\n", bit, buf_bit, buf[i]); \
        printf("input:\n");                                               \
        for (j = 0; j < size; j++) {                                      \
          printf("%d,", input[j]);                                        \
        }                                                                 \
        printf("\n");                                                     \
        assert(0, "vp10_inv_txfm1d.c: range_check overflow");             \
      }                                                                   \
    }                                                                     \
  }
#else
#define range_check(stage, input, buf, size, bit) \
  {                                               \
    (void) stage;                                 \
    (void) input;                                 \
    (void) buf;                                   \
    (void) size;                                  \
    (void) bit;                                   \
  }
#endif

void vp10_idct4_new(const int32_t *input, int32_t *output,
                    const int8_t *cos_bit, const int8_t *stage_range) {
  const int32_t size = 4;
  const int32_t *cospi;

  int32_t stage = 0;
  int32_t *bf0, *bf1;
  int32_t step[4];

  // stage 0;
  range_check(stage, input, input, size, stage_range[stage]);

  // stage 1;
  stage++;
  bf1 = output;
  bf1[0] = input[0];
  bf1[1] = input[2];
  bf1[2] = input[1];
  bf1[3] = input[3];
  range_check(stage, input, bf1, size, stage_range[stage]);

  // stage 2
  stage++;
  cospi = cospi_arr[cos_bit[stage] - cos_bit_min];
  bf0 = output;
  bf1 = step;
  bf1[0] = half_btf(cospi[32], bf0[0], cospi[32], bf0[1], cos_bit[stage]);
  bf1[1] = half_btf(cospi[32], bf0[0], -cospi[32], bf0[1], cos_bit[stage]);
  bf1[2] = half_btf(cospi[48], bf0[2], -cospi[16], bf0[3], cos_bit[stage]);
  bf1[3] = half_btf(cospi[16], bf0[2], cospi[48], bf0[3], cos_bit[stage]);
  range_check(stage, input, bf1, size, stage_range[stage]);

  // stage 3
  stage++;
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[3];
  bf1[1] = bf0[1] + bf0[2];
  bf1[2] = bf0[1] - bf0[2];
  bf1[3] = bf0[0] - bf0[3];
  range_check(stage, input, bf1, size, stage_range[stage]);
}
