/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>

#include "./vp10_rtcd.h"
#include "./vpx_dsp_rtcd.h"
#include "vp10/common/blockd.h"
#include "vp10/common/enums.h"
#include "vp10/common/idct.h"
#include "vpx_dsp/inv_txfm.h"
#include "vpx_ports/mem.h"

#if CONFIG_EXT_TX
void idst4_c(const tran_low_t *input, tran_low_t *output) {
  // {sin(pi/5), sin(pi*2/5)} * sqrt(2/5) * sqrt(2)
  static const int32_t sinvalue_lookup[] = {
    141124871, 228344838,
  };
  int64_t sum;
  int64_t s03 = (input[0] + input[3]);
  int64_t d03 = (input[0] - input[3]);
  int64_t s12 = (input[1] + input[2]);
  int64_t d12 = (input[1] - input[2]);
  sum = s03 * sinvalue_lookup[0] + s12 * sinvalue_lookup[1];
  output[0] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = d03 * sinvalue_lookup[1] + d12 * sinvalue_lookup[0];
  output[1] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = s03 * sinvalue_lookup[1] - s12 * sinvalue_lookup[0];
  output[2] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = d03 * sinvalue_lookup[0] - d12 * sinvalue_lookup[1];
  output[3] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
}

void idst8_c(const tran_low_t *input, tran_low_t *output) {
  // {sin(pi/9), sin(pi*2/9), ..., sin(pi*4/9)} * sqrt(2/9) * 2
  static const int32_t sinvalue_lookup[] = {
    86559612, 162678858, 219176632, 249238470
  };
  int64_t sum;
  int64_t s07 = (input[0] + input[7]);
  int64_t d07 = (input[0] - input[7]);
  int64_t s16 = (input[1] + input[6]);
  int64_t d16 = (input[1] - input[6]);
  int64_t s25 = (input[2] + input[5]);
  int64_t d25 = (input[2] - input[5]);
  int64_t s34 = (input[3] + input[4]);
  int64_t d34 = (input[3] - input[4]);
  sum = s07 * sinvalue_lookup[0] + s16 * sinvalue_lookup[1] +
        s25 * sinvalue_lookup[2] + s34 * sinvalue_lookup[3];
  output[0] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = d07 * sinvalue_lookup[1] + d16 * sinvalue_lookup[3] +
        d25 * sinvalue_lookup[2] + d34 * sinvalue_lookup[0];
  output[1] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = (s07 + s16 - s34)* sinvalue_lookup[2];
  output[2] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = d07 * sinvalue_lookup[3] + d16 * sinvalue_lookup[0] -
        d25 * sinvalue_lookup[2] - d34 * sinvalue_lookup[1];
  output[3] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = s07 * sinvalue_lookup[3] - s16 * sinvalue_lookup[0] -
        s25 * sinvalue_lookup[2] + s34 * sinvalue_lookup[1];
  output[4] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = (d07 - d16 + d34)* sinvalue_lookup[2];
  output[5] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = s07 * sinvalue_lookup[1] - s16 * sinvalue_lookup[3] +
        s25 * sinvalue_lookup[2] - s34 * sinvalue_lookup[0];
  output[6] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = d07 * sinvalue_lookup[0] - d16 * sinvalue_lookup[1] +
        d25 * sinvalue_lookup[2] - d34 * sinvalue_lookup[3];
  output[7] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
}

void idst16_c(const tran_low_t *input, tran_low_t *output) {
  // {sin(pi/17), sin(pi*2/17, ..., sin(pi*8/17)} * sqrt(2/17) * 2 * sqrt(2)
  static const int32_t sinvalue_lookup[] = {
    47852167, 94074787, 137093803, 175444254,
    207820161, 233119001, 250479254, 259309736
  };
  int64_t sum;
  int64_t s015 = (input[0] + input[15]);
  int64_t d015 = (input[0] - input[15]);
  int64_t s114 = (input[1] + input[14]);
  int64_t d114 = (input[1] - input[14]);
  int64_t s213 = (input[2] + input[13]);
  int64_t d213 = (input[2] - input[13]);
  int64_t s312 = (input[3] + input[12]);
  int64_t d312 = (input[3] - input[12]);
  int64_t s411 = (input[4] + input[11]);
  int64_t d411 = (input[4] - input[11]);
  int64_t s510 = (input[5] + input[10]);
  int64_t d510 = (input[5] - input[10]);
  int64_t s69  = (input[6] + input[9]);
  int64_t d69  = (input[6] - input[9]);
  int64_t s78  = (input[7] + input[8]);
  int64_t d78  = (input[7] - input[8]);
  sum = s015 * sinvalue_lookup[0] + s114 * sinvalue_lookup[1] +
        s213 * sinvalue_lookup[2] + s312 * sinvalue_lookup[3] +
        s411 * sinvalue_lookup[4] + s510 * sinvalue_lookup[5] +
        s69  * sinvalue_lookup[6] + s78  * sinvalue_lookup[7];
  output[0]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = d015 * sinvalue_lookup[1] + d114 * sinvalue_lookup[3] +
        d213 * sinvalue_lookup[5] + d312 * sinvalue_lookup[7] +
        d411 * sinvalue_lookup[6] + d510 * sinvalue_lookup[4] +
        d69  * sinvalue_lookup[2] + d78  * sinvalue_lookup[0];
  output[1]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = s015 * sinvalue_lookup[2] + s114 * sinvalue_lookup[5] +
        s213 * sinvalue_lookup[7] + s312 * sinvalue_lookup[4] +
        s411 * sinvalue_lookup[1] - s510 * sinvalue_lookup[0] -
        s69  * sinvalue_lookup[3] - s78  * sinvalue_lookup[6];
  output[2]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = d015 * sinvalue_lookup[3] + d114 * sinvalue_lookup[7] +
        d213 * sinvalue_lookup[4] + d312 * sinvalue_lookup[0] -
        d411 * sinvalue_lookup[2] - d510 * sinvalue_lookup[6] -
        d69  * sinvalue_lookup[5] - d78  * sinvalue_lookup[1];
  output[3]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = s015 * sinvalue_lookup[4] + s114 * sinvalue_lookup[6] +
        s213 * sinvalue_lookup[1] - s312 * sinvalue_lookup[2] -
        s411 * sinvalue_lookup[7] - s510 * sinvalue_lookup[3] +
        s69  * sinvalue_lookup[0] + s78  * sinvalue_lookup[5];
  output[4]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = d015 * sinvalue_lookup[5] + d114 * sinvalue_lookup[4] -
        d213 * sinvalue_lookup[0] - d312 * sinvalue_lookup[6] -
        d411 * sinvalue_lookup[3] + d510 * sinvalue_lookup[1] +
        d69  * sinvalue_lookup[7] + d78  * sinvalue_lookup[2];
  output[5]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = s015 * sinvalue_lookup[6] + s114 * sinvalue_lookup[2] -
        s213 * sinvalue_lookup[3] - s312 * sinvalue_lookup[5] +
        s411 * sinvalue_lookup[0] + s510 * sinvalue_lookup[7] +
        s69  * sinvalue_lookup[1] - s78  * sinvalue_lookup[4];
  output[6]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = d015 * sinvalue_lookup[7] + d114 * sinvalue_lookup[0] -
        d213 * sinvalue_lookup[6] - d312 * sinvalue_lookup[1] +
        d411 * sinvalue_lookup[5] + d510 * sinvalue_lookup[2] -
        d69  * sinvalue_lookup[4] - d78  * sinvalue_lookup[3];
  output[7]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = s015 * sinvalue_lookup[7] - s114 * sinvalue_lookup[0] -
        s213 * sinvalue_lookup[6] + s312 * sinvalue_lookup[1] +
        s411 * sinvalue_lookup[5] - s510 * sinvalue_lookup[2] -
        s69  * sinvalue_lookup[4] + s78  * sinvalue_lookup[3];
  output[8]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = d015 * sinvalue_lookup[6] - d114 * sinvalue_lookup[2] -
        d213 * sinvalue_lookup[3] + d312 * sinvalue_lookup[5] +
        d411 * sinvalue_lookup[0] - d510 * sinvalue_lookup[7] +
        d69  * sinvalue_lookup[1] + d78  * sinvalue_lookup[4];
  output[9]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = s015 * sinvalue_lookup[5] - s114 * sinvalue_lookup[4] -
        s213 * sinvalue_lookup[0] + s312 * sinvalue_lookup[6] -
        s411 * sinvalue_lookup[3] - s510 * sinvalue_lookup[1] +
        s69  * sinvalue_lookup[7] - s78  * sinvalue_lookup[2];
  output[10] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = d015 * sinvalue_lookup[4] - d114 * sinvalue_lookup[6] +
        d213 * sinvalue_lookup[1] + d312 * sinvalue_lookup[2] -
        d411 * sinvalue_lookup[7] + d510 * sinvalue_lookup[3] +
        d69  * sinvalue_lookup[0] - d78  * sinvalue_lookup[5];
  output[11] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = s015 * sinvalue_lookup[3] - s114 * sinvalue_lookup[7] +
        s213 * sinvalue_lookup[4] - s312 * sinvalue_lookup[0] -
        s411 * sinvalue_lookup[2] + s510 * sinvalue_lookup[6] -
        s69  * sinvalue_lookup[5] + s78  * sinvalue_lookup[1];
  output[12] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = d015 * sinvalue_lookup[2] - d114 * sinvalue_lookup[5] +
        d213 * sinvalue_lookup[7] - d312 * sinvalue_lookup[4] +
        d411 * sinvalue_lookup[1] + d510 * sinvalue_lookup[0] -
        d69  * sinvalue_lookup[3] + d78  * sinvalue_lookup[6];
  output[13] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = s015 * sinvalue_lookup[1] - s114 * sinvalue_lookup[3] +
        s213 * sinvalue_lookup[5] - s312 * sinvalue_lookup[7] +
        s411 * sinvalue_lookup[6] - s510 * sinvalue_lookup[4] +
        s69  * sinvalue_lookup[2] - s78  * sinvalue_lookup[0];
  output[14] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
  sum = d015 * sinvalue_lookup[0] - d114 * sinvalue_lookup[1] +
        d213 * sinvalue_lookup[2] - d312 * sinvalue_lookup[3] +
        d411 * sinvalue_lookup[4] - d510 * sinvalue_lookup[5] +
        d69  * sinvalue_lookup[6] - d78  * sinvalue_lookup[7];
  output[15] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), 8);
}

static void fliplr(uint8_t *dest, int stride, int l) {
  int i, j;
  for (i = 0; i < l; ++i) {
    for (j = 0; j < l / 2; ++j) {
      const uint8_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[i * stride + l - 1 - j];
      dest[i * stride + l - 1 - j] = tmp;
    }
  }
}

static void flipud(uint8_t *dest, int stride, int l) {
  int i, j;
  for (j = 0; j < l; ++j) {
    for (i = 0; i < l / 2; ++i) {
      const uint8_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[(l - 1 - i) * stride + j];
      dest[(l - 1 - i) * stride + j] = tmp;
    }
  }
}

static void fliplrud(uint8_t *dest, int stride, int l) {
  int i, j;
  for (i = 0; i < l / 2; ++i) {
    for (j = 0; j < l; ++j) {
      const uint8_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[(l - 1 - i) * stride + l - 1 - j];
      dest[(l - 1 - i) * stride + l - 1 - j] = tmp;
    }
  }
}

// Inverse identiy transform and add.
static void inv_idtx_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                           int bs) {
  int r, c;
  const int shift = bs < 32 ? 3 : 2;
  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c)
      dest[c] = clip_pixel_add(dest[c], input[c] >> shift);
    dest += stride;
    input += bs;
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
void highbd_idst4_c(const tran_low_t *input, tran_low_t *output, int bd) {
  // {sin(pi/5), sin(pi*2/5)} * sqrt(2/5) * sqrt(2)
  static const int32_t sinvalue_lookup[] = {
    141124871, 228344838,
  };
  int64_t sum;
  int64_t s03 = (input[0] + input[3]);
  int64_t d03 = (input[0] - input[3]);
  int64_t s12 = (input[1] + input[2]);
  int64_t d12 = (input[1] - input[2]);

#if !CONFIG_EMULATE_HARDWARE
  (void)bd;
#endif

  sum = s03 * sinvalue_lookup[0] + s12 * sinvalue_lookup[1];
  output[0] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = d03 * sinvalue_lookup[1] + d12 * sinvalue_lookup[0];
  output[1] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = s03 * sinvalue_lookup[1] - s12 * sinvalue_lookup[0];
  output[2] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = d03 * sinvalue_lookup[0] - d12 * sinvalue_lookup[1];
  output[3] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
}

void highbd_idst8_c(const tran_low_t *input, tran_low_t *output, int bd) {
  // {sin(pi/9), sin(pi*2/9), ..., sin(pi*4/9)} * sqrt(2/9) * 2
  static const int32_t sinvalue_lookup[] = {
    86559612, 162678858, 219176632, 249238470
  };
  int64_t sum;
  int64_t s07 = (input[0] + input[7]);
  int64_t d07 = (input[0] - input[7]);
  int64_t s16 = (input[1] + input[6]);
  int64_t d16 = (input[1] - input[6]);
  int64_t s25 = (input[2] + input[5]);
  int64_t d25 = (input[2] - input[5]);
  int64_t s34 = (input[3] + input[4]);
  int64_t d34 = (input[3] - input[4]);

#if !CONFIG_EMULATE_HARDWARE
  (void)bd;
#endif

  sum = s07 * sinvalue_lookup[0] + s16 * sinvalue_lookup[1] +
        s25 * sinvalue_lookup[2] + s34 * sinvalue_lookup[3];
  output[0] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = d07 * sinvalue_lookup[1] + d16 * sinvalue_lookup[3] +
        d25 * sinvalue_lookup[2] + d34 * sinvalue_lookup[0];
  output[1] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = (s07 + s16 - s34)* sinvalue_lookup[2];
  output[2] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = d07 * sinvalue_lookup[3] + d16 * sinvalue_lookup[0] -
        d25 * sinvalue_lookup[2] - d34 * sinvalue_lookup[1];
  output[3] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = s07 * sinvalue_lookup[3] - s16 * sinvalue_lookup[0] -
        s25 * sinvalue_lookup[2] + s34 * sinvalue_lookup[1];
  output[4] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = (d07 - d16 + d34)* sinvalue_lookup[2];
  output[5] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = s07 * sinvalue_lookup[1] - s16 * sinvalue_lookup[3] +
        s25 * sinvalue_lookup[2] - s34 * sinvalue_lookup[0];
  output[6] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = d07 * sinvalue_lookup[0] - d16 * sinvalue_lookup[1] +
        d25 * sinvalue_lookup[2] - d34 * sinvalue_lookup[3];
  output[7] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
}

void highbd_idst16_c(const tran_low_t *input, tran_low_t *output, int bd) {
  // {sin(pi/17), sin(pi*2/17, ..., sin(pi*8/17)} * sqrt(2/17) * 2 * sqrt(2)
  static const int32_t sinvalue_lookup[] = {
    47852167, 94074787, 137093803, 175444254,
    207820161, 233119001, 250479254, 259309736
  };
  int64_t sum;
  int64_t s015 = (input[0] + input[15]);
  int64_t d015 = (input[0] - input[15]);
  int64_t s114 = (input[1] + input[14]);
  int64_t d114 = (input[1] - input[14]);
  int64_t s213 = (input[2] + input[13]);
  int64_t d213 = (input[2] - input[13]);
  int64_t s312 = (input[3] + input[12]);
  int64_t d312 = (input[3] - input[12]);
  int64_t s411 = (input[4] + input[11]);
  int64_t d411 = (input[4] - input[11]);
  int64_t s510 = (input[5] + input[10]);
  int64_t d510 = (input[5] - input[10]);
  int64_t s69  = (input[6] + input[9]);
  int64_t d69  = (input[6] - input[9]);
  int64_t s78  = (input[7] + input[8]);
  int64_t d78  = (input[7] - input[8]);

#if !CONFIG_EMULATE_HARDWARE
  (void)bd;
#endif

  sum = s015 * sinvalue_lookup[0] + s114 * sinvalue_lookup[1] +
        s213 * sinvalue_lookup[2] + s312 * sinvalue_lookup[3] +
        s411 * sinvalue_lookup[4] + s510 * sinvalue_lookup[5] +
        s69  * sinvalue_lookup[6] + s78  * sinvalue_lookup[7];
  output[0]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = d015 * sinvalue_lookup[1] + d114 * sinvalue_lookup[3] +
        d213 * sinvalue_lookup[5] + d312 * sinvalue_lookup[7] +
        d411 * sinvalue_lookup[6] + d510 * sinvalue_lookup[4] +
        d69  * sinvalue_lookup[2] + d78  * sinvalue_lookup[0];
  output[1]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = s015 * sinvalue_lookup[2] + s114 * sinvalue_lookup[5] +
        s213 * sinvalue_lookup[7] + s312 * sinvalue_lookup[4] +
        s411 * sinvalue_lookup[1] - s510 * sinvalue_lookup[0] -
        s69  * sinvalue_lookup[3] - s78  * sinvalue_lookup[6];
  output[2]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = d015 * sinvalue_lookup[3] + d114 * sinvalue_lookup[7] +
        d213 * sinvalue_lookup[4] + d312 * sinvalue_lookup[0] -
        d411 * sinvalue_lookup[2] - d510 * sinvalue_lookup[6] -
        d69  * sinvalue_lookup[5] - d78  * sinvalue_lookup[1];
  output[3]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = s015 * sinvalue_lookup[4] + s114 * sinvalue_lookup[6] +
        s213 * sinvalue_lookup[1] - s312 * sinvalue_lookup[2] -
        s411 * sinvalue_lookup[7] - s510 * sinvalue_lookup[3] +
        s69  * sinvalue_lookup[0] + s78  * sinvalue_lookup[5];
  output[4]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = d015 * sinvalue_lookup[5] + d114 * sinvalue_lookup[4] -
        d213 * sinvalue_lookup[0] - d312 * sinvalue_lookup[6] -
        d411 * sinvalue_lookup[3] + d510 * sinvalue_lookup[1] +
        d69  * sinvalue_lookup[7] + d78  * sinvalue_lookup[2];
  output[5]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = s015 * sinvalue_lookup[6] + s114 * sinvalue_lookup[2] -
        s213 * sinvalue_lookup[3] - s312 * sinvalue_lookup[5] +
        s411 * sinvalue_lookup[0] + s510 * sinvalue_lookup[7] +
        s69  * sinvalue_lookup[1] - s78  * sinvalue_lookup[4];
  output[6]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = d015 * sinvalue_lookup[7] + d114 * sinvalue_lookup[0] -
        d213 * sinvalue_lookup[6] - d312 * sinvalue_lookup[1] +
        d411 * sinvalue_lookup[5] + d510 * sinvalue_lookup[2] -
        d69  * sinvalue_lookup[4] - d78  * sinvalue_lookup[3];
  output[7]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = s015 * sinvalue_lookup[7] - s114 * sinvalue_lookup[0] -
        s213 * sinvalue_lookup[6] + s312 * sinvalue_lookup[1] +
        s411 * sinvalue_lookup[5] - s510 * sinvalue_lookup[2] -
        s69  * sinvalue_lookup[4] + s78  * sinvalue_lookup[3];
  output[8]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = d015 * sinvalue_lookup[6] - d114 * sinvalue_lookup[2] -
        d213 * sinvalue_lookup[3] + d312 * sinvalue_lookup[5] +
        d411 * sinvalue_lookup[0] - d510 * sinvalue_lookup[7] +
        d69  * sinvalue_lookup[1] + d78  * sinvalue_lookup[4];
  output[9]  = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = s015 * sinvalue_lookup[5] - s114 * sinvalue_lookup[4] -
        s213 * sinvalue_lookup[0] + s312 * sinvalue_lookup[6] -
        s411 * sinvalue_lookup[3] - s510 * sinvalue_lookup[1] +
        s69  * sinvalue_lookup[7] - s78  * sinvalue_lookup[2];
  output[10] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = d015 * sinvalue_lookup[4] - d114 * sinvalue_lookup[6] +
        d213 * sinvalue_lookup[1] + d312 * sinvalue_lookup[2] -
        d411 * sinvalue_lookup[7] + d510 * sinvalue_lookup[3] +
        d69  * sinvalue_lookup[0] - d78  * sinvalue_lookup[5];
  output[11] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = s015 * sinvalue_lookup[3] - s114 * sinvalue_lookup[7] +
        s213 * sinvalue_lookup[4] - s312 * sinvalue_lookup[0] -
        s411 * sinvalue_lookup[2] + s510 * sinvalue_lookup[6] -
        s69  * sinvalue_lookup[5] + s78  * sinvalue_lookup[1];
  output[12] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = d015 * sinvalue_lookup[2] - d114 * sinvalue_lookup[5] +
        d213 * sinvalue_lookup[7] - d312 * sinvalue_lookup[4] +
        d411 * sinvalue_lookup[1] + d510 * sinvalue_lookup[0] -
        d69  * sinvalue_lookup[3] + d78  * sinvalue_lookup[6];
  output[13] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = s015 * sinvalue_lookup[1] - s114 * sinvalue_lookup[3] +
        s213 * sinvalue_lookup[5] - s312 * sinvalue_lookup[7] +
        s411 * sinvalue_lookup[6] - s510 * sinvalue_lookup[4] +
        s69  * sinvalue_lookup[2] - s78  * sinvalue_lookup[0];
  output[14] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
  sum = d015 * sinvalue_lookup[0] - d114 * sinvalue_lookup[1] +
        d213 * sinvalue_lookup[2] - d312 * sinvalue_lookup[3] +
        d411 * sinvalue_lookup[4] - d510 * sinvalue_lookup[5] +
        d69  * sinvalue_lookup[6] - d78  * sinvalue_lookup[7];
  output[15] = WRAPLOW(ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS)), bd);
}

static void fliplr16(uint16_t *dest, int stride, int l) {
  int i, j;
  for (i = 0; i < l; ++i) {
    for (j = 0; j < l / 2; ++j) {
      const uint16_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[i * stride + l - 1 - j];
      dest[i * stride + l - 1 - j] = tmp;
    }
  }
}

static void flipud16(uint16_t *dest, int stride, int l) {
  int i, j;
  for (j = 0; j < l; ++j) {
    for (i = 0; i < l / 2; ++i) {
      const uint16_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[(l - 1 - i) * stride + j];
      dest[(l - 1 - i) * stride + j] = tmp;
    }
  }
}

static void fliplrud16(uint16_t *dest, int stride, int l) {
  int i, j;
  for (i = 0; i < l / 2; ++i) {
    for (j = 0; j < l; ++j) {
      const uint16_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[(l - 1 - i) * stride + l - 1 - j];
      dest[(l - 1 - i) * stride + l - 1 - j] = tmp;
    }
  }
}

static void highbd_inv_idtx_add_c(const tran_low_t *input, uint8_t *dest8,
                                  int stride, int bs, int bd) {
  int r, c;
  const int shift = bs < 32 ? 3 : 2;
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c)
      dest[c] = highbd_clip_pixel_add(dest[c], input[c] >> shift, bd);
    dest += stride;
    input += bs;
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_EXT_TX

void vp10_iht4x4_16_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                          int tx_type) {
  const transform_2d IHT_4[] = {
    { idct4_c, idct4_c  },   // DCT_DCT  = 0
    { iadst4_c, idct4_c  },  // ADST_DCT = 1
    { idct4_c, iadst4_c },   // DCT_ADST = 2
    { iadst4_c, iadst4_c },  // ADST_ADST = 3
#if CONFIG_EXT_TX
    { iadst4_c, idct4_c },   // FLIPADST_DCT = 4
    { idct4_c,  iadst4_c },  // DCT_FLIPADST = 5
    { iadst4_c, iadst4_c },  // FLIPADST_FLIPADST = 6
    { iadst4_c, iadst4_c },  // ADST_FLIPADST = 7
    { iadst4_c, iadst4_c },  // FLIPADST_ADST = 8
    { idst4_c,  idst4_c },   // DST_DST = 9
    { idst4_c,  idct4_c  },  // DST_DCT = 10
    { idct4_c,  idst4_c  },  // DCT_DST = 11
    { idst4_c,  iadst4_c },  // DST_ADST = 12
    { iadst4_c, idst4_c  },  // ADST_DST = 13
    { idst4_c,  iadst4_c },  // DST_FLIPADST = 14
    { iadst4_c, idst4_c  },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t out[4 * 4];
  tran_low_t *outptr = out;
  tran_low_t temp_in[4], temp_out[4];

  // inverse transform row vectors
  for (i = 0; i < 4; ++i) {
    IHT_4[tx_type].rows(input, outptr);
    input  += 4;
    outptr += 4;
  }

  // inverse transform column vectors
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      temp_in[j] = out[j * 4 + i];
    IHT_4[tx_type].cols(temp_in, temp_out);
    for (j = 0; j < 4; ++j) {
      dest[j * stride + i] = clip_pixel_add(dest[j * stride + i],
                                            ROUND_POWER_OF_TWO(temp_out[j], 4));
    }
  }
}

void vp10_iht8x8_64_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                         int tx_type) {
  static const transform_2d IHT_8[] = {
    { idct8_c,  idct8_c  },  // DCT_DCT  = 0
    { iadst8_c, idct8_c  },  // ADST_DCT = 1
    { idct8_c,  iadst8_c },  // DCT_ADST = 2
    { iadst8_c, iadst8_c },  // ADST_ADST = 3
#if CONFIG_EXT_TX
    { iadst8_c, idct8_c },   // FLIPADST_DCT = 4
    { idct8_c,  iadst8_c },  // DCT_FLIPADST = 5
    { iadst8_c, iadst8_c },  // FLIPADST_FLIPADST = 6
    { iadst8_c, iadst8_c },  // ADST_FLIPADST = 7
    { iadst8_c, iadst8_c },  // FLIPADST_ADST = 8
    { idst8_c,  idst8_c },   // DST_DST = 9
    { idst8_c,  idct8_c  },  // DST_DCT = 10
    { idct8_c,  idst8_c  },  // DCT_DST = 11
    { idst8_c,  iadst8_c },  // DST_ADST = 12
    { iadst8_c, idst8_c  },  // ADST_DST = 13
    { idst8_c,  iadst8_c },  // DST_FLIPADST = 14
    { iadst8_c, idst8_c  },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t out[8 * 8];
  tran_low_t *outptr = out;
  tran_low_t temp_in[8], temp_out[8];
  const transform_2d ht = IHT_8[tx_type];

  // inverse transform row vectors
  for (i = 0; i < 8; ++i) {
    ht.rows(input, outptr);
    input += 8;
    outptr += 8;
  }

  // inverse transform column vectors
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = out[j * 8 + i];
    ht.cols(temp_in, temp_out);
    for (j = 0; j < 8; ++j) {
      dest[j * stride + i] = clip_pixel_add(dest[j * stride + i],
                                            ROUND_POWER_OF_TWO(temp_out[j], 5));
    }
  }
}

void vp10_iht16x16_256_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                            int tx_type) {
  static const transform_2d IHT_16[] = {
    { idct16_c,  idct16_c  },  // DCT_DCT  = 0
    { iadst16_c, idct16_c  },  // ADST_DCT = 1
    { idct16_c,  iadst16_c },  // DCT_ADST = 2
    { iadst16_c, iadst16_c },  // ADST_ADST = 3
#if CONFIG_EXT_TX
    { iadst16_c, idct16_c  },  // FLIPADST_DCT = 4
    { idct16_c,  iadst16_c },  // DCT_FLIPADST = 5
    { iadst16_c, iadst16_c },  // FLIPADST_FLIPADST = 6
    { iadst16_c, iadst16_c },  // ADST_FLIPADST = 7
    { iadst16_c, iadst16_c },  // FLIPADST_ADST = 8
    { idst16_c,  idst16_c  },  // DST_DST = 9
    { idst16_c,  idct16_c  },  // DST_DCT = 10
    { idct16_c,  idst16_c  },  // DCT_DST = 11
    { idst16_c,  iadst16_c },  // DST_ADST = 12
    { iadst16_c, idst16_c  },  // ADST_DST = 13
    { idst16_c,  iadst16_c },  // DST_FLIPADST = 14
    { iadst16_c, idst16_c  },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t out[16 * 16];
  tran_low_t *outptr = out;
  tran_low_t temp_in[16], temp_out[16];
  const transform_2d ht = IHT_16[tx_type];

  // Rows
  for (i = 0; i < 16; ++i) {
    ht.rows(input, outptr);
    input += 16;
    outptr += 16;
  }

  // Columns
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j)
      temp_in[j] = out[j * 16 + i];
    ht.cols(temp_in, temp_out);
    for (j = 0; j < 16; ++j) {
      dest[j * stride + i] = clip_pixel_add(dest[j * stride + i],
                                            ROUND_POWER_OF_TWO(temp_out[j], 6));
    }
  }
}

// idct
void vp10_idct4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob) {
  if (eob > 1)
    vpx_idct4x4_16_add(input, dest, stride);
  else
    vpx_idct4x4_1_add(input, dest, stride);
}


void vp10_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob) {
  if (eob > 1)
    vpx_iwht4x4_16_add(input, dest, stride);
  else
    vpx_iwht4x4_1_add(input, dest, stride);
}

void vp10_idct8x8_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob) {
  // If dc is 1, then input[0] is the reconstructed value, do not need
  // dequantization. Also, when dc is 1, dc is counted in eobs, namely eobs >=1.

  // The calculation can be simplified if there are not many non-zero dct
  // coefficients. Use eobs to decide what to do.
  // TODO(yunqingwang): "eobs = 1" case is also handled in vp10_short_idct8x8_c.
  // Combine that with code here.
  if (eob == 1)
    // DC only DCT coefficient
    vpx_idct8x8_1_add(input, dest, stride);
  else if (eob <= 12)
    vpx_idct8x8_12_add(input, dest, stride);
  else
    vpx_idct8x8_64_add(input, dest, stride);
}

void vp10_idct16x16_add(const tran_low_t *input, uint8_t *dest, int stride,
                       int eob) {
  /* The calculation can be simplified if there are not many non-zero dct
   * coefficients. Use eobs to separate different cases. */
  if (eob == 1)
    /* DC only DCT coefficient. */
    vpx_idct16x16_1_add(input, dest, stride);
  else if (eob <= 10)
    vpx_idct16x16_10_add(input, dest, stride);
  else
    vpx_idct16x16_256_add(input, dest, stride);
}

void vp10_idct32x32_add(const tran_low_t *input, uint8_t *dest, int stride,
                       int eob) {
  if (eob == 1)
    vpx_idct32x32_1_add(input, dest, stride);
  else if (eob <= 34)
    // non-zero coeff only in upper-left 8x8
    vpx_idct32x32_34_add(input, dest, stride);
  else
    vpx_idct32x32_1024_add(input, dest, stride);
}

void vp10_inv_txfm_add_4x4(const tran_low_t *input, uint8_t *dest,
                           int stride, int eob, TX_TYPE tx_type, int lossless) {
  if (lossless) {
    assert(tx_type == DCT_DCT);
    vp10_iwht4x4_add(input, dest, stride, eob);
  } else {
    switch (tx_type) {
      case DCT_DCT:
        vp10_idct4x4_add(input, dest, stride, eob);
        break;
      case ADST_DCT:
      case DCT_ADST:
      case ADST_ADST:
        vp10_iht4x4_16_add(input, dest, stride, tx_type);
        break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      flipud(dest, stride, 4);
      vp10_iht4x4_16_add(input, dest, stride, ADST_DCT);
      flipud(dest, stride, 4);
        break;
    case DCT_FLIPADST:
      fliplr(dest, stride, 4);
      vp10_iht4x4_16_add(input, dest, stride, DCT_ADST);
      fliplr(dest, stride, 4);
      break;
    case FLIPADST_FLIPADST:
      fliplrud(dest, stride, 4);
      vp10_iht4x4_16_add(input, dest, stride, ADST_ADST);
      fliplrud(dest, stride, 4);
      break;
    case ADST_FLIPADST:
      fliplr(dest, stride, 4);
      vp10_iht4x4_16_add(input, dest, stride, ADST_ADST);
      fliplr(dest, stride, 4);
      break;
    case FLIPADST_ADST:
      flipud(dest, stride, 4);
      vp10_iht4x4_16_add(input, dest, stride, ADST_ADST);
      flipud(dest, stride, 4);
      break;
    case DST_DST:
    case DST_DCT:
    case DCT_DST:
    case DST_ADST:
    case ADST_DST:
      // Use C version since DST only exists in C code
      vp10_iht4x4_16_add_c(input, dest, stride, tx_type);
      break;
    case FLIPADST_DST:
      flipud(dest, stride, 4);
      vp10_iht4x4_16_add_c(input, dest, stride, ADST_DST);
      flipud(dest, stride, 4);
      break;
    case DST_FLIPADST:
      fliplr(dest, stride, 4);
      vp10_iht4x4_16_add_c(input, dest, stride, DST_ADST);
      fliplr(dest, stride, 4);
      break;
    case IDTX:
      inv_idtx_add_c(input, dest, stride, 4);
      break;
#endif  // CONFIG_EXT_TX
      default:
        assert(0);
        break;
    }
  }
}

void vp10_inv_txfm_add_8x8(const tran_low_t *input, uint8_t *dest,
                           int stride, int eob, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
      vp10_idct8x8_add(input, dest, stride, eob);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_iht8x8_64_add(input, dest, stride, tx_type);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      flipud(dest, stride, 8);
      vp10_iht8x8_64_add(input, dest, stride, ADST_DCT);
      flipud(dest, stride, 8);
      break;
    case DCT_FLIPADST:
      fliplr(dest, stride, 8);
      vp10_iht8x8_64_add(input, dest, stride, DCT_ADST);
      fliplr(dest, stride, 8);
      break;
    case FLIPADST_FLIPADST:
      fliplrud(dest, stride, 8);
      vp10_iht8x8_64_add(input, dest, stride, ADST_ADST);
      fliplrud(dest, stride, 8);
      break;
    case ADST_FLIPADST:
      fliplr(dest, stride, 8);
      vp10_iht8x8_64_add(input, dest, stride, ADST_ADST);
      fliplr(dest, stride, 8);
      break;
    case FLIPADST_ADST:
      flipud(dest, stride, 8);
      vp10_iht8x8_64_add(input, dest, stride, ADST_ADST);
      flipud(dest, stride, 8);
      break;
    case DST_DST:
    case DST_DCT:
    case DCT_DST:
    case DST_ADST:
    case ADST_DST:
      // Use C version since DST only exists in C code
      vp10_iht8x8_64_add_c(input, dest, stride, tx_type);
      break;
    case FLIPADST_DST:
      flipud(dest, stride, 8);
      vp10_iht8x8_64_add_c(input, dest, stride, ADST_DST);
      flipud(dest, stride, 8);
      break;
    case DST_FLIPADST:
      fliplr(dest, stride, 8);
      vp10_iht8x8_64_add_c(input, dest, stride, DST_ADST);
      fliplr(dest, stride, 8);
      break;
    case IDTX:
      inv_idtx_add_c(input, dest, stride, 8);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}

void vp10_inv_txfm_add_16x16(const tran_low_t *input, uint8_t *dest,
                             int stride, int eob, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
      vp10_idct16x16_add(input, dest, stride, eob);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_iht16x16_256_add(input, dest, stride, tx_type);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      flipud(dest, stride, 16);
      vp10_iht16x16_256_add(input, dest, stride, ADST_DCT);
      flipud(dest, stride, 16);
      break;
    case DCT_FLIPADST:
      fliplr(dest, stride, 16);
      vp10_iht16x16_256_add(input, dest, stride, DCT_ADST);
      fliplr(dest, stride, 16);
      break;
    case FLIPADST_FLIPADST:
      fliplrud(dest, stride, 16);
      vp10_iht16x16_256_add(input, dest, stride, ADST_ADST);
      fliplrud(dest, stride, 16);
      break;
    case ADST_FLIPADST:
      fliplr(dest, stride, 16);
      vp10_iht16x16_256_add(input, dest, stride, ADST_ADST);
      fliplr(dest, stride, 16);
      break;
    case FLIPADST_ADST:
      flipud(dest, stride, 16);
      vp10_iht16x16_256_add(input, dest, stride, ADST_ADST);
      flipud(dest, stride, 16);
      break;
    case DST_DST:
    case DST_DCT:
    case DCT_DST:
    case DST_ADST:
    case ADST_DST:
      // Use C version since DST only exists in C code
      vp10_iht16x16_256_add_c(input, dest, stride, tx_type);
      break;
    case FLIPADST_DST:
      flipud(dest, stride, 16);
      vp10_iht16x16_256_add_c(input, dest, stride, ADST_DST);
      flipud(dest, stride, 16);
      break;
    case DST_FLIPADST:
      fliplr(dest, stride, 16);
      vp10_iht16x16_256_add_c(input, dest, stride, DST_ADST);
      fliplr(dest, stride, 16);
      break;
    case IDTX:
      inv_idtx_add_c(input, dest, stride, 16);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}

void vp10_inv_txfm_add_32x32(const tran_low_t *input, uint8_t *dest,
                             int stride, int eob, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
      vp10_idct32x32_add(input, dest, stride, eob);
      break;
#if CONFIG_EXT_TX
    case IDTX:
      inv_idtx_add_c(input, dest, stride, 32);
      break;
#endif  // CONFIG_EXT_TX
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      assert(0);
      break;
    default:
      assert(0);
      break;
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
void vp10_highbd_iht4x4_16_add_c(const tran_low_t *input, uint8_t *dest8,
                                int stride, int tx_type, int bd) {
  const highbd_transform_2d IHT_4[] = {
    { vpx_highbd_idct4_c,  vpx_highbd_idct4_c  },  // DCT_DCT  = 0
    { vpx_highbd_iadst4_c, vpx_highbd_idct4_c  },  // ADST_DCT = 1
    { vpx_highbd_idct4_c,  vpx_highbd_iadst4_c },  // DCT_ADST = 2
    { vpx_highbd_iadst4_c, vpx_highbd_iadst4_c },  // ADST_ADST = 3
#if CONFIG_EXT_TX
    { vpx_highbd_iadst4_c, vpx_highbd_idct4_c  },  // FLIPADST_DCT = 4
    { vpx_highbd_idct4_c,  vpx_highbd_iadst4_c },  // DCT_FLIPADST = 5
    { vpx_highbd_iadst4_c, vpx_highbd_iadst4_c },  // FLIPADST_FLIPADST = 6
    { vpx_highbd_iadst4_c, vpx_highbd_iadst4_c },  // ADST_FLIPADST = 7
    { vpx_highbd_iadst4_c, vpx_highbd_iadst4_c },  // FLIPADST_ADST = 8
    { highbd_idst4_c,      highbd_idst4_c      },  // DST_DST = 9
    { highbd_idst4_c,      vpx_highbd_idct4_c  },  // DST_DCT = 10
    { vpx_highbd_idct4_c,  highbd_idst4_c      },  // DCT_DST = 11
    { highbd_idst4_c,      vpx_highbd_iadst4_c },  // DST_ADST = 12
    { vpx_highbd_iadst4_c, highbd_idst4_c      },  // ADST_DST = 13
    { highbd_idst4_c,      vpx_highbd_iadst4_c },  // DST_FLIPADST = 14
    { vpx_highbd_iadst4_c, highbd_idst4_c      },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
  };
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  int i, j;
  tran_low_t out[4 * 4];
  tran_low_t *outptr = out;
  tran_low_t temp_in[4], temp_out[4];

  // Inverse transform row vectors.
  for (i = 0; i < 4; ++i) {
    IHT_4[tx_type].rows(input, outptr, bd);
    input  += 4;
    outptr += 4;
  }

  // Inverse transform column vectors.
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      temp_in[j] = out[j * 4 + i];
    IHT_4[tx_type].cols(temp_in, temp_out, bd);
    for (j = 0; j < 4; ++j) {
      dest[j * stride + i] = highbd_clip_pixel_add(
          dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 4), bd);
    }
  }
}

void vp10_highbd_iht8x8_64_add_c(const tran_low_t *input, uint8_t *dest8,
                                int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_8[] = {
    { vpx_highbd_idct8_c,  vpx_highbd_idct8_c  },  // DCT_DCT  = 0
    { vpx_highbd_iadst8_c, vpx_highbd_idct8_c  },  // ADST_DCT = 1
    { vpx_highbd_idct8_c,  vpx_highbd_iadst8_c },  // DCT_ADST = 2
    { vpx_highbd_iadst8_c, vpx_highbd_iadst8_c },  // ADST_ADST = 3
#if CONFIG_EXT_TX
    { vpx_highbd_iadst8_c, vpx_highbd_idct8_c  },  // FLIPADST_DCT = 4
    { vpx_highbd_idct8_c,  vpx_highbd_iadst8_c },  // DCT_FLIPADST = 5
    { vpx_highbd_iadst8_c, vpx_highbd_iadst8_c },  // FLIPADST_FLIPADST = 6
    { vpx_highbd_iadst8_c, vpx_highbd_iadst8_c },  // ADST_FLIPADST = 7
    { vpx_highbd_iadst8_c, vpx_highbd_iadst8_c },  // FLIPADST_ADST = 8
    { highbd_idst8_c,      highbd_idst8_c      },  // DST_DST = 9
    { highbd_idst8_c,      vpx_highbd_idct8_c  },  // DST_DCT = 10
    { vpx_highbd_idct8_c,  highbd_idst8_c      },  // DCT_DST = 11
    { highbd_idst8_c,      vpx_highbd_iadst8_c },  // DST_ADST = 12
    { vpx_highbd_iadst8_c, highbd_idst8_c      },  // ADST_DST = 13
    { highbd_idst8_c,      vpx_highbd_iadst8_c },  // DST_FLIPADST = 14
    { vpx_highbd_iadst8_c, highbd_idst8_c      },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t out[8 * 8];
  tran_low_t *outptr = out;
  tran_low_t temp_in[8], temp_out[8];
  const highbd_transform_2d ht = HIGH_IHT_8[tx_type];
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  // Inverse transform row vectors.
  for (i = 0; i < 8; ++i) {
    ht.rows(input, outptr, bd);
    input += 8;
    outptr += 8;
  }

  // Inverse transform column vectors.
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = out[j * 8 + i];
    ht.cols(temp_in, temp_out, bd);
    for (j = 0; j < 8; ++j) {
      dest[j * stride + i] = highbd_clip_pixel_add(
          dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 5), bd);
    }
  }
}

void vp10_highbd_iht16x16_256_add_c(const tran_low_t *input, uint8_t *dest8,
                                   int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_16[] = {
    { vpx_highbd_idct16_c,  vpx_highbd_idct16_c  },  // DCT_DCT  = 0
    { vpx_highbd_iadst16_c, vpx_highbd_idct16_c  },  // ADST_DCT = 1
    { vpx_highbd_idct16_c,  vpx_highbd_iadst16_c },  // DCT_ADST = 2
    { vpx_highbd_iadst16_c, vpx_highbd_iadst16_c },  // ADST_ADST = 3
#if CONFIG_EXT_TX
    { vpx_highbd_iadst16_c, vpx_highbd_idct16_c  },  // FLIPADST_DCT = 4
    { vpx_highbd_idct16_c,  vpx_highbd_iadst16_c },  // DCT_FLIPADST = 5
    { vpx_highbd_iadst16_c, vpx_highbd_iadst16_c },  // FLIPADST_FLIPADST = 6
    { vpx_highbd_iadst16_c, vpx_highbd_iadst16_c },  // ADST_FLIPADST = 7
    { vpx_highbd_iadst16_c, vpx_highbd_iadst16_c },  // FLIPADST_ADST = 8
    { highbd_idst16_c,      highbd_idst16_c      },  // DST_DST = 9
    { highbd_idst16_c,      vpx_highbd_idct16_c  },  // DST_DCT = 10
    { vpx_highbd_idct16_c,  highbd_idst16_c      },  // DCT_DST = 11
    { highbd_idst16_c,      vpx_highbd_iadst16_c },  // DST_ADST = 12
    { vpx_highbd_iadst16_c, highbd_idst16_c      },  // ADST_DST = 13
    { highbd_idst16_c,      vpx_highbd_iadst16_c },  // DST_FLIPADST = 14
    { vpx_highbd_iadst16_c, highbd_idst16_c      },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t out[16 * 16];
  tran_low_t *outptr = out;
  tran_low_t temp_in[16], temp_out[16];
  const highbd_transform_2d ht = HIGH_IHT_16[tx_type];
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  // Rows
  for (i = 0; i < 16; ++i) {
    ht.rows(input, outptr, bd);
    input += 16;
    outptr += 16;
  }

  // Columns
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j)
      temp_in[j] = out[j * 16 + i];
    ht.cols(temp_in, temp_out, bd);
    for (j = 0; j < 16; ++j) {
      dest[j * stride + i] = highbd_clip_pixel_add(
          dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 6), bd);
    }
  }
}

// idct
void vp10_highbd_idct4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd) {
  if (eob > 1)
    vpx_highbd_idct4x4_16_add(input, dest, stride, bd);
  else
    vpx_highbd_idct4x4_1_add(input, dest, stride, bd);
}


void vp10_highbd_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd) {
  if (eob > 1)
    vpx_highbd_iwht4x4_16_add(input, dest, stride, bd);
  else
    vpx_highbd_iwht4x4_1_add(input, dest, stride, bd);
}

void vp10_highbd_idct8x8_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd) {
  // If dc is 1, then input[0] is the reconstructed value, do not need
  // dequantization. Also, when dc is 1, dc is counted in eobs, namely eobs >=1.

  // The calculation can be simplified if there are not many non-zero dct
  // coefficients. Use eobs to decide what to do.
  // TODO(yunqingwang): "eobs = 1" case is also handled in vp10_short_idct8x8_c.
  // Combine that with code here.
  // DC only DCT coefficient
  if (eob == 1) {
    vpx_highbd_idct8x8_1_add(input, dest, stride, bd);
  } else if (eob <= 10) {
    vpx_highbd_idct8x8_10_add(input, dest, stride, bd);
  } else {
    vpx_highbd_idct8x8_64_add(input, dest, stride, bd);
  }
}

void vp10_highbd_idct16x16_add(const tran_low_t *input, uint8_t *dest,
                              int stride, int eob, int bd) {
  // The calculation can be simplified if there are not many non-zero dct
  // coefficients. Use eobs to separate different cases.
  // DC only DCT coefficient.
  if (eob == 1) {
    vpx_highbd_idct16x16_1_add(input, dest, stride, bd);
  } else if (eob <= 10) {
    vpx_highbd_idct16x16_10_add(input, dest, stride, bd);
  } else {
    vpx_highbd_idct16x16_256_add(input, dest, stride, bd);
  }
}

void vp10_highbd_idct32x32_add(const tran_low_t *input, uint8_t *dest,
                              int stride, int eob, int bd) {
  // Non-zero coeff only in upper-left 8x8
  if (eob == 1) {
    vpx_highbd_idct32x32_1_add(input, dest, stride, bd);
  } else if (eob <= 34) {
    vpx_highbd_idct32x32_34_add(input, dest, stride, bd);
  } else {
    vpx_highbd_idct32x32_1024_add(input, dest, stride, bd);
  }
}

void vp10_highbd_inv_txfm_add_4x4(const tran_low_t *input, uint8_t *dest,
                                  int stride, int eob, int bd, TX_TYPE tx_type,
                                  int lossless) {
  if (lossless) {
    assert(tx_type == DCT_DCT);
    vp10_highbd_iwht4x4_add(input, dest, stride, eob, bd);
  } else {
    switch (tx_type) {
      case DCT_DCT:
        vp10_highbd_idct4x4_add(input, dest, stride, eob, bd);
        break;
      case ADST_DCT:
      case DCT_ADST:
      case ADST_ADST:
         vp10_highbd_iht4x4_16_add(input, dest, stride, tx_type, bd);
         break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add(input, dest, stride, ADST_DCT, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
         break;
    case DCT_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add(input, dest, stride, DCT_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      break;
    case FLIPADST_FLIPADST:
      fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add(input, dest, stride, ADST_ADST, bd);
      fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      break;
    case ADST_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add(input, dest, stride, ADST_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      break;
    case FLIPADST_ADST:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add(input, dest, stride, ADST_ADST, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      break;
    case DST_DST:
    case DST_DCT:
    case DCT_DST:
    case DST_ADST:
    case ADST_DST:
      // Use C version since DST only exists in C code
      vp10_highbd_iht4x4_16_add_c(input, dest, stride, tx_type, bd);
      break;
    case FLIPADST_DST:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add_c(input, dest, stride, ADST_DST, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      break;
    case DST_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add_c(input, dest, stride, DST_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      break;
    case IDTX:
      highbd_inv_idtx_add_c(input, dest, stride, 4, bd);
      break;
#endif  // CONFIG_EXT_TX
      default:
        assert(0);
        break;
    }
  }
}

void vp10_highbd_inv_txfm_add_8x8(const tran_low_t *input, uint8_t *dest,
                                  int stride, int eob, int bd,
                                  TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
      vp10_highbd_idct8x8_add(input, dest, stride, eob, bd);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_highbd_iht8x8_64_add(input, dest, stride, tx_type, bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add(input, dest, stride, ADST_DCT, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      break;
    case DCT_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add(input, dest, stride, DCT_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      break;
    case FLIPADST_FLIPADST:
      fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add(input, dest, stride, ADST_ADST, bd);
      fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      break;
    case ADST_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add(input, dest, stride, ADST_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      break;
    case FLIPADST_ADST:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add(input, dest, stride, ADST_ADST, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      break;
    case DST_DST:
    case DST_DCT:
    case DCT_DST:
    case DST_ADST:
    case ADST_DST:
      // Use C version since DST only exists in C code
      vp10_highbd_iht8x8_64_add_c(input, dest, stride, tx_type, bd);
      break;
    case FLIPADST_DST:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add_c(input, dest, stride, ADST_DST, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      break;
    case DST_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add_c(input, dest, stride, DST_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      break;
    case IDTX:
      highbd_inv_idtx_add_c(input, dest, stride, 8, bd);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}

void vp10_highbd_inv_txfm_add_16x16(const tran_low_t *input, uint8_t *dest,
                                    int stride, int eob, int bd,
                                    TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
      vp10_highbd_idct16x16_add(input, dest, stride, eob, bd);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_highbd_iht16x16_256_add(input, dest, stride, tx_type, bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add(input, dest, stride, ADST_DCT, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      break;
    case DCT_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add(input, dest, stride, DCT_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      break;
    case FLIPADST_FLIPADST:
      fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add(input, dest, stride, ADST_ADST, bd);
      fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      break;
    case ADST_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add(input, dest, stride, ADST_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      break;
    case FLIPADST_ADST:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add(input, dest, stride, ADST_ADST, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      break;
    case DST_DST:
    case DST_DCT:
    case DCT_DST:
    case DST_ADST:
    case ADST_DST:
      // Use C version since DST only exists in C code
      vp10_highbd_iht16x16_256_add_c(input, dest, stride, tx_type, bd);
      break;
    case FLIPADST_DST:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add_c(input, dest, stride, ADST_DST, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      break;
    case DST_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add_c(input, dest, stride, DST_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      break;
    case IDTX:
      highbd_inv_idtx_add_c(input, dest, stride, 16, bd);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}

void vp10_highbd_inv_txfm_add_32x32(const tran_low_t *input, uint8_t *dest,
                                    int stride, int eob, int bd,
                                    TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
      vp10_highbd_idct32x32_add(input, dest, stride, eob, bd);
      break;
#if CONFIG_EXT_TX
    case IDTX:
      highbd_inv_idtx_add_c(input, dest, stride, 32, bd);
      break;
#endif  // CONFIG_EXT_TX
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      assert(0);
      break;
    default:
      assert(0);
      break;
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
