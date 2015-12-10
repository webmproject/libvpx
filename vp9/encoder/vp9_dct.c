/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <math.h>

#include "./vpx_config.h"
#include "./vp9_rtcd.h"

#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_idct.h"
#include "vp9/common/vp9_systemdependent.h"
#include "vp9/encoder/vp9_dct.h"

static INLINE tran_high_t fdct_round_shift(tran_high_t input) {
  tran_high_t rv = ROUND_POWER_OF_TWO(input, DCT_CONST_BITS);
  // TODO(debargha, peter.derivaz): Find new bounds for this assert
  // and make the bounds consts.
  // assert(INT16_MIN <= rv && rv <= INT16_MAX);
  return rv;
}

#if CONFIG_EXT_TX
void vp9_fklt4(const tran_low_t *input, tran_low_t *output) {
  vp9_fgentx4(input, output, Tx4);
}

void vp9_fklt8(const tran_low_t *input, tran_low_t *output) {
  vp9_fgentx8(input, output, Tx8);
}

void vp9_fklt16(const tran_low_t *input, tran_low_t *output) {
  vp9_fgentx16(input, output, Tx16);
}

void vp9_fdst4(const tran_low_t *input, tran_low_t *output) {
#if USE_DST2
  // vp9_fgentx4(input, output, Tx4);
  tran_high_t step[4];
  tran_high_t temp1, temp2;

  step[0] = input[0] - input[3];
  step[1] = -input[1] + input[2];
  step[2] = -input[1] - input[2];
  step[3] = input[0] + input[3];

  temp1 = (step[0] + step[1]) * cospi_16_64;
  temp2 = (step[0] - step[1]) * cospi_16_64;
  output[3] = fdct_round_shift(temp1);
  output[1] = fdct_round_shift(temp2);
  temp1 = step[2] * cospi_24_64 + step[3] * cospi_8_64;
  temp2 = -step[2] * cospi_8_64 + step[3] * cospi_24_64;
  output[2] = fdct_round_shift(temp1);
  output[0] = fdct_round_shift(temp2);
#else
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
  output[0] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = d03 * sinvalue_lookup[1] + d12 * sinvalue_lookup[0];
  output[1] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = s03 * sinvalue_lookup[1] - s12 * sinvalue_lookup[0];
  output[2] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = d03 * sinvalue_lookup[0] - d12 * sinvalue_lookup[1];
  output[3] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
#endif
}

void vp9_fdst8(const tran_low_t *input, tran_low_t *output) {
#if USE_DST2
  // vp9_fgentx8(input, output, Tx8);
  tran_high_t s0, s1, s2, s3, s4, s5, s6, s7;  // canbe16
  tran_high_t t0, t1, t2, t3;                  // needs32
  tran_high_t x0, x1, x2, x3;                  // canbe16

  // stage 1
  s0 = input[0] - input[7];
  s1 = -input[1] + input[6];
  s2 = input[2] - input[5];
  s3 = -input[3] + input[4];
  s4 = -input[3] - input[4];
  s5 = input[2] + input[5];
  s6 = -input[1] - input[6];
  s7 = input[0] + input[7];

  x0 = s0 + s3;
  x1 = s1 + s2;
  x2 = s1 - s2;
  x3 = s0 - s3;
  t0 = (x0 + x1) * cospi_16_64;
  t1 = (x0 - x1) * cospi_16_64;
  t2 =  x2 * cospi_24_64 + x3 *  cospi_8_64;
  t3 = -x2 * cospi_8_64  + x3 * cospi_24_64;
  output[7] = fdct_round_shift(t0);
  output[5] = fdct_round_shift(t2);
  output[3] = fdct_round_shift(t1);
  output[1] = fdct_round_shift(t3);

  // Stage 2
  t0 = (s6 - s5) * cospi_16_64;
  t1 = (s6 + s5) * cospi_16_64;
  t2 = fdct_round_shift(t0);
  t3 = fdct_round_shift(t1);

  // Stage 3
  x0 = s4 + t2;
  x1 = s4 - t2;
  x2 = s7 - t3;
  x3 = s7 + t3;

  // Stage 4
  t0 = x0 * cospi_28_64 + x3 *   cospi_4_64;
  t1 = x1 * cospi_12_64 + x2 *  cospi_20_64;
  t2 = x2 * cospi_12_64 + x1 * -cospi_20_64;
  t3 = x3 * cospi_28_64 + x0 *  -cospi_4_64;
  output[6] = fdct_round_shift(t0);
  output[4] = fdct_round_shift(t2);
  output[2] = fdct_round_shift(t1);
  output[0] = fdct_round_shift(t3);
#else
  // {sin(pi/9), sin(pi*2/9), ..., sin(pi*4/9)} * sqrt(2/9) * 2
  static const int sinvalue_lookup[] = {
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
  output[0] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = d07 * sinvalue_lookup[1] + d16 * sinvalue_lookup[3] +
        d25 * sinvalue_lookup[2] + d34 * sinvalue_lookup[0];
  output[1] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = (s07 + s16 - s34)* sinvalue_lookup[2];
  output[2] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = d07 * sinvalue_lookup[3] + d16 * sinvalue_lookup[0] -
        d25 * sinvalue_lookup[2] - d34 * sinvalue_lookup[1];
  output[3] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = s07 * sinvalue_lookup[3] - s16 * sinvalue_lookup[0] -
        s25 * sinvalue_lookup[2] + s34 * sinvalue_lookup[1];
  output[4] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = (d07 - d16 + d34)* sinvalue_lookup[2];
  output[5] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = s07 * sinvalue_lookup[1] - s16 * sinvalue_lookup[3] +
        s25 * sinvalue_lookup[2] - s34 * sinvalue_lookup[0];
  output[6] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = d07 * sinvalue_lookup[0] - d16 * sinvalue_lookup[1] +
        d25 * sinvalue_lookup[2] - d34 * sinvalue_lookup[3];
  output[7] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
#endif
}

void vp9_fdst16(const tran_low_t *input, tran_low_t *output) {
#if USE_DST2
  // vp9_fgentx16(input, output, Tx16);
  tran_high_t step1[8];      // canbe16
  tran_high_t step2[8];      // canbe16
  tran_high_t step3[8];      // canbe16
  tran_high_t in[8];      // canbe16
  tran_high_t temp1, temp2;  // needs32

  // step 1
  in[0] = input[0] - input[15];
  in[1] = -input[1] + input[14];
  in[2] = input[2] - input[13];
  in[3] = -input[3] + input[12];
  in[4] = input[4] - input[11];
  in[5] = -input[5] + input[10];
  in[6] = input[6] - input[ 9];
  in[7] = -input[7] + input[ 8];

  step1[0] = -input[7] - input[ 8];
  step1[1] = input[6] + input[ 9];
  step1[2] = -input[5] - input[10];
  step1[3] = input[4] + input[11];
  step1[4] = -input[3] - input[12];
  step1[5] = input[2] + input[13];
  step1[6] = -input[1] - input[14];
  step1[7] = input[0] + input[15];

  // fdct8(step, step);
  {
    tran_high_t s0, s1, s2, s3, s4, s5, s6, s7;  // canbe16
    tran_high_t t0, t1, t2, t3;                  // needs32
    tran_high_t x0, x1, x2, x3;                  // canbe16

    // stage 1
    s0 = in[0] + in[7];
    s1 = in[1] + in[6];
    s2 = in[2] + in[5];
    s3 = in[3] + in[4];
    s4 = in[3] - in[4];
    s5 = in[2] - in[5];
    s6 = in[1] - in[6];
    s7 = in[0] - in[7];

    // fdct4(step, step);
    x0 = s0 + s3;
    x1 = s1 + s2;
    x2 = s1 - s2;
    x3 = s0 - s3;
    t0 = (x0 + x1) * cospi_16_64;
    t1 = (x0 - x1) * cospi_16_64;
    t2 = x3 * cospi_8_64  + x2 * cospi_24_64;
    t3 = x3 * cospi_24_64 - x2 * cospi_8_64;
    output[15] = fdct_round_shift(t0);
    output[11] = fdct_round_shift(t2);
    output[7] = fdct_round_shift(t1);
    output[3] = fdct_round_shift(t3);

    // Stage 2
    t0 = (s6 - s5) * cospi_16_64;
    t1 = (s6 + s5) * cospi_16_64;
    t2 = fdct_round_shift(t0);
    t3 = fdct_round_shift(t1);

    // Stage 3
    x0 = s4 + t2;
    x1 = s4 - t2;
    x2 = s7 - t3;
    x3 = s7 + t3;

    // Stage 4
    t0 = x0 * cospi_28_64 + x3 *   cospi_4_64;
    t1 = x1 * cospi_12_64 + x2 *  cospi_20_64;
    t2 = x2 * cospi_12_64 + x1 * -cospi_20_64;
    t3 = x3 * cospi_28_64 + x0 *  -cospi_4_64;
    output[13] = fdct_round_shift(t0);
    output[9] = fdct_round_shift(t2);
    output[5] = fdct_round_shift(t1);
    output[1] = fdct_round_shift(t3);
  }

  // step 2
  temp1 = (step1[5] - step1[2]) * cospi_16_64;
  temp2 = (step1[4] - step1[3]) * cospi_16_64;
  step2[2] = fdct_round_shift(temp1);
  step2[3] = fdct_round_shift(temp2);
  temp1 = (step1[4] + step1[3]) * cospi_16_64;
  temp2 = (step1[5] + step1[2]) * cospi_16_64;
  step2[4] = fdct_round_shift(temp1);
  step2[5] = fdct_round_shift(temp2);

  // step 3
  step3[0] = step1[0] + step2[3];
  step3[1] = step1[1] + step2[2];
  step3[2] = step1[1] - step2[2];
  step3[3] = step1[0] - step2[3];
  step3[4] = step1[7] - step2[4];
  step3[5] = step1[6] - step2[5];
  step3[6] = step1[6] + step2[5];
  step3[7] = step1[7] + step2[4];

  // step 4
  temp1 = step3[1] *  -cospi_8_64 + step3[6] * cospi_24_64;
  temp2 = step3[2] * cospi_24_64 + step3[5] *  cospi_8_64;
  step2[1] = fdct_round_shift(temp1);
  step2[2] = fdct_round_shift(temp2);
  temp1 = step3[2] * cospi_8_64 - step3[5] * cospi_24_64;
  temp2 = step3[1] * cospi_24_64 + step3[6] *  cospi_8_64;
  step2[5] = fdct_round_shift(temp1);
  step2[6] = fdct_round_shift(temp2);

  // step 5
  step1[0] = step3[0] + step2[1];
  step1[1] = step3[0] - step2[1];
  step1[2] = step3[3] + step2[2];
  step1[3] = step3[3] - step2[2];
  step1[4] = step3[4] - step2[5];
  step1[5] = step3[4] + step2[5];
  step1[6] = step3[7] - step2[6];
  step1[7] = step3[7] + step2[6];

  // step 6
  temp1 = step1[0] * cospi_30_64 + step1[7] *  cospi_2_64;
  temp2 = step1[1] * cospi_14_64 + step1[6] * cospi_18_64;
  output[14] = fdct_round_shift(temp1);
  output[6] = fdct_round_shift(temp2);

  temp1 = step1[2] * cospi_22_64 + step1[5] * cospi_10_64;
  temp2 = step1[3] *  cospi_6_64 + step1[4] * cospi_26_64;
  output[10] = fdct_round_shift(temp1);
  output[2] = fdct_round_shift(temp2);

  temp1 = step1[3] * -cospi_26_64 + step1[4] *  cospi_6_64;
  temp2 = step1[2] * -cospi_10_64 + step1[5] * cospi_22_64;
  output[12] = fdct_round_shift(temp1);
  output[4] = fdct_round_shift(temp2);

  temp1 = step1[1] * -cospi_18_64 + step1[6] * cospi_14_64;
  temp2 = step1[0] *  -cospi_2_64 + step1[7] * cospi_30_64;
  output[8] = fdct_round_shift(temp1);
  output[0] = fdct_round_shift(temp2);
#else
  // {sin(pi/17), sin(pi*2/17, ..., sin(pi*8/17)} * sqrt(2/17) * 2 * sqrt(2)
  static const int sinvalue_lookup[] = {
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
  output[0]  = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = d015 * sinvalue_lookup[1] + d114 * sinvalue_lookup[3] +
        d213 * sinvalue_lookup[5] + d312 * sinvalue_lookup[7] +
        d411 * sinvalue_lookup[6] + d510 * sinvalue_lookup[4] +
        d69  * sinvalue_lookup[2] + d78  * sinvalue_lookup[0];
  output[1]  = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = s015 * sinvalue_lookup[2] + s114 * sinvalue_lookup[5] +
        s213 * sinvalue_lookup[7] + s312 * sinvalue_lookup[4] +
        s411 * sinvalue_lookup[1] - s510 * sinvalue_lookup[0] -
        s69  * sinvalue_lookup[3] - s78  * sinvalue_lookup[6];
  output[2]  = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = d015 * sinvalue_lookup[3] + d114 * sinvalue_lookup[7] +
        d213 * sinvalue_lookup[4] + d312 * sinvalue_lookup[0] -
        d411 * sinvalue_lookup[2] - d510 * sinvalue_lookup[6] -
        d69  * sinvalue_lookup[5] - d78  * sinvalue_lookup[1];
  output[3]  = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = s015 * sinvalue_lookup[4] + s114 * sinvalue_lookup[6] +
        s213 * sinvalue_lookup[1] - s312 * sinvalue_lookup[2] -
        s411 * sinvalue_lookup[7] - s510 * sinvalue_lookup[3] +
        s69  * sinvalue_lookup[0] + s78  * sinvalue_lookup[5];
  output[4]  = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = d015 * sinvalue_lookup[5] + d114 * sinvalue_lookup[4] -
        d213 * sinvalue_lookup[0] - d312 * sinvalue_lookup[6] -
        d411 * sinvalue_lookup[3] + d510 * sinvalue_lookup[1] +
        d69  * sinvalue_lookup[7] + d78  * sinvalue_lookup[2];
  output[5]  = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = s015 * sinvalue_lookup[6] + s114 * sinvalue_lookup[2] -
        s213 * sinvalue_lookup[3] - s312 * sinvalue_lookup[5] +
        s411 * sinvalue_lookup[0] + s510 * sinvalue_lookup[7] +
        s69  * sinvalue_lookup[1] - s78  * sinvalue_lookup[4];
  output[6]  = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = d015 * sinvalue_lookup[7] + d114 * sinvalue_lookup[0] -
        d213 * sinvalue_lookup[6] - d312 * sinvalue_lookup[1] +
        d411 * sinvalue_lookup[5] + d510 * sinvalue_lookup[2] -
        d69  * sinvalue_lookup[4] - d78  * sinvalue_lookup[3];
  output[7]  = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = s015 * sinvalue_lookup[7] - s114 * sinvalue_lookup[0] -
        s213 * sinvalue_lookup[6] + s312 * sinvalue_lookup[1] +
        s411 * sinvalue_lookup[5] - s510 * sinvalue_lookup[2] -
        s69  * sinvalue_lookup[4] + s78  * sinvalue_lookup[3];
  output[8]  = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = d015 * sinvalue_lookup[6] - d114 * sinvalue_lookup[2] -
        d213 * sinvalue_lookup[3] + d312 * sinvalue_lookup[5] +
        d411 * sinvalue_lookup[0] - d510 * sinvalue_lookup[7] +
        d69  * sinvalue_lookup[1] + d78  * sinvalue_lookup[4];
  output[9]  = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = s015 * sinvalue_lookup[5] - s114 * sinvalue_lookup[4] -
        s213 * sinvalue_lookup[0] + s312 * sinvalue_lookup[6] -
        s411 * sinvalue_lookup[3] - s510 * sinvalue_lookup[1] +
        s69  * sinvalue_lookup[7] - s78  * sinvalue_lookup[2];
  output[10] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = d015 * sinvalue_lookup[4] - d114 * sinvalue_lookup[6] +
        d213 * sinvalue_lookup[1] + d312 * sinvalue_lookup[2] -
        d411 * sinvalue_lookup[7] + d510 * sinvalue_lookup[3] +
        d69  * sinvalue_lookup[0] - d78  * sinvalue_lookup[5];
  output[11] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = s015 * sinvalue_lookup[3] - s114 * sinvalue_lookup[7] +
        s213 * sinvalue_lookup[4] - s312 * sinvalue_lookup[0] -
        s411 * sinvalue_lookup[2] + s510 * sinvalue_lookup[6] -
        s69  * sinvalue_lookup[5] + s78  * sinvalue_lookup[1];
  output[12] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = d015 * sinvalue_lookup[2] - d114 * sinvalue_lookup[5] +
        d213 * sinvalue_lookup[7] - d312 * sinvalue_lookup[4] +
        d411 * sinvalue_lookup[1] + d510 * sinvalue_lookup[0] -
        d69  * sinvalue_lookup[3] + d78  * sinvalue_lookup[6];
  output[13] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = s015 * sinvalue_lookup[1] - s114 * sinvalue_lookup[3] +
        s213 * sinvalue_lookup[5] - s312 * sinvalue_lookup[7] +
        s411 * sinvalue_lookup[6] - s510 * sinvalue_lookup[4] +
        s69  * sinvalue_lookup[2] - s78  * sinvalue_lookup[0];
  output[14] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
  sum = d015 * sinvalue_lookup[0] - d114 * sinvalue_lookup[1] +
        d213 * sinvalue_lookup[2] - d312 * sinvalue_lookup[3] +
        d411 * sinvalue_lookup[4] - d510 * sinvalue_lookup[5] +
        d69  * sinvalue_lookup[6] - d78  * sinvalue_lookup[7];
  output[15] = ROUND_POWER_OF_TWO(sum, (2 * DCT_CONST_BITS));
#endif
}
#endif  // CONFIG_EXT_TX

void vp9_fdct4(const tran_low_t *input, tran_low_t *output) {
  tran_high_t step[4];
  tran_high_t temp1, temp2;

  step[0] = input[0] + input[3];
  step[1] = input[1] + input[2];
  step[2] = input[1] - input[2];
  step[3] = input[0] - input[3];

  temp1 = (step[0] + step[1]) * cospi_16_64;
  temp2 = (step[0] - step[1]) * cospi_16_64;
  output[0] = fdct_round_shift(temp1);
  output[2] = fdct_round_shift(temp2);
  temp1 = step[2] * cospi_24_64 + step[3] * cospi_8_64;
  temp2 = -step[2] * cospi_8_64 + step[3] * cospi_24_64;
  output[1] = fdct_round_shift(temp1);
  output[3] = fdct_round_shift(temp2);
}

void vp9_fdct4x4_1_c(const int16_t *input, tran_low_t *output, int stride) {
  int r, c;
  tran_low_t sum = 0;
  for (r = 0; r < 4; ++r)
    for (c = 0; c < 4; ++c)
      sum += input[r * stride + c];

  output[0] = sum << 1;
  output[1] = 0;
}

void vp9_fdct4x4_c(const int16_t *input, tran_low_t *output, int stride) {
  // The 2D transform is done with two passes which are actually pretty
  // similar. In the first one, we transform the columns and transpose
  // the results. In the second one, we transform the rows. To achieve that,
  // as the first pass results are transposed, we transpose the columns (that
  // is the transposed rows) and transpose the results (so that it goes back
  // in normal/row positions).
  int pass;
  // We need an intermediate buffer between passes.
  tran_low_t intermediate[4 * 4];
  const int16_t *in_pass0 = input;
  const tran_low_t *in = NULL;
  tran_low_t *out = intermediate;
  // Do the two transform/transpose passes
  for (pass = 0; pass < 2; ++pass) {
    tran_high_t input[4];      // canbe16
    tran_high_t step[4];       // canbe16
    tran_high_t temp1, temp2;  // needs32
    int i;
    for (i = 0; i < 4; ++i) {
      // Load inputs.
      if (0 == pass) {
        input[0] = in_pass0[0 * stride] * 16;
        input[1] = in_pass0[1 * stride] * 16;
        input[2] = in_pass0[2 * stride] * 16;
        input[3] = in_pass0[3 * stride] * 16;
        if (i == 0 && input[0]) {
          input[0] += 1;
        }
      } else {
        input[0] = in[0 * 4];
        input[1] = in[1 * 4];
        input[2] = in[2 * 4];
        input[3] = in[3 * 4];
      }
      // Transform.
      step[0] = input[0] + input[3];
      step[1] = input[1] + input[2];
      step[2] = input[1] - input[2];
      step[3] = input[0] - input[3];
      temp1 = (step[0] + step[1]) * cospi_16_64;
      temp2 = (step[0] - step[1]) * cospi_16_64;
      out[0] = fdct_round_shift(temp1);
      out[2] = fdct_round_shift(temp2);
      temp1 = step[2] * cospi_24_64 + step[3] * cospi_8_64;
      temp2 = -step[2] * cospi_8_64 + step[3] * cospi_24_64;
      out[1] = fdct_round_shift(temp1);
      out[3] = fdct_round_shift(temp2);
      // Do next column (which is a transposed row in second/horizontal pass)
      in_pass0++;
      in++;
      out += 4;
    }
    // Setup in/out for next pass.
    in = intermediate;
    out = output;
  }

  {
    int i, j;
    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j)
        output[j + i * 4] = (output[j + i * 4] + 1) >> 2;
    }
  }
}

void vp9_fadst4(const tran_low_t *input, tran_low_t *output) {
  tran_high_t x0, x1, x2, x3;
  tran_high_t s0, s1, s2, s3, s4, s5, s6, s7;

  x0 = input[0];
  x1 = input[1];
  x2 = input[2];
  x3 = input[3];

  if (!(x0 | x1 | x2 | x3)) {
    output[0] = output[1] = output[2] = output[3] = 0;
    return;
  }

  s0 = sinpi_1_9 * x0;
  s1 = sinpi_4_9 * x0;
  s2 = sinpi_2_9 * x1;
  s3 = sinpi_1_9 * x1;
  s4 = sinpi_3_9 * x2;
  s5 = sinpi_4_9 * x3;
  s6 = sinpi_2_9 * x3;
  s7 = x0 + x1 - x3;

  x0 = s0 + s2 + s5;
  x1 = sinpi_3_9 * s7;
  x2 = s1 - s3 + s6;
  x3 = s4;

  s0 = x0 + x3;
  s1 = x1;
  s2 = x2 - x3;
  s3 = x2 - x0 + x3;

  // 1-D transform scaling factor is sqrt(2).
  output[0] = fdct_round_shift(s0);
  output[1] = fdct_round_shift(s1);
  output[2] = fdct_round_shift(s2);
  output[3] = fdct_round_shift(s3);
}

#if CONFIG_EXT_TX
static void copy_block(const int16_t *src, int src_stride, int l,
                       int16_t *dest, int dest_stride) {
  int i;
  for (i = 0; i < l; ++i) {
    memcpy(dest + dest_stride * i, src + src_stride * i,
           l * sizeof(int16_t));
  }
}

static void fliplr(int16_t *dest, int stride, int l) {
  int i, j;
  for (i = 0; i < l; ++i) {
    for (j = 0; j < l / 2; ++j) {
      const int16_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[i * stride + l - 1 - j];
      dest[i * stride + l - 1 - j] = tmp;
    }
  }
}

static void flipud(int16_t *dest, int stride, int l) {
  int i, j;
  for (j = 0; j < l; ++j) {
    for (i = 0; i < l / 2; ++i) {
      const int16_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[(l - 1 - i) * stride + j];
      dest[(l - 1 - i) * stride + j] = tmp;
    }
  }
}

static void fliplrud(int16_t *dest, int stride, int l) {
  int i, j;
  for (i = 0; i < l / 2; ++i) {
    for (j = 0; j < l; ++j) {
      const int16_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[(l - 1 - i) * stride + l - 1 - j];
      dest[(l - 1 - i) * stride + l - 1 - j] = tmp;
    }
  }
}

static void copy_fliplr(const int16_t *src, int src_stride, int l,
                          int16_t *dest, int dest_stride) {
  copy_block(src, src_stride, l, dest, dest_stride);
  fliplr(dest, dest_stride, l);
}

static void copy_flipud(const int16_t *src, int src_stride, int l,
                          int16_t *dest, int dest_stride) {
  copy_block(src, src_stride, l, dest, dest_stride);
  flipud(dest, dest_stride, l);
}

static void copy_fliplrud(const int16_t *src, int src_stride, int l,
                            int16_t *dest, int dest_stride) {
  copy_block(src, src_stride, l, dest, dest_stride);
  fliplrud(dest, dest_stride, l);
}

static void maybe_flip_input(const int16_t **src, int *src_stride, int l,
                             int16_t *buff, int tx_type) {
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      break;
    case FLIPADST_DCT:
    case FLIPADST_ADST:
      copy_flipud(*src, *src_stride, l, buff, l);
      *src = buff;
      *src_stride = l;
      break;
    case DCT_FLIPADST:
    case ADST_FLIPADST:
      copy_fliplr(*src, *src_stride, l, buff, l);
      *src = buff;
      *src_stride = l;
      break;
    case FLIPADST_FLIPADST:
      copy_fliplrud(*src, *src_stride, l, buff, l);
      *src = buff;
      *src_stride = l;
      break;
    case DST_DST:
    case DCT_DST:
    case DST_DCT:
    case DST_ADST:
    case ADST_DST:
      break;
    case FLIPADST_DST:
      copy_flipud(*src, *src_stride, l, buff, l);
      *src = buff;
      *src_stride = l;
      break;
    case DST_FLIPADST:
      copy_fliplr(*src, *src_stride, l, buff, l);
      *src = buff;
      *src_stride = l;
      break;
    default:
      assert(0);
      break;
  }
}
#endif  // CONFIG_EXT_TX

void vp9_fht4x4_c(const int16_t *input, tran_low_t *output,
                  int stride, int tx_type) {
  if (tx_type == DCT_DCT) {
    vp9_fdct4x4_c(input, output, stride);
  } else {
    tran_low_t out[4 * 4];
    tran_low_t *outptr = &out[0];
    int i, j;
    tran_low_t temp_in[4], temp_out[4];
    const transform_2d ht = FHT_4[tx_type];

#if CONFIG_EXT_TX
    int16_t flipped_input[4 * 4];
    maybe_flip_input(&input, &stride, 4, flipped_input, tx_type);
#endif

    // Columns
    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j)
        temp_in[j] = input[j * stride + i] * 16;
      if (i == 0 && temp_in[0])
        temp_in[0] += 1;
      ht.cols(temp_in, temp_out);
      for (j = 0; j < 4; ++j)
        outptr[j * 4 + i] = temp_out[j];
    }

    // Rows
    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j)
        temp_in[j] = out[j + i * 4];
      ht.rows(temp_in, temp_out);
      for (j = 0; j < 4; ++j)
        output[j + i * 4] = (temp_out[j] + 1) >> 2;
    }
  }
}

void vp9_fdct8(const tran_low_t *input, tran_low_t *output) {
  tran_high_t s0, s1, s2, s3, s4, s5, s6, s7;  // canbe16
  tran_high_t t0, t1, t2, t3;                  // needs32
  tran_high_t x0, x1, x2, x3;                  // canbe16

  // stage 1
  s0 = input[0] + input[7];
  s1 = input[1] + input[6];
  s2 = input[2] + input[5];
  s3 = input[3] + input[4];
  s4 = input[3] - input[4];
  s5 = input[2] - input[5];
  s6 = input[1] - input[6];
  s7 = input[0] - input[7];

  // fdct4(step, step);
  x0 = s0 + s3;
  x1 = s1 + s2;
  x2 = s1 - s2;
  x3 = s0 - s3;
  t0 = (x0 + x1) * cospi_16_64;
  t1 = (x0 - x1) * cospi_16_64;
  t2 =  x2 * cospi_24_64 + x3 *  cospi_8_64;
  t3 = -x2 * cospi_8_64  + x3 * cospi_24_64;
  output[0] = fdct_round_shift(t0);
  output[2] = fdct_round_shift(t2);
  output[4] = fdct_round_shift(t1);
  output[6] = fdct_round_shift(t3);

  // Stage 2
  t0 = (s6 - s5) * cospi_16_64;
  t1 = (s6 + s5) * cospi_16_64;
  t2 = fdct_round_shift(t0);
  t3 = fdct_round_shift(t1);

  // Stage 3
  x0 = s4 + t2;
  x1 = s4 - t2;
  x2 = s7 - t3;
  x3 = s7 + t3;

  // Stage 4
  t0 = x0 * cospi_28_64 + x3 *   cospi_4_64;
  t1 = x1 * cospi_12_64 + x2 *  cospi_20_64;
  t2 = x2 * cospi_12_64 + x1 * -cospi_20_64;
  t3 = x3 * cospi_28_64 + x0 *  -cospi_4_64;
  output[1] = fdct_round_shift(t0);
  output[3] = fdct_round_shift(t2);
  output[5] = fdct_round_shift(t1);
  output[7] = fdct_round_shift(t3);
}

void vp9_fdct8x8_1_c(const int16_t *input, tran_low_t *output, int stride) {
  int r, c;
  tran_low_t sum = 0;
  for (r = 0; r < 8; ++r)
    for (c = 0; c < 8; ++c)
      sum += input[r * stride + c];

  output[0] = sum;
  output[1] = 0;
}

void vp9_fdct8x8_c(const int16_t *input, tran_low_t *final_output, int stride) {
  int i, j;
  tran_low_t intermediate[64];

  // Transform columns
  {
    tran_low_t *output = intermediate;
    tran_high_t s0, s1, s2, s3, s4, s5, s6, s7;  // canbe16
    tran_high_t t0, t1, t2, t3;                  // needs32
    tran_high_t x0, x1, x2, x3;                  // canbe16

    int i;
    for (i = 0; i < 8; i++) {
      // stage 1
      s0 = (input[0 * stride] + input[7 * stride]) * 4;
      s1 = (input[1 * stride] + input[6 * stride]) * 4;
      s2 = (input[2 * stride] + input[5 * stride]) * 4;
      s3 = (input[3 * stride] + input[4 * stride]) * 4;
      s4 = (input[3 * stride] - input[4 * stride]) * 4;
      s5 = (input[2 * stride] - input[5 * stride]) * 4;
      s6 = (input[1 * stride] - input[6 * stride]) * 4;
      s7 = (input[0 * stride] - input[7 * stride]) * 4;

      // fdct4(step, step);
      x0 = s0 + s3;
      x1 = s1 + s2;
      x2 = s1 - s2;
      x3 = s0 - s3;
      t0 = (x0 + x1) * cospi_16_64;
      t1 = (x0 - x1) * cospi_16_64;
      t2 =  x2 * cospi_24_64 + x3 *  cospi_8_64;
      t3 = -x2 * cospi_8_64  + x3 * cospi_24_64;
      output[0 * 8] = fdct_round_shift(t0);
      output[2 * 8] = fdct_round_shift(t2);
      output[4 * 8] = fdct_round_shift(t1);
      output[6 * 8] = fdct_round_shift(t3);

      // Stage 2
      t0 = (s6 - s5) * cospi_16_64;
      t1 = (s6 + s5) * cospi_16_64;
      t2 = fdct_round_shift(t0);
      t3 = fdct_round_shift(t1);

      // Stage 3
      x0 = s4 + t2;
      x1 = s4 - t2;
      x2 = s7 - t3;
      x3 = s7 + t3;

      // Stage 4
      t0 = x0 * cospi_28_64 + x3 *   cospi_4_64;
      t1 = x1 * cospi_12_64 + x2 *  cospi_20_64;
      t2 = x2 * cospi_12_64 + x1 * -cospi_20_64;
      t3 = x3 * cospi_28_64 + x0 *  -cospi_4_64;
      output[1 * 8] = fdct_round_shift(t0);
      output[3 * 8] = fdct_round_shift(t2);
      output[5 * 8] = fdct_round_shift(t1);
      output[7 * 8] = fdct_round_shift(t3);
      input++;
      output++;
    }
  }

  // Rows
  for (i = 0; i < 8; ++i) {
    vp9_fdct8(&intermediate[i * 8], &final_output[i * 8]);
    for (j = 0; j < 8; ++j)
      final_output[j + i * 8] /= 2;
  }
}

void vp9_fdct16x16_1_c(const int16_t *input, tran_low_t *output, int stride) {
  int r, c;
  tran_low_t sum = 0;
  for (r = 0; r < 16; ++r)
    for (c = 0; c < 16; ++c)
      sum += input[r * stride + c];

  output[0] = sum >> 1;
  output[1] = 0;
}

void vp9_fdct16x16_c(const int16_t *input, tran_low_t *output, int stride) {
  // The 2D transform is done with two passes which are actually pretty
  // similar. In the first one, we transform the columns and transpose
  // the results. In the second one, we transform the rows. To achieve that,
  // as the first pass results are transposed, we transpose the columns (that
  // is the transposed rows) and transpose the results (so that it goes back
  // in normal/row positions).
  int pass;
  // We need an intermediate buffer between passes.
  tran_low_t intermediate[256];
  const int16_t *in_pass0 = input;
  const tran_low_t *in = NULL;
  tran_low_t *out = intermediate;
  // Do the two transform/transpose passes
  for (pass = 0; pass < 2; ++pass) {
    tran_high_t step1[8];      // canbe16
    tran_high_t step2[8];      // canbe16
    tran_high_t step3[8];      // canbe16
    tran_high_t input[8];      // canbe16
    tran_high_t temp1, temp2;  // needs32
    int i;
    for (i = 0; i < 16; i++) {
      if (0 == pass) {
        // Calculate input for the first 8 results.
        input[0] = (in_pass0[0 * stride] + in_pass0[15 * stride]) * 4;
        input[1] = (in_pass0[1 * stride] + in_pass0[14 * stride]) * 4;
        input[2] = (in_pass0[2 * stride] + in_pass0[13 * stride]) * 4;
        input[3] = (in_pass0[3 * stride] + in_pass0[12 * stride]) * 4;
        input[4] = (in_pass0[4 * stride] + in_pass0[11 * stride]) * 4;
        input[5] = (in_pass0[5 * stride] + in_pass0[10 * stride]) * 4;
        input[6] = (in_pass0[6 * stride] + in_pass0[ 9 * stride]) * 4;
        input[7] = (in_pass0[7 * stride] + in_pass0[ 8 * stride]) * 4;
        // Calculate input for the next 8 results.
        step1[0] = (in_pass0[7 * stride] - in_pass0[ 8 * stride]) * 4;
        step1[1] = (in_pass0[6 * stride] - in_pass0[ 9 * stride]) * 4;
        step1[2] = (in_pass0[5 * stride] - in_pass0[10 * stride]) * 4;
        step1[3] = (in_pass0[4 * stride] - in_pass0[11 * stride]) * 4;
        step1[4] = (in_pass0[3 * stride] - in_pass0[12 * stride]) * 4;
        step1[5] = (in_pass0[2 * stride] - in_pass0[13 * stride]) * 4;
        step1[6] = (in_pass0[1 * stride] - in_pass0[14 * stride]) * 4;
        step1[7] = (in_pass0[0 * stride] - in_pass0[15 * stride]) * 4;
      } else {
        // Calculate input for the first 8 results.
        input[0] = ((in[0 * 16] + 1) >> 2) + ((in[15 * 16] + 1) >> 2);
        input[1] = ((in[1 * 16] + 1) >> 2) + ((in[14 * 16] + 1) >> 2);
        input[2] = ((in[2 * 16] + 1) >> 2) + ((in[13 * 16] + 1) >> 2);
        input[3] = ((in[3 * 16] + 1) >> 2) + ((in[12 * 16] + 1) >> 2);
        input[4] = ((in[4 * 16] + 1) >> 2) + ((in[11 * 16] + 1) >> 2);
        input[5] = ((in[5 * 16] + 1) >> 2) + ((in[10 * 16] + 1) >> 2);
        input[6] = ((in[6 * 16] + 1) >> 2) + ((in[ 9 * 16] + 1) >> 2);
        input[7] = ((in[7 * 16] + 1) >> 2) + ((in[ 8 * 16] + 1) >> 2);
        // Calculate input for the next 8 results.
        step1[0] = ((in[7 * 16] + 1) >> 2) - ((in[ 8 * 16] + 1) >> 2);
        step1[1] = ((in[6 * 16] + 1) >> 2) - ((in[ 9 * 16] + 1) >> 2);
        step1[2] = ((in[5 * 16] + 1) >> 2) - ((in[10 * 16] + 1) >> 2);
        step1[3] = ((in[4 * 16] + 1) >> 2) - ((in[11 * 16] + 1) >> 2);
        step1[4] = ((in[3 * 16] + 1) >> 2) - ((in[12 * 16] + 1) >> 2);
        step1[5] = ((in[2 * 16] + 1) >> 2) - ((in[13 * 16] + 1) >> 2);
        step1[6] = ((in[1 * 16] + 1) >> 2) - ((in[14 * 16] + 1) >> 2);
        step1[7] = ((in[0 * 16] + 1) >> 2) - ((in[15 * 16] + 1) >> 2);
      }
      // Work on the first eight values; fdct8(input, even_results);
      {
        tran_high_t s0, s1, s2, s3, s4, s5, s6, s7;  // canbe16
        tran_high_t t0, t1, t2, t3;                  // needs32
        tran_high_t x0, x1, x2, x3;                  // canbe16

        // stage 1
        s0 = input[0] + input[7];
        s1 = input[1] + input[6];
        s2 = input[2] + input[5];
        s3 = input[3] + input[4];
        s4 = input[3] - input[4];
        s5 = input[2] - input[5];
        s6 = input[1] - input[6];
        s7 = input[0] - input[7];

        // fdct4(step, step);
        x0 = s0 + s3;
        x1 = s1 + s2;
        x2 = s1 - s2;
        x3 = s0 - s3;
        t0 = (x0 + x1) * cospi_16_64;
        t1 = (x0 - x1) * cospi_16_64;
        t2 = x3 * cospi_8_64  + x2 * cospi_24_64;
        t3 = x3 * cospi_24_64 - x2 * cospi_8_64;
        out[0] = fdct_round_shift(t0);
        out[4] = fdct_round_shift(t2);
        out[8] = fdct_round_shift(t1);
        out[12] = fdct_round_shift(t3);

        // Stage 2
        t0 = (s6 - s5) * cospi_16_64;
        t1 = (s6 + s5) * cospi_16_64;
        t2 = fdct_round_shift(t0);
        t3 = fdct_round_shift(t1);

        // Stage 3
        x0 = s4 + t2;
        x1 = s4 - t2;
        x2 = s7 - t3;
        x3 = s7 + t3;

        // Stage 4
        t0 = x0 * cospi_28_64 + x3 *   cospi_4_64;
        t1 = x1 * cospi_12_64 + x2 *  cospi_20_64;
        t2 = x2 * cospi_12_64 + x1 * -cospi_20_64;
        t3 = x3 * cospi_28_64 + x0 *  -cospi_4_64;
        out[2] = fdct_round_shift(t0);
        out[6] = fdct_round_shift(t2);
        out[10] = fdct_round_shift(t1);
        out[14] = fdct_round_shift(t3);
      }
      // Work on the next eight values; step1 -> odd_results
      {
        // step 2
        temp1 = (step1[5] - step1[2]) * cospi_16_64;
        temp2 = (step1[4] - step1[3]) * cospi_16_64;
        step2[2] = fdct_round_shift(temp1);
        step2[3] = fdct_round_shift(temp2);
        temp1 = (step1[4] + step1[3]) * cospi_16_64;
        temp2 = (step1[5] + step1[2]) * cospi_16_64;
        step2[4] = fdct_round_shift(temp1);
        step2[5] = fdct_round_shift(temp2);
        // step 3
        step3[0] = step1[0] + step2[3];
        step3[1] = step1[1] + step2[2];
        step3[2] = step1[1] - step2[2];
        step3[3] = step1[0] - step2[3];
        step3[4] = step1[7] - step2[4];
        step3[5] = step1[6] - step2[5];
        step3[6] = step1[6] + step2[5];
        step3[7] = step1[7] + step2[4];
        // step 4
        temp1 = step3[1] *  -cospi_8_64 + step3[6] * cospi_24_64;
        temp2 = step3[2] * cospi_24_64 + step3[5] *  cospi_8_64;
        step2[1] = fdct_round_shift(temp1);
        step2[2] = fdct_round_shift(temp2);
        temp1 = step3[2] * cospi_8_64 - step3[5] * cospi_24_64;
        temp2 = step3[1] * cospi_24_64 + step3[6] *  cospi_8_64;
        step2[5] = fdct_round_shift(temp1);
        step2[6] = fdct_round_shift(temp2);
        // step 5
        step1[0] = step3[0] + step2[1];
        step1[1] = step3[0] - step2[1];
        step1[2] = step3[3] + step2[2];
        step1[3] = step3[3] - step2[2];
        step1[4] = step3[4] - step2[5];
        step1[5] = step3[4] + step2[5];
        step1[6] = step3[7] - step2[6];
        step1[7] = step3[7] + step2[6];
        // step 6
        temp1 = step1[0] * cospi_30_64 + step1[7] *  cospi_2_64;
        temp2 = step1[1] * cospi_14_64 + step1[6] * cospi_18_64;
        out[1] = fdct_round_shift(temp1);
        out[9] = fdct_round_shift(temp2);
        temp1 = step1[2] * cospi_22_64 + step1[5] * cospi_10_64;
        temp2 = step1[3] *  cospi_6_64 + step1[4] * cospi_26_64;
        out[5] = fdct_round_shift(temp1);
        out[13] = fdct_round_shift(temp2);
        temp1 = step1[3] * -cospi_26_64 + step1[4] *  cospi_6_64;
        temp2 = step1[2] * -cospi_10_64 + step1[5] * cospi_22_64;
        out[3] = fdct_round_shift(temp1);
        out[11] = fdct_round_shift(temp2);
        temp1 = step1[1] * -cospi_18_64 + step1[6] * cospi_14_64;
        temp2 = step1[0] *  -cospi_2_64 + step1[7] * cospi_30_64;
        out[7] = fdct_round_shift(temp1);
        out[15] = fdct_round_shift(temp2);
      }
      // Do next column (which is a transposed row in second/horizontal pass)
      in++;
      in_pass0++;
      out += 16;
    }
    // Setup in/out for next pass.
    in = intermediate;
    out = output;
  }
}

#if CONFIG_WAVELETS
// The difference between this one and the function above is scaling
// of the input. This function does not scale so that the actual 2D
// transform is unitary. The function above scales the transform to be
// 8 times unitary.
void vp9_fdct16x16_noscale_c(const int16_t *input, tran_low_t *output,
                             int stride) {
  // The 2D transform is done with two passes which are actually pretty
  // similar. In the first one, we transform the columns and transpose
  // the results. In the second one, we transform the rows. To achieve that,
  // as the first pass results are transposed, we transpose the columns (that
  // is the transposed rows) and transpose the results (so that it goes back
  // in normal/row positions).
  int pass;
  // We need an intermediate buffer between passes.
  tran_low_t intermediate[256];
  const int16_t *in_pass0 = input;
  const tran_low_t *in = NULL;
  tran_low_t *out = intermediate;
  // Do the two transform/transpose passes
  for (pass = 0; pass < 2; ++pass) {
    tran_high_t step1[8];      // canbe16
    tran_high_t step2[8];      // canbe16
    tran_high_t step3[8];      // canbe16
    tran_high_t input[8];      // canbe16
    tran_high_t temp1, temp2;  // needs32
    int i;
    for (i = 0; i < 16; i++) {
      if (0 == pass) {
        // Calculate input for the first 8 results.
        input[0] = (in_pass0[0 * stride] + in_pass0[15 * stride]) >> 1;
        input[1] = (in_pass0[1 * stride] + in_pass0[14 * stride]) >> 1;
        input[2] = (in_pass0[2 * stride] + in_pass0[13 * stride]) >> 1;
        input[3] = (in_pass0[3 * stride] + in_pass0[12 * stride]) >> 1;
        input[4] = (in_pass0[4 * stride] + in_pass0[11 * stride]) >> 1;
        input[5] = (in_pass0[5 * stride] + in_pass0[10 * stride]) >> 1;
        input[6] = (in_pass0[6 * stride] + in_pass0[ 9 * stride]) >> 1;
        input[7] = (in_pass0[7 * stride] + in_pass0[ 8 * stride]) >> 1;
        // Calculate input for the next 8 results.
        step1[0] = (in_pass0[7 * stride] - in_pass0[ 8 * stride]) >> 1;
        step1[1] = (in_pass0[6 * stride] - in_pass0[ 9 * stride]) >> 1;
        step1[2] = (in_pass0[5 * stride] - in_pass0[10 * stride]) >> 1;
        step1[3] = (in_pass0[4 * stride] - in_pass0[11 * stride]) >> 1;
        step1[4] = (in_pass0[3 * stride] - in_pass0[12 * stride]) >> 1;
        step1[5] = (in_pass0[2 * stride] - in_pass0[13 * stride]) >> 1;
        step1[6] = (in_pass0[1 * stride] - in_pass0[14 * stride]) >> 1;
        step1[7] = (in_pass0[0 * stride] - in_pass0[15 * stride]) >> 1;
      } else {
        // Calculate input for the first 8 results.
        input[0] = ((in[0 * 16] + 1) >> 2) + ((in[15 * 16] + 1) >> 2);
        input[1] = ((in[1 * 16] + 1) >> 2) + ((in[14 * 16] + 1) >> 2);
        input[2] = ((in[2 * 16] + 1) >> 2) + ((in[13 * 16] + 1) >> 2);
        input[3] = ((in[3 * 16] + 1) >> 2) + ((in[12 * 16] + 1) >> 2);
        input[4] = ((in[4 * 16] + 1) >> 2) + ((in[11 * 16] + 1) >> 2);
        input[5] = ((in[5 * 16] + 1) >> 2) + ((in[10 * 16] + 1) >> 2);
        input[6] = ((in[6 * 16] + 1) >> 2) + ((in[ 9 * 16] + 1) >> 2);
        input[7] = ((in[7 * 16] + 1) >> 2) + ((in[ 8 * 16] + 1) >> 2);
        // Calculate input for the next 8 results.
        step1[0] = ((in[7 * 16] + 1) >> 2) - ((in[ 8 * 16] + 1) >> 2);
        step1[1] = ((in[6 * 16] + 1) >> 2) - ((in[ 9 * 16] + 1) >> 2);
        step1[2] = ((in[5 * 16] + 1) >> 2) - ((in[10 * 16] + 1) >> 2);
        step1[3] = ((in[4 * 16] + 1) >> 2) - ((in[11 * 16] + 1) >> 2);
        step1[4] = ((in[3 * 16] + 1) >> 2) - ((in[12 * 16] + 1) >> 2);
        step1[5] = ((in[2 * 16] + 1) >> 2) - ((in[13 * 16] + 1) >> 2);
        step1[6] = ((in[1 * 16] + 1) >> 2) - ((in[14 * 16] + 1) >> 2);
        step1[7] = ((in[0 * 16] + 1) >> 2) - ((in[15 * 16] + 1) >> 2);
      }
      // Work on the first eight values; fdct8(input, even_results);
      {
        tran_high_t s0, s1, s2, s3, s4, s5, s6, s7;  // canbe16
        tran_high_t t0, t1, t2, t3;                  // needs32
        tran_high_t x0, x1, x2, x3;                  // canbe16

        // stage 1
        s0 = input[0] + input[7];
        s1 = input[1] + input[6];
        s2 = input[2] + input[5];
        s3 = input[3] + input[4];
        s4 = input[3] - input[4];
        s5 = input[2] - input[5];
        s6 = input[1] - input[6];
        s7 = input[0] - input[7];

        // fdct4(step, step);
        x0 = s0 + s3;
        x1 = s1 + s2;
        x2 = s1 - s2;
        x3 = s0 - s3;
        t0 = (x0 + x1) * cospi_16_64;
        t1 = (x0 - x1) * cospi_16_64;
        t2 = x3 * cospi_8_64  + x2 * cospi_24_64;
        t3 = x3 * cospi_24_64 - x2 * cospi_8_64;
        out[0] = fdct_round_shift(t0);
        out[4] = fdct_round_shift(t2);
        out[8] = fdct_round_shift(t1);
        out[12] = fdct_round_shift(t3);

        // Stage 2
        t0 = (s6 - s5) * cospi_16_64;
        t1 = (s6 + s5) * cospi_16_64;
        t2 = fdct_round_shift(t0);
        t3 = fdct_round_shift(t1);

        // Stage 3
        x0 = s4 + t2;
        x1 = s4 - t2;
        x2 = s7 - t3;
        x3 = s7 + t3;

        // Stage 4
        t0 = x0 * cospi_28_64 + x3 *   cospi_4_64;
        t1 = x1 * cospi_12_64 + x2 *  cospi_20_64;
        t2 = x2 * cospi_12_64 + x1 * -cospi_20_64;
        t3 = x3 * cospi_28_64 + x0 *  -cospi_4_64;
        out[2] = fdct_round_shift(t0);
        out[6] = fdct_round_shift(t2);
        out[10] = fdct_round_shift(t1);
        out[14] = fdct_round_shift(t3);
      }
      // Work on the next eight values; step1 -> odd_results
      {
        // step 2
        temp1 = (step1[5] - step1[2]) * cospi_16_64;
        temp2 = (step1[4] - step1[3]) * cospi_16_64;
        step2[2] = fdct_round_shift(temp1);
        step2[3] = fdct_round_shift(temp2);
        temp1 = (step1[4] + step1[3]) * cospi_16_64;
        temp2 = (step1[5] + step1[2]) * cospi_16_64;
        step2[4] = fdct_round_shift(temp1);
        step2[5] = fdct_round_shift(temp2);
        // step 3
        step3[0] = step1[0] + step2[3];
        step3[1] = step1[1] + step2[2];
        step3[2] = step1[1] - step2[2];
        step3[3] = step1[0] - step2[3];
        step3[4] = step1[7] - step2[4];
        step3[5] = step1[6] - step2[5];
        step3[6] = step1[6] + step2[5];
        step3[7] = step1[7] + step2[4];
        // step 4
        temp1 = step3[1] *  -cospi_8_64 + step3[6] * cospi_24_64;
        temp2 = step3[2] * cospi_24_64 + step3[5] *  cospi_8_64;
        step2[1] = fdct_round_shift(temp1);
        step2[2] = fdct_round_shift(temp2);
        temp1 = step3[2] * cospi_8_64 - step3[5] * cospi_24_64;
        temp2 = step3[1] * cospi_24_64 + step3[6] *  cospi_8_64;
        step2[5] = fdct_round_shift(temp1);
        step2[6] = fdct_round_shift(temp2);
        // step 5
        step1[0] = step3[0] + step2[1];
        step1[1] = step3[0] - step2[1];
        step1[2] = step3[3] + step2[2];
        step1[3] = step3[3] - step2[2];
        step1[4] = step3[4] - step2[5];
        step1[5] = step3[4] + step2[5];
        step1[6] = step3[7] - step2[6];
        step1[7] = step3[7] + step2[6];
        // step 6
        temp1 = step1[0] * cospi_30_64 + step1[7] *  cospi_2_64;
        temp2 = step1[1] * cospi_14_64 + step1[6] * cospi_18_64;
        out[1] = fdct_round_shift(temp1);
        out[9] = fdct_round_shift(temp2);
        temp1 = step1[2] * cospi_22_64 + step1[5] * cospi_10_64;
        temp2 = step1[3] *  cospi_6_64 + step1[4] * cospi_26_64;
        out[5] = fdct_round_shift(temp1);
        out[13] = fdct_round_shift(temp2);
        temp1 = step1[3] * -cospi_26_64 + step1[4] *  cospi_6_64;
        temp2 = step1[2] * -cospi_10_64 + step1[5] * cospi_22_64;
        out[3] = fdct_round_shift(temp1);
        out[11] = fdct_round_shift(temp2);
        temp1 = step1[1] * -cospi_18_64 + step1[6] * cospi_14_64;
        temp2 = step1[0] *  -cospi_2_64 + step1[7] * cospi_30_64;
        out[7] = fdct_round_shift(temp1);
        out[15] = fdct_round_shift(temp2);
      }
      // Do next column (which is a transposed row in second/horizontal pass)
      in++;
      in_pass0++;
      out += 16;
    }
    // Setup in/out for next pass.
    in = intermediate;
    out = output;
  }
}
#endif  // CONFIG_WAVELETS

void vp9_fadst8(const tran_low_t *input, tran_low_t *output) {
  tran_high_t s0, s1, s2, s3, s4, s5, s6, s7;

  tran_high_t x0 = input[7];
  tran_high_t x1 = input[0];
  tran_high_t x2 = input[5];
  tran_high_t x3 = input[2];
  tran_high_t x4 = input[3];
  tran_high_t x5 = input[4];
  tran_high_t x6 = input[1];
  tran_high_t x7 = input[6];

  // stage 1
  s0 = cospi_2_64  * x0 + cospi_30_64 * x1;
  s1 = cospi_30_64 * x0 - cospi_2_64  * x1;
  s2 = cospi_10_64 * x2 + cospi_22_64 * x3;
  s3 = cospi_22_64 * x2 - cospi_10_64 * x3;
  s4 = cospi_18_64 * x4 + cospi_14_64 * x5;
  s5 = cospi_14_64 * x4 - cospi_18_64 * x5;
  s6 = cospi_26_64 * x6 + cospi_6_64  * x7;
  s7 = cospi_6_64  * x6 - cospi_26_64 * x7;

  x0 = fdct_round_shift(s0 + s4);
  x1 = fdct_round_shift(s1 + s5);
  x2 = fdct_round_shift(s2 + s6);
  x3 = fdct_round_shift(s3 + s7);
  x4 = fdct_round_shift(s0 - s4);
  x5 = fdct_round_shift(s1 - s5);
  x6 = fdct_round_shift(s2 - s6);
  x7 = fdct_round_shift(s3 - s7);

  // stage 2
  s0 = x0;
  s1 = x1;
  s2 = x2;
  s3 = x3;
  s4 = cospi_8_64  * x4 + cospi_24_64 * x5;
  s5 = cospi_24_64 * x4 - cospi_8_64  * x5;
  s6 = - cospi_24_64 * x6 + cospi_8_64  * x7;
  s7 =   cospi_8_64  * x6 + cospi_24_64 * x7;

  x0 = s0 + s2;
  x1 = s1 + s3;
  x2 = s0 - s2;
  x3 = s1 - s3;
  x4 = fdct_round_shift(s4 + s6);
  x5 = fdct_round_shift(s5 + s7);
  x6 = fdct_round_shift(s4 - s6);
  x7 = fdct_round_shift(s5 - s7);

  // stage 3
  s2 = cospi_16_64 * (x2 + x3);
  s3 = cospi_16_64 * (x2 - x3);
  s6 = cospi_16_64 * (x6 + x7);
  s7 = cospi_16_64 * (x6 - x7);

  x2 = fdct_round_shift(s2);
  x3 = fdct_round_shift(s3);
  x6 = fdct_round_shift(s6);
  x7 = fdct_round_shift(s7);

  output[0] =   x0;
  output[1] = - x4;
  output[2] =   x6;
  output[3] = - x2;
  output[4] =   x3;
  output[5] = - x7;
  output[6] =   x5;
  output[7] = - x1;
}

void vp9_fht8x8_c(const int16_t *input, tran_low_t *output,
                  int stride, int tx_type) {
  if (tx_type == DCT_DCT) {
    vp9_fdct8x8_c(input, output, stride);
  } else {
    tran_low_t out[64];
    tran_low_t *outptr = &out[0];
    int i, j;
    tran_low_t temp_in[8], temp_out[8];
    const transform_2d ht = FHT_8[tx_type];

#if CONFIG_EXT_TX
    int16_t flipped_input[8 * 8];
    maybe_flip_input(&input, &stride, 8, flipped_input, tx_type);
#endif

    // Columns
    for (i = 0; i < 8; ++i) {
      for (j = 0; j < 8; ++j)
        temp_in[j] = input[j * stride + i] * 4;
      ht.cols(temp_in, temp_out);
      for (j = 0; j < 8; ++j)
        outptr[j * 8 + i] = temp_out[j];
    }

    // Rows
    for (i = 0; i < 8; ++i) {
      for (j = 0; j < 8; ++j)
        temp_in[j] = out[j + i * 8];
      ht.rows(temp_in, temp_out);
      for (j = 0; j < 8; ++j)
        output[j + i * 8] = (temp_out[j] + (temp_out[j] < 0)) >> 1;
    }
  }
}

/* 4-point reversible, orthonormal Walsh-Hadamard in 3.5 adds, 0.5 shifts per
   pixel. */
void vp9_fwht4x4_c(const int16_t *input, tran_low_t *output, int stride) {
  int i;
  tran_high_t a1, b1, c1, d1, e1;
  const int16_t *ip_pass0 = input;
  const tran_low_t *ip = NULL;
  tran_low_t *op = output;

  for (i = 0; i < 4; i++) {
    a1 = ip_pass0[0 * stride];
    b1 = ip_pass0[1 * stride];
    c1 = ip_pass0[2 * stride];
    d1 = ip_pass0[3 * stride];

    a1 += b1;
    d1 = d1 - c1;
    e1 = (a1 - d1) >> 1;
    b1 = e1 - b1;
    c1 = e1 - c1;
    a1 -= c1;
    d1 += b1;
    op[0] = a1;
    op[4] = c1;
    op[8] = d1;
    op[12] = b1;

    ip_pass0++;
    op++;
  }
  ip = output;
  op = output;

  for (i = 0; i < 4; i++) {
    a1 = ip[0];
    b1 = ip[1];
    c1 = ip[2];
    d1 = ip[3];

    a1 += b1;
    d1 -= c1;
    e1 = (a1 - d1) >> 1;
    b1 = e1 - b1;
    c1 = e1 - c1;
    a1 -= c1;
    d1 += b1;
    op[0] = a1 * UNIT_QUANT_FACTOR;
    op[1] = c1 * UNIT_QUANT_FACTOR;
    op[2] = d1 * UNIT_QUANT_FACTOR;
    op[3] = b1 * UNIT_QUANT_FACTOR;

    ip += 4;
    op += 4;
  }
}

// Rewrote to use same algorithm as others.
void vp9_fdct16(const tran_low_t in[16], tran_low_t out[16]) {
  tran_high_t step1[8];      // canbe16
  tran_high_t step2[8];      // canbe16
  tran_high_t step3[8];      // canbe16
  tran_high_t input[8];      // canbe16
  tran_high_t temp1, temp2;  // needs32

  // step 1
  input[0] = in[0] + in[15];
  input[1] = in[1] + in[14];
  input[2] = in[2] + in[13];
  input[3] = in[3] + in[12];
  input[4] = in[4] + in[11];
  input[5] = in[5] + in[10];
  input[6] = in[6] + in[ 9];
  input[7] = in[7] + in[ 8];

  step1[0] = in[7] - in[ 8];
  step1[1] = in[6] - in[ 9];
  step1[2] = in[5] - in[10];
  step1[3] = in[4] - in[11];
  step1[4] = in[3] - in[12];
  step1[5] = in[2] - in[13];
  step1[6] = in[1] - in[14];
  step1[7] = in[0] - in[15];

  // fdct8(step, step);
  {
    tran_high_t s0, s1, s2, s3, s4, s5, s6, s7;  // canbe16
    tran_high_t t0, t1, t2, t3;                  // needs32
    tran_high_t x0, x1, x2, x3;                  // canbe16

    // stage 1
    s0 = input[0] + input[7];
    s1 = input[1] + input[6];
    s2 = input[2] + input[5];
    s3 = input[3] + input[4];
    s4 = input[3] - input[4];
    s5 = input[2] - input[5];
    s6 = input[1] - input[6];
    s7 = input[0] - input[7];

    // fdct4(step, step);
    x0 = s0 + s3;
    x1 = s1 + s2;
    x2 = s1 - s2;
    x3 = s0 - s3;
    t0 = (x0 + x1) * cospi_16_64;
    t1 = (x0 - x1) * cospi_16_64;
    t2 = x3 * cospi_8_64  + x2 * cospi_24_64;
    t3 = x3 * cospi_24_64 - x2 * cospi_8_64;
    out[0] = fdct_round_shift(t0);
    out[4] = fdct_round_shift(t2);
    out[8] = fdct_round_shift(t1);
    out[12] = fdct_round_shift(t3);

    // Stage 2
    t0 = (s6 - s5) * cospi_16_64;
    t1 = (s6 + s5) * cospi_16_64;
    t2 = fdct_round_shift(t0);
    t3 = fdct_round_shift(t1);

    // Stage 3
    x0 = s4 + t2;
    x1 = s4 - t2;
    x2 = s7 - t3;
    x3 = s7 + t3;

    // Stage 4
    t0 = x0 * cospi_28_64 + x3 *   cospi_4_64;
    t1 = x1 * cospi_12_64 + x2 *  cospi_20_64;
    t2 = x2 * cospi_12_64 + x1 * -cospi_20_64;
    t3 = x3 * cospi_28_64 + x0 *  -cospi_4_64;
    out[2] = fdct_round_shift(t0);
    out[6] = fdct_round_shift(t2);
    out[10] = fdct_round_shift(t1);
    out[14] = fdct_round_shift(t3);
  }

  // step 2
  temp1 = (step1[5] - step1[2]) * cospi_16_64;
  temp2 = (step1[4] - step1[3]) * cospi_16_64;
  step2[2] = fdct_round_shift(temp1);
  step2[3] = fdct_round_shift(temp2);
  temp1 = (step1[4] + step1[3]) * cospi_16_64;
  temp2 = (step1[5] + step1[2]) * cospi_16_64;
  step2[4] = fdct_round_shift(temp1);
  step2[5] = fdct_round_shift(temp2);

  // step 3
  step3[0] = step1[0] + step2[3];
  step3[1] = step1[1] + step2[2];
  step3[2] = step1[1] - step2[2];
  step3[3] = step1[0] - step2[3];
  step3[4] = step1[7] - step2[4];
  step3[5] = step1[6] - step2[5];
  step3[6] = step1[6] + step2[5];
  step3[7] = step1[7] + step2[4];

  // step 4
  temp1 = step3[1] *  -cospi_8_64 + step3[6] * cospi_24_64;
  temp2 = step3[2] * cospi_24_64 + step3[5] *  cospi_8_64;
  step2[1] = fdct_round_shift(temp1);
  step2[2] = fdct_round_shift(temp2);
  temp1 = step3[2] * cospi_8_64 - step3[5] * cospi_24_64;
  temp2 = step3[1] * cospi_24_64 + step3[6] *  cospi_8_64;
  step2[5] = fdct_round_shift(temp1);
  step2[6] = fdct_round_shift(temp2);

  // step 5
  step1[0] = step3[0] + step2[1];
  step1[1] = step3[0] - step2[1];
  step1[2] = step3[3] + step2[2];
  step1[3] = step3[3] - step2[2];
  step1[4] = step3[4] - step2[5];
  step1[5] = step3[4] + step2[5];
  step1[6] = step3[7] - step2[6];
  step1[7] = step3[7] + step2[6];

  // step 6
  temp1 = step1[0] * cospi_30_64 + step1[7] *  cospi_2_64;
  temp2 = step1[1] * cospi_14_64 + step1[6] * cospi_18_64;
  out[1] = fdct_round_shift(temp1);
  out[9] = fdct_round_shift(temp2);

  temp1 = step1[2] * cospi_22_64 + step1[5] * cospi_10_64;
  temp2 = step1[3] *  cospi_6_64 + step1[4] * cospi_26_64;
  out[5] = fdct_round_shift(temp1);
  out[13] = fdct_round_shift(temp2);

  temp1 = step1[3] * -cospi_26_64 + step1[4] *  cospi_6_64;
  temp2 = step1[2] * -cospi_10_64 + step1[5] * cospi_22_64;
  out[3] = fdct_round_shift(temp1);
  out[11] = fdct_round_shift(temp2);

  temp1 = step1[1] * -cospi_18_64 + step1[6] * cospi_14_64;
  temp2 = step1[0] *  -cospi_2_64 + step1[7] * cospi_30_64;
  out[7] = fdct_round_shift(temp1);
  out[15] = fdct_round_shift(temp2);
}

void vp9_fadst16(const tran_low_t *input, tran_low_t *output) {
  tran_high_t s0, s1, s2, s3, s4, s5, s6, s7, s8;
  tran_high_t s9, s10, s11, s12, s13, s14, s15;

  tran_high_t x0 = input[15];
  tran_high_t x1 = input[0];
  tran_high_t x2 = input[13];
  tran_high_t x3 = input[2];
  tran_high_t x4 = input[11];
  tran_high_t x5 = input[4];
  tran_high_t x6 = input[9];
  tran_high_t x7 = input[6];
  tran_high_t x8 = input[7];
  tran_high_t x9 = input[8];
  tran_high_t x10 = input[5];
  tran_high_t x11 = input[10];
  tran_high_t x12 = input[3];
  tran_high_t x13 = input[12];
  tran_high_t x14 = input[1];
  tran_high_t x15 = input[14];

  // stage 1
  s0 = x0 * cospi_1_64  + x1 * cospi_31_64;
  s1 = x0 * cospi_31_64 - x1 * cospi_1_64;
  s2 = x2 * cospi_5_64  + x3 * cospi_27_64;
  s3 = x2 * cospi_27_64 - x3 * cospi_5_64;
  s4 = x4 * cospi_9_64  + x5 * cospi_23_64;
  s5 = x4 * cospi_23_64 - x5 * cospi_9_64;
  s6 = x6 * cospi_13_64 + x7 * cospi_19_64;
  s7 = x6 * cospi_19_64 - x7 * cospi_13_64;
  s8 = x8 * cospi_17_64 + x9 * cospi_15_64;
  s9 = x8 * cospi_15_64 - x9 * cospi_17_64;
  s10 = x10 * cospi_21_64 + x11 * cospi_11_64;
  s11 = x10 * cospi_11_64 - x11 * cospi_21_64;
  s12 = x12 * cospi_25_64 + x13 * cospi_7_64;
  s13 = x12 * cospi_7_64  - x13 * cospi_25_64;
  s14 = x14 * cospi_29_64 + x15 * cospi_3_64;
  s15 = x14 * cospi_3_64  - x15 * cospi_29_64;

  x0 = fdct_round_shift(s0 + s8);
  x1 = fdct_round_shift(s1 + s9);
  x2 = fdct_round_shift(s2 + s10);
  x3 = fdct_round_shift(s3 + s11);
  x4 = fdct_round_shift(s4 + s12);
  x5 = fdct_round_shift(s5 + s13);
  x6 = fdct_round_shift(s6 + s14);
  x7 = fdct_round_shift(s7 + s15);
  x8  = fdct_round_shift(s0 - s8);
  x9  = fdct_round_shift(s1 - s9);
  x10 = fdct_round_shift(s2 - s10);
  x11 = fdct_round_shift(s3 - s11);
  x12 = fdct_round_shift(s4 - s12);
  x13 = fdct_round_shift(s5 - s13);
  x14 = fdct_round_shift(s6 - s14);
  x15 = fdct_round_shift(s7 - s15);

  // stage 2
  s0 = x0;
  s1 = x1;
  s2 = x2;
  s3 = x3;
  s4 = x4;
  s5 = x5;
  s6 = x6;
  s7 = x7;
  s8 =    x8 * cospi_4_64   + x9 * cospi_28_64;
  s9 =    x8 * cospi_28_64  - x9 * cospi_4_64;
  s10 =   x10 * cospi_20_64 + x11 * cospi_12_64;
  s11 =   x10 * cospi_12_64 - x11 * cospi_20_64;
  s12 = - x12 * cospi_28_64 + x13 * cospi_4_64;
  s13 =   x12 * cospi_4_64  + x13 * cospi_28_64;
  s14 = - x14 * cospi_12_64 + x15 * cospi_20_64;
  s15 =   x14 * cospi_20_64 + x15 * cospi_12_64;

  x0 = s0 + s4;
  x1 = s1 + s5;
  x2 = s2 + s6;
  x3 = s3 + s7;
  x4 = s0 - s4;
  x5 = s1 - s5;
  x6 = s2 - s6;
  x7 = s3 - s7;
  x8 = fdct_round_shift(s8 + s12);
  x9 = fdct_round_shift(s9 + s13);
  x10 = fdct_round_shift(s10 + s14);
  x11 = fdct_round_shift(s11 + s15);
  x12 = fdct_round_shift(s8 - s12);
  x13 = fdct_round_shift(s9 - s13);
  x14 = fdct_round_shift(s10 - s14);
  x15 = fdct_round_shift(s11 - s15);

  // stage 3
  s0 = x0;
  s1 = x1;
  s2 = x2;
  s3 = x3;
  s4 = x4 * cospi_8_64  + x5 * cospi_24_64;
  s5 = x4 * cospi_24_64 - x5 * cospi_8_64;
  s6 = - x6 * cospi_24_64 + x7 * cospi_8_64;
  s7 =   x6 * cospi_8_64  + x7 * cospi_24_64;
  s8 = x8;
  s9 = x9;
  s10 = x10;
  s11 = x11;
  s12 = x12 * cospi_8_64  + x13 * cospi_24_64;
  s13 = x12 * cospi_24_64 - x13 * cospi_8_64;
  s14 = - x14 * cospi_24_64 + x15 * cospi_8_64;
  s15 =   x14 * cospi_8_64  + x15 * cospi_24_64;

  x0 = s0 + s2;
  x1 = s1 + s3;
  x2 = s0 - s2;
  x3 = s1 - s3;
  x4 = fdct_round_shift(s4 + s6);
  x5 = fdct_round_shift(s5 + s7);
  x6 = fdct_round_shift(s4 - s6);
  x7 = fdct_round_shift(s5 - s7);
  x8 = s8 + s10;
  x9 = s9 + s11;
  x10 = s8 - s10;
  x11 = s9 - s11;
  x12 = fdct_round_shift(s12 + s14);
  x13 = fdct_round_shift(s13 + s15);
  x14 = fdct_round_shift(s12 - s14);
  x15 = fdct_round_shift(s13 - s15);

  // stage 4
  s2 = (- cospi_16_64) * (x2 + x3);
  s3 = cospi_16_64 * (x2 - x3);
  s6 = cospi_16_64 * (x6 + x7);
  s7 = cospi_16_64 * (- x6 + x7);
  s10 = cospi_16_64 * (x10 + x11);
  s11 = cospi_16_64 * (- x10 + x11);
  s14 = (- cospi_16_64) * (x14 + x15);
  s15 = cospi_16_64 * (x14 - x15);

  x2 = fdct_round_shift(s2);
  x3 = fdct_round_shift(s3);
  x6 = fdct_round_shift(s6);
  x7 = fdct_round_shift(s7);
  x10 = fdct_round_shift(s10);
  x11 = fdct_round_shift(s11);
  x14 = fdct_round_shift(s14);
  x15 = fdct_round_shift(s15);

  output[0] = x0;
  output[1] = - x8;
  output[2] = x12;
  output[3] = - x4;
  output[4] = x6;
  output[5] = x14;
  output[6] = x10;
  output[7] = x2;
  output[8] = x3;
  output[9] =  x11;
  output[10] = x15;
  output[11] = x7;
  output[12] = x5;
  output[13] = - x13;
  output[14] = x9;
  output[15] = - x1;
}

void vp9_fht16x16_c(const int16_t *input, tran_low_t *output,
                    int stride, int tx_type) {
  if (tx_type == DCT_DCT) {
    vp9_fdct16x16_c(input, output, stride);
  } else {
    tran_low_t out[256];
    tran_low_t *outptr = &out[0];
    int i, j;
    tran_low_t temp_in[16], temp_out[16];
    const transform_2d ht = FHT_16[tx_type];

#if CONFIG_EXT_TX
    int16_t flipped_input[16 * 16];
    maybe_flip_input(&input, &stride, 16, flipped_input, tx_type);
#endif

    // Columns
    for (i = 0; i < 16; ++i) {
      for (j = 0; j < 16; ++j)
        temp_in[j] = input[j * stride + i] * 4;
      ht.cols(temp_in, temp_out);
      for (j = 0; j < 16; ++j)
        outptr[j * 16 + i] = (temp_out[j] + 1 + (temp_out[j] < 0)) >> 2;
    }

    // Rows
    for (i = 0; i < 16; ++i) {
      for (j = 0; j < 16; ++j)
        temp_in[j] = out[j + i * 16];
      ht.rows(temp_in, temp_out);
      for (j = 0; j < 16; ++j)
        output[j + i * 16] = temp_out[j];
    }
  }
}

static INLINE tran_high_t dct_32_round(tran_high_t input) {
  tran_high_t rv = ROUND_POWER_OF_TWO(input, DCT_CONST_BITS);
  // TODO(debargha, peter.derivaz): Find new bounds for this assert,
  // and make the bounds consts.
  // assert(-131072 <= rv && rv <= 131071);
  return rv;
}

static INLINE tran_high_t half_round_shift(tran_high_t input) {
  tran_high_t rv = (input + 1 + (input < 0)) >> 2;
  return rv;
}

void vp9_fdct32(const tran_high_t *input, tran_high_t *output, int round) {
  tran_high_t step[32];
  // Stage 1
  step[0] = input[0] + input[(32 - 1)];
  step[1] = input[1] + input[(32 - 2)];
  step[2] = input[2] + input[(32 - 3)];
  step[3] = input[3] + input[(32 - 4)];
  step[4] = input[4] + input[(32 - 5)];
  step[5] = input[5] + input[(32 - 6)];
  step[6] = input[6] + input[(32 - 7)];
  step[7] = input[7] + input[(32 - 8)];
  step[8] = input[8] + input[(32 - 9)];
  step[9] = input[9] + input[(32 - 10)];
  step[10] = input[10] + input[(32 - 11)];
  step[11] = input[11] + input[(32 - 12)];
  step[12] = input[12] + input[(32 - 13)];
  step[13] = input[13] + input[(32 - 14)];
  step[14] = input[14] + input[(32 - 15)];
  step[15] = input[15] + input[(32 - 16)];
  step[16] = -input[16] + input[(32 - 17)];
  step[17] = -input[17] + input[(32 - 18)];
  step[18] = -input[18] + input[(32 - 19)];
  step[19] = -input[19] + input[(32 - 20)];
  step[20] = -input[20] + input[(32 - 21)];
  step[21] = -input[21] + input[(32 - 22)];
  step[22] = -input[22] + input[(32 - 23)];
  step[23] = -input[23] + input[(32 - 24)];
  step[24] = -input[24] + input[(32 - 25)];
  step[25] = -input[25] + input[(32 - 26)];
  step[26] = -input[26] + input[(32 - 27)];
  step[27] = -input[27] + input[(32 - 28)];
  step[28] = -input[28] + input[(32 - 29)];
  step[29] = -input[29] + input[(32 - 30)];
  step[30] = -input[30] + input[(32 - 31)];
  step[31] = -input[31] + input[(32 - 32)];

  // Stage 2
  output[0] = step[0] + step[16 - 1];
  output[1] = step[1] + step[16 - 2];
  output[2] = step[2] + step[16 - 3];
  output[3] = step[3] + step[16 - 4];
  output[4] = step[4] + step[16 - 5];
  output[5] = step[5] + step[16 - 6];
  output[6] = step[6] + step[16 - 7];
  output[7] = step[7] + step[16 - 8];
  output[8] = -step[8] + step[16 - 9];
  output[9] = -step[9] + step[16 - 10];
  output[10] = -step[10] + step[16 - 11];
  output[11] = -step[11] + step[16 - 12];
  output[12] = -step[12] + step[16 - 13];
  output[13] = -step[13] + step[16 - 14];
  output[14] = -step[14] + step[16 - 15];
  output[15] = -step[15] + step[16 - 16];

  output[16] = step[16];
  output[17] = step[17];
  output[18] = step[18];
  output[19] = step[19];

  output[20] = dct_32_round((-step[20] + step[27]) * cospi_16_64);
  output[21] = dct_32_round((-step[21] + step[26]) * cospi_16_64);
  output[22] = dct_32_round((-step[22] + step[25]) * cospi_16_64);
  output[23] = dct_32_round((-step[23] + step[24]) * cospi_16_64);

  output[24] = dct_32_round((step[24] + step[23]) * cospi_16_64);
  output[25] = dct_32_round((step[25] + step[22]) * cospi_16_64);
  output[26] = dct_32_round((step[26] + step[21]) * cospi_16_64);
  output[27] = dct_32_round((step[27] + step[20]) * cospi_16_64);

  output[28] = step[28];
  output[29] = step[29];
  output[30] = step[30];
  output[31] = step[31];

  // dump the magnitude by 4, hence the intermediate values are within
  // the range of 16 bits.
  if (round) {
    output[0] = half_round_shift(output[0]);
    output[1] = half_round_shift(output[1]);
    output[2] = half_round_shift(output[2]);
    output[3] = half_round_shift(output[3]);
    output[4] = half_round_shift(output[4]);
    output[5] = half_round_shift(output[5]);
    output[6] = half_round_shift(output[6]);
    output[7] = half_round_shift(output[7]);
    output[8] = half_round_shift(output[8]);
    output[9] = half_round_shift(output[9]);
    output[10] = half_round_shift(output[10]);
    output[11] = half_round_shift(output[11]);
    output[12] = half_round_shift(output[12]);
    output[13] = half_round_shift(output[13]);
    output[14] = half_round_shift(output[14]);
    output[15] = half_round_shift(output[15]);

    output[16] = half_round_shift(output[16]);
    output[17] = half_round_shift(output[17]);
    output[18] = half_round_shift(output[18]);
    output[19] = half_round_shift(output[19]);
    output[20] = half_round_shift(output[20]);
    output[21] = half_round_shift(output[21]);
    output[22] = half_round_shift(output[22]);
    output[23] = half_round_shift(output[23]);
    output[24] = half_round_shift(output[24]);
    output[25] = half_round_shift(output[25]);
    output[26] = half_round_shift(output[26]);
    output[27] = half_round_shift(output[27]);
    output[28] = half_round_shift(output[28]);
    output[29] = half_round_shift(output[29]);
    output[30] = half_round_shift(output[30]);
    output[31] = half_round_shift(output[31]);
  }

  // Stage 3
  step[0] = output[0] + output[(8 - 1)];
  step[1] = output[1] + output[(8 - 2)];
  step[2] = output[2] + output[(8 - 3)];
  step[3] = output[3] + output[(8 - 4)];
  step[4] = -output[4] + output[(8 - 5)];
  step[5] = -output[5] + output[(8 - 6)];
  step[6] = -output[6] + output[(8 - 7)];
  step[7] = -output[7] + output[(8 - 8)];
  step[8] = output[8];
  step[9] = output[9];
  step[10] = dct_32_round((-output[10] + output[13]) * cospi_16_64);
  step[11] = dct_32_round((-output[11] + output[12]) * cospi_16_64);
  step[12] = dct_32_round((output[12] + output[11]) * cospi_16_64);
  step[13] = dct_32_round((output[13] + output[10]) * cospi_16_64);
  step[14] = output[14];
  step[15] = output[15];

  step[16] = output[16] + output[23];
  step[17] = output[17] + output[22];
  step[18] = output[18] + output[21];
  step[19] = output[19] + output[20];
  step[20] = -output[20] + output[19];
  step[21] = -output[21] + output[18];
  step[22] = -output[22] + output[17];
  step[23] = -output[23] + output[16];
  step[24] = -output[24] + output[31];
  step[25] = -output[25] + output[30];
  step[26] = -output[26] + output[29];
  step[27] = -output[27] + output[28];
  step[28] = output[28] + output[27];
  step[29] = output[29] + output[26];
  step[30] = output[30] + output[25];
  step[31] = output[31] + output[24];

  // Stage 4
  output[0] = step[0] + step[3];
  output[1] = step[1] + step[2];
  output[2] = -step[2] + step[1];
  output[3] = -step[3] + step[0];
  output[4] = step[4];
  output[5] = dct_32_round((-step[5] + step[6]) * cospi_16_64);
  output[6] = dct_32_round((step[6] + step[5]) * cospi_16_64);
  output[7] = step[7];
  output[8] = step[8] + step[11];
  output[9] = step[9] + step[10];
  output[10] = -step[10] + step[9];
  output[11] = -step[11] + step[8];
  output[12] = -step[12] + step[15];
  output[13] = -step[13] + step[14];
  output[14] = step[14] + step[13];
  output[15] = step[15] + step[12];

  output[16] = step[16];
  output[17] = step[17];
  output[18] = dct_32_round(step[18] * -cospi_8_64 + step[29] * cospi_24_64);
  output[19] = dct_32_round(step[19] * -cospi_8_64 + step[28] * cospi_24_64);
  output[20] = dct_32_round(step[20] * -cospi_24_64 + step[27] * -cospi_8_64);
  output[21] = dct_32_round(step[21] * -cospi_24_64 + step[26] * -cospi_8_64);
  output[22] = step[22];
  output[23] = step[23];
  output[24] = step[24];
  output[25] = step[25];
  output[26] = dct_32_round(step[26] * cospi_24_64 + step[21] * -cospi_8_64);
  output[27] = dct_32_round(step[27] * cospi_24_64 + step[20] * -cospi_8_64);
  output[28] = dct_32_round(step[28] * cospi_8_64 + step[19] * cospi_24_64);
  output[29] = dct_32_round(step[29] * cospi_8_64 + step[18] * cospi_24_64);
  output[30] = step[30];
  output[31] = step[31];

  // Stage 5
  step[0] = dct_32_round((output[0] + output[1]) * cospi_16_64);
  step[1] = dct_32_round((-output[1] + output[0]) * cospi_16_64);
  step[2] = dct_32_round(output[2] * cospi_24_64 + output[3] * cospi_8_64);
  step[3] = dct_32_round(output[3] * cospi_24_64 - output[2] * cospi_8_64);
  step[4] = output[4] + output[5];
  step[5] = -output[5] + output[4];
  step[6] = -output[6] + output[7];
  step[7] = output[7] + output[6];
  step[8] = output[8];
  step[9] = dct_32_round(output[9] * -cospi_8_64 + output[14] * cospi_24_64);
  step[10] = dct_32_round(output[10] * -cospi_24_64 + output[13] * -cospi_8_64);
  step[11] = output[11];
  step[12] = output[12];
  step[13] = dct_32_round(output[13] * cospi_24_64 + output[10] * -cospi_8_64);
  step[14] = dct_32_round(output[14] * cospi_8_64 + output[9] * cospi_24_64);
  step[15] = output[15];

  step[16] = output[16] + output[19];
  step[17] = output[17] + output[18];
  step[18] = -output[18] + output[17];
  step[19] = -output[19] + output[16];
  step[20] = -output[20] + output[23];
  step[21] = -output[21] + output[22];
  step[22] = output[22] + output[21];
  step[23] = output[23] + output[20];
  step[24] = output[24] + output[27];
  step[25] = output[25] + output[26];
  step[26] = -output[26] + output[25];
  step[27] = -output[27] + output[24];
  step[28] = -output[28] + output[31];
  step[29] = -output[29] + output[30];
  step[30] = output[30] + output[29];
  step[31] = output[31] + output[28];

  // Stage 6
  output[0] = step[0];
  output[1] = step[1];
  output[2] = step[2];
  output[3] = step[3];
  output[4] = dct_32_round(step[4] * cospi_28_64 + step[7] * cospi_4_64);
  output[5] = dct_32_round(step[5] * cospi_12_64 + step[6] * cospi_20_64);
  output[6] = dct_32_round(step[6] * cospi_12_64 + step[5] * -cospi_20_64);
  output[7] = dct_32_round(step[7] * cospi_28_64 + step[4] * -cospi_4_64);
  output[8] = step[8] + step[9];
  output[9] = -step[9] + step[8];
  output[10] = -step[10] + step[11];
  output[11] = step[11] + step[10];
  output[12] = step[12] + step[13];
  output[13] = -step[13] + step[12];
  output[14] = -step[14] + step[15];
  output[15] = step[15] + step[14];

  output[16] = step[16];
  output[17] = dct_32_round(step[17] * -cospi_4_64 + step[30] * cospi_28_64);
  output[18] = dct_32_round(step[18] * -cospi_28_64 + step[29] * -cospi_4_64);
  output[19] = step[19];
  output[20] = step[20];
  output[21] = dct_32_round(step[21] * -cospi_20_64 + step[26] * cospi_12_64);
  output[22] = dct_32_round(step[22] * -cospi_12_64 + step[25] * -cospi_20_64);
  output[23] = step[23];
  output[24] = step[24];
  output[25] = dct_32_round(step[25] * cospi_12_64 + step[22] * -cospi_20_64);
  output[26] = dct_32_round(step[26] * cospi_20_64 + step[21] * cospi_12_64);
  output[27] = step[27];
  output[28] = step[28];
  output[29] = dct_32_round(step[29] * cospi_28_64 + step[18] * -cospi_4_64);
  output[30] = dct_32_round(step[30] * cospi_4_64 + step[17] * cospi_28_64);
  output[31] = step[31];

  // Stage 7
  step[0] = output[0];
  step[1] = output[1];
  step[2] = output[2];
  step[3] = output[3];
  step[4] = output[4];
  step[5] = output[5];
  step[6] = output[6];
  step[7] = output[7];
  step[8] = dct_32_round(output[8] * cospi_30_64 + output[15] * cospi_2_64);
  step[9] = dct_32_round(output[9] * cospi_14_64 + output[14] * cospi_18_64);
  step[10] = dct_32_round(output[10] * cospi_22_64 + output[13] * cospi_10_64);
  step[11] = dct_32_round(output[11] * cospi_6_64 + output[12] * cospi_26_64);
  step[12] = dct_32_round(output[12] * cospi_6_64 + output[11] * -cospi_26_64);
  step[13] = dct_32_round(output[13] * cospi_22_64 + output[10] * -cospi_10_64);
  step[14] = dct_32_round(output[14] * cospi_14_64 + output[9] * -cospi_18_64);
  step[15] = dct_32_round(output[15] * cospi_30_64 + output[8] * -cospi_2_64);

  step[16] = output[16] + output[17];
  step[17] = -output[17] + output[16];
  step[18] = -output[18] + output[19];
  step[19] = output[19] + output[18];
  step[20] = output[20] + output[21];
  step[21] = -output[21] + output[20];
  step[22] = -output[22] + output[23];
  step[23] = output[23] + output[22];
  step[24] = output[24] + output[25];
  step[25] = -output[25] + output[24];
  step[26] = -output[26] + output[27];
  step[27] = output[27] + output[26];
  step[28] = output[28] + output[29];
  step[29] = -output[29] + output[28];
  step[30] = -output[30] + output[31];
  step[31] = output[31] + output[30];

  // Final stage --- outputs indices are bit-reversed.
  output[0]  = step[0];
  output[16] = step[1];
  output[8]  = step[2];
  output[24] = step[3];
  output[4]  = step[4];
  output[20] = step[5];
  output[12] = step[6];
  output[28] = step[7];
  output[2]  = step[8];
  output[18] = step[9];
  output[10] = step[10];
  output[26] = step[11];
  output[6]  = step[12];
  output[22] = step[13];
  output[14] = step[14];
  output[30] = step[15];

  output[1]  = dct_32_round(step[16] * cospi_31_64 + step[31] * cospi_1_64);
  output[17] = dct_32_round(step[17] * cospi_15_64 + step[30] * cospi_17_64);
  output[9]  = dct_32_round(step[18] * cospi_23_64 + step[29] * cospi_9_64);
  output[25] = dct_32_round(step[19] * cospi_7_64 + step[28] * cospi_25_64);
  output[5]  = dct_32_round(step[20] * cospi_27_64 + step[27] * cospi_5_64);
  output[21] = dct_32_round(step[21] * cospi_11_64 + step[26] * cospi_21_64);
  output[13] = dct_32_round(step[22] * cospi_19_64 + step[25] * cospi_13_64);
  output[29] = dct_32_round(step[23] * cospi_3_64 + step[24] * cospi_29_64);
  output[3]  = dct_32_round(step[24] * cospi_3_64 + step[23] * -cospi_29_64);
  output[19] = dct_32_round(step[25] * cospi_19_64 + step[22] * -cospi_13_64);
  output[11] = dct_32_round(step[26] * cospi_11_64 + step[21] * -cospi_21_64);
  output[27] = dct_32_round(step[27] * cospi_27_64 + step[20] * -cospi_5_64);
  output[7]  = dct_32_round(step[28] * cospi_7_64 + step[19] * -cospi_25_64);
  output[23] = dct_32_round(step[29] * cospi_23_64 + step[18] * -cospi_9_64);
  output[15] = dct_32_round(step[30] * cospi_15_64 + step[17] * -cospi_17_64);
  output[31] = dct_32_round(step[31] * cospi_31_64 + step[16] * -cospi_1_64);
}

void vp9_fdct32x32_1_c(const int16_t *input, tran_low_t *output, int stride) {
  int r, c;
  tran_low_t sum = 0;
  for (r = 0; r < 32; ++r)
    for (c = 0; c < 32; ++c)
      sum += input[r * stride + c];

  output[0] = sum >> 3;
  output[1] = 0;
}

void vp9_fdct32x32_c(const int16_t *input, tran_low_t *out, int stride) {
  int i, j;
  tran_high_t output[32 * 32];

  // Columns
  for (i = 0; i < 32; ++i) {
    tran_high_t temp_in[32], temp_out[32];
    for (j = 0; j < 32; ++j)
      temp_in[j] = input[j * stride + i] * 4;
    vp9_fdct32(temp_in, temp_out, 0);
    for (j = 0; j < 32; ++j)
      output[j * 32 + i] = (temp_out[j] + 1 + (temp_out[j] > 0)) >> 2;
  }

  // Rows
  for (i = 0; i < 32; ++i) {
    tran_high_t temp_in[32], temp_out[32];
    for (j = 0; j < 32; ++j)
      temp_in[j] = output[j + i * 32];
    vp9_fdct32(temp_in, temp_out, 0);
    for (j = 0; j < 32; ++j)
      out[j + i * 32] = (tran_low_t)
          ((temp_out[j] + 1 + (temp_out[j] < 0)) >> 2);
  }
}

#if CONFIG_WAVELETS
void vp9_fdct32x32_noscale_c(const int16_t *input, tran_low_t *out,
                             int stride) {
  int i, j;
  tran_high_t output[32 * 32];

  // Columns
  for (i = 0; i < 32; ++i) {
    tran_high_t temp_in[32], temp_out[32];
    for (j = 0; j < 32; ++j)
      temp_in[j] = input[j * stride + i];
    vp9_fdct32(temp_in, temp_out, 0);
    for (j = 0; j < 32; ++j)
      output[j * 32 + i] = (temp_out[j] + 1 + (temp_out[j] > 0)) >> 2;
  }

  // Rows
  for (i = 0; i < 32; ++i) {
    tran_high_t temp_in[32], temp_out[32];
    for (j = 0; j < 32; ++j)
      temp_in[j] = output[j + i * 32];
    vp9_fdct32(temp_in, temp_out, 0);
    for (j = 0; j < 32; ++j)
      out[j + i * 32] = (tran_low_t)
          ((temp_out[j] + 1 + (temp_out[j] < 0)) >> 2);
  }
}
#endif  // CONFIG_WAVELETS

// Note that although we use dct_32_round in dct32 computation flow,
// this 2d fdct32x32 for rate-distortion optimization loop is operating
// within 16 bits precision.
void vp9_fdct32x32_rd_c(const int16_t *input, tran_low_t *out, int stride) {
  int i, j;
  tran_high_t output[32 * 32];

  // Columns
  for (i = 0; i < 32; ++i) {
    tran_high_t temp_in[32], temp_out[32];
    for (j = 0; j < 32; ++j)
      temp_in[j] = input[j * stride + i] * 4;
    vp9_fdct32(temp_in, temp_out, 0);
    for (j = 0; j < 32; ++j)
      // TODO(cd): see quality impact of only doing
      //           output[j * 32 + i] = (temp_out[j] + 1) >> 2;
      //           PS: also change code in vp9/encoder/x86/vp9_dct_sse2.c
      output[j * 32 + i] = (temp_out[j] + 1 + (temp_out[j] > 0)) >> 2;
  }

  // Rows
  for (i = 0; i < 32; ++i) {
    tran_high_t temp_in[32], temp_out[32];
    for (j = 0; j < 32; ++j)
      temp_in[j] = output[j + i * 32];
    vp9_fdct32(temp_in, temp_out, 1);
    for (j = 0; j < 32; ++j)
      out[j + i * 32] = (tran_low_t)temp_out[j];
  }
}

#if CONFIG_TX_SKIP
void vp9_tx_identity_rect(const int16_t *input, tran_low_t *out,
                          int row, int col,
                          int stride_in, int stride_out, int shift) {
  int r, c;
  for (r = 0; r < row; r++)
    for (c = 0; c < col; c++) {
      out[stride_out * r + c] = input[stride_in * r + c] * (1 << shift);
    }
}

void vp9_tx_identity(const int16_t *input, tran_low_t *out, int stride,
                     int bs, int shift) {
  vp9_tx_identity_rect(input, out, bs, bs, stride, bs, shift);
}
#endif

#if CONFIG_TX64X64
// TODO(debargha): Using a floating point implementation for now.
// Should re-use the 32x32 integer dct we already have.
static void dct32_1d(double *input, double *output, int stride) {
  static const double C1 = 0.998795456205;  // cos(pi * 1 / 64)
  static const double C2 = 0.995184726672;  // cos(pi * 2 / 64)
  static const double C3 = 0.989176509965;  // cos(pi * 3 / 64)
  static const double C4 = 0.980785280403;  // cos(pi * 4 / 64)
  static const double C5 = 0.970031253195;  // cos(pi * 5 / 64)
  static const double C6 = 0.956940335732;  // cos(pi * 6 / 64)
  static const double C7 = 0.941544065183;  // cos(pi * 7 / 64)
  static const double C8 = 0.923879532511;  // cos(pi * 8 / 64)
  static const double C9 = 0.903989293123;  // cos(pi * 9 / 64)
  static const double C10 = 0.881921264348;  // cos(pi * 10 / 64)
  static const double C11 = 0.857728610000;  // cos(pi * 11 / 64)
  static const double C12 = 0.831469612303;  // cos(pi * 12 / 64)
  static const double C13 = 0.803207531481;  // cos(pi * 13 / 64)
  static const double C14 = 0.773010453363;  // cos(pi * 14 / 64)
  static const double C15 = 0.740951125355;  // cos(pi * 15 / 64)
  static const double C16 = 0.707106781187;  // cos(pi * 16 / 64)
  static const double C17 = 0.671558954847;  // cos(pi * 17 / 64)
  static const double C18 = 0.634393284164;  // cos(pi * 18 / 64)
  static const double C19 = 0.595699304492;  // cos(pi * 19 / 64)
  static const double C20 = 0.555570233020;  // cos(pi * 20 / 64)
  static const double C21 = 0.514102744193;  // cos(pi * 21 / 64)
  static const double C22 = 0.471396736826;  // cos(pi * 22 / 64)
  static const double C23 = 0.427555093430;  // cos(pi * 23 / 64)
  static const double C24 = 0.382683432365;  // cos(pi * 24 / 64)
  static const double C25 = 0.336889853392;  // cos(pi * 25 / 64)
  static const double C26 = 0.290284677254;  // cos(pi * 26 / 64)
  static const double C27 = 0.242980179903;  // cos(pi * 27 / 64)
  static const double C28 = 0.195090322016;  // cos(pi * 28 / 64)
  static const double C29 = 0.146730474455;  // cos(pi * 29 / 64)
  static const double C30 = 0.098017140330;  // cos(pi * 30 / 64)
  static const double C31 = 0.049067674327;  // cos(pi * 31 / 64)

  double step[32];

  // Stage 1
  step[0] = input[stride*0] + input[stride*(32 - 1)];
  step[1] = input[stride*1] + input[stride*(32 - 2)];
  step[2] = input[stride*2] + input[stride*(32 - 3)];
  step[3] = input[stride*3] + input[stride*(32 - 4)];
  step[4] = input[stride*4] + input[stride*(32 - 5)];
  step[5] = input[stride*5] + input[stride*(32 - 6)];
  step[6] = input[stride*6] + input[stride*(32 - 7)];
  step[7] = input[stride*7] + input[stride*(32 - 8)];
  step[8] = input[stride*8] + input[stride*(32 - 9)];
  step[9] = input[stride*9] + input[stride*(32 - 10)];
  step[10] = input[stride*10] + input[stride*(32 - 11)];
  step[11] = input[stride*11] + input[stride*(32 - 12)];
  step[12] = input[stride*12] + input[stride*(32 - 13)];
  step[13] = input[stride*13] + input[stride*(32 - 14)];
  step[14] = input[stride*14] + input[stride*(32 - 15)];
  step[15] = input[stride*15] + input[stride*(32 - 16)];
  step[16] = -input[stride*16] + input[stride*(32 - 17)];
  step[17] = -input[stride*17] + input[stride*(32 - 18)];
  step[18] = -input[stride*18] + input[stride*(32 - 19)];
  step[19] = -input[stride*19] + input[stride*(32 - 20)];
  step[20] = -input[stride*20] + input[stride*(32 - 21)];
  step[21] = -input[stride*21] + input[stride*(32 - 22)];
  step[22] = -input[stride*22] + input[stride*(32 - 23)];
  step[23] = -input[stride*23] + input[stride*(32 - 24)];
  step[24] = -input[stride*24] + input[stride*(32 - 25)];
  step[25] = -input[stride*25] + input[stride*(32 - 26)];
  step[26] = -input[stride*26] + input[stride*(32 - 27)];
  step[27] = -input[stride*27] + input[stride*(32 - 28)];
  step[28] = -input[stride*28] + input[stride*(32 - 29)];
  step[29] = -input[stride*29] + input[stride*(32 - 30)];
  step[30] = -input[stride*30] + input[stride*(32 - 31)];
  step[31] = -input[stride*31] + input[stride*(32 - 32)];

  // Stage 2
  output[stride*0] = step[0] + step[16 - 1];
  output[stride*1] = step[1] + step[16 - 2];
  output[stride*2] = step[2] + step[16 - 3];
  output[stride*3] = step[3] + step[16 - 4];
  output[stride*4] = step[4] + step[16 - 5];
  output[stride*5] = step[5] + step[16 - 6];
  output[stride*6] = step[6] + step[16 - 7];
  output[stride*7] = step[7] + step[16 - 8];
  output[stride*8] = -step[8] + step[16 - 9];
  output[stride*9] = -step[9] + step[16 - 10];
  output[stride*10] = -step[10] + step[16 - 11];
  output[stride*11] = -step[11] + step[16 - 12];
  output[stride*12] = -step[12] + step[16 - 13];
  output[stride*13] = -step[13] + step[16 - 14];
  output[stride*14] = -step[14] + step[16 - 15];
  output[stride*15] = -step[15] + step[16 - 16];

  output[stride*16] = step[16];
  output[stride*17] = step[17];
  output[stride*18] = step[18];
  output[stride*19] = step[19];

  output[stride*20] = (-step[20] + step[27])*C16;
  output[stride*21] = (-step[21] + step[26])*C16;
  output[stride*22] = (-step[22] + step[25])*C16;
  output[stride*23] = (-step[23] + step[24])*C16;

  output[stride*24] = (step[24] + step[23])*C16;
  output[stride*25] = (step[25] + step[22])*C16;
  output[stride*26] = (step[26] + step[21])*C16;
  output[stride*27] = (step[27] + step[20])*C16;

  output[stride*28] = step[28];
  output[stride*29] = step[29];
  output[stride*30] = step[30];
  output[stride*31] = step[31];

  // Stage 3
  step[0] = output[stride*0] + output[stride*(8 - 1)];
  step[1] = output[stride*1] + output[stride*(8 - 2)];
  step[2] = output[stride*2] + output[stride*(8 - 3)];
  step[3] = output[stride*3] + output[stride*(8 - 4)];
  step[4] = -output[stride*4] + output[stride*(8 - 5)];
  step[5] = -output[stride*5] + output[stride*(8 - 6)];
  step[6] = -output[stride*6] + output[stride*(8 - 7)];
  step[7] = -output[stride*7] + output[stride*(8 - 8)];
  step[8] = output[stride*8];
  step[9] = output[stride*9];
  step[10] = (-output[stride*10] + output[stride*13])*C16;
  step[11] = (-output[stride*11] + output[stride*12])*C16;
  step[12] = (output[stride*12] + output[stride*11])*C16;
  step[13] = (output[stride*13] + output[stride*10])*C16;
  step[14] = output[stride*14];
  step[15] = output[stride*15];

  step[16] = output[stride*16] + output[stride*23];
  step[17] = output[stride*17] + output[stride*22];
  step[18] = output[stride*18] + output[stride*21];
  step[19] = output[stride*19] + output[stride*20];
  step[20] = -output[stride*20] + output[stride*19];
  step[21] = -output[stride*21] + output[stride*18];
  step[22] = -output[stride*22] + output[stride*17];
  step[23] = -output[stride*23] + output[stride*16];
  step[24] = -output[stride*24] + output[stride*31];
  step[25] = -output[stride*25] + output[stride*30];
  step[26] = -output[stride*26] + output[stride*29];
  step[27] = -output[stride*27] + output[stride*28];
  step[28] = output[stride*28] + output[stride*27];
  step[29] = output[stride*29] + output[stride*26];
  step[30] = output[stride*30] + output[stride*25];
  step[31] = output[stride*31] + output[stride*24];

  // Stage 4
  output[stride*0] = step[0] + step[3];
  output[stride*1] = step[1] + step[2];
  output[stride*2] = -step[2] + step[1];
  output[stride*3] = -step[3] + step[0];
  output[stride*4] = step[4];
  output[stride*5] = (-step[5] + step[6])*C16;
  output[stride*6] = (step[6] + step[5])*C16;
  output[stride*7] = step[7];
  output[stride*8] = step[8] + step[11];
  output[stride*9] = step[9] + step[10];
  output[stride*10] = -step[10] + step[9];
  output[stride*11] = -step[11] + step[8];
  output[stride*12] = -step[12] + step[15];
  output[stride*13] = -step[13] + step[14];
  output[stride*14] = step[14] + step[13];
  output[stride*15] = step[15] + step[12];

  output[stride*16] = step[16];
  output[stride*17] = step[17];
  output[stride*18] = step[18]*-C8 + step[29]*C24;
  output[stride*19] = step[19]*-C8 + step[28]*C24;
  output[stride*20] = step[20]*-C24 + step[27]*-C8;
  output[stride*21] = step[21]*-C24 + step[26]*-C8;
  output[stride*22] = step[22];
  output[stride*23] = step[23];
  output[stride*24] = step[24];
  output[stride*25] = step[25];
  output[stride*26] = step[26]*C24 + step[21]*-C8;
  output[stride*27] = step[27]*C24 + step[20]*-C8;
  output[stride*28] = step[28]*C8 + step[19]*C24;
  output[stride*29] = step[29]*C8 + step[18]*C24;
  output[stride*30] = step[30];
  output[stride*31] = step[31];

  // Stage 5
  step[0] = (output[stride*0] + output[stride*1]) * C16;
  step[1] = (-output[stride*1] + output[stride*0]) * C16;
  step[2] = output[stride*2]*C24 + output[stride*3] * C8;
  step[3] = output[stride*3]*C24 - output[stride*2] * C8;
  step[4] = output[stride*4] + output[stride*5];
  step[5] = -output[stride*5] + output[stride*4];
  step[6] = -output[stride*6] + output[stride*7];
  step[7] = output[stride*7] + output[stride*6];
  step[8] = output[stride*8];
  step[9] = output[stride*9]*-C8 + output[stride*14]*C24;
  step[10] = output[stride*10]*-C24 + output[stride*13]*-C8;
  step[11] = output[stride*11];
  step[12] = output[stride*12];
  step[13] = output[stride*13]*C24 + output[stride*10]*-C8;
  step[14] = output[stride*14]*C8 + output[stride*9]*C24;
  step[15] = output[stride*15];

  step[16] = output[stride*16] + output[stride*19];
  step[17] = output[stride*17] + output[stride*18];
  step[18] = -output[stride*18] + output[stride*17];
  step[19] = -output[stride*19] + output[stride*16];
  step[20] = -output[stride*20] + output[stride*23];
  step[21] = -output[stride*21] + output[stride*22];
  step[22] = output[stride*22] + output[stride*21];
  step[23] = output[stride*23] + output[stride*20];
  step[24] = output[stride*24] + output[stride*27];
  step[25] = output[stride*25] + output[stride*26];
  step[26] = -output[stride*26] + output[stride*25];
  step[27] = -output[stride*27] + output[stride*24];
  step[28] = -output[stride*28] + output[stride*31];
  step[29] = -output[stride*29] + output[stride*30];
  step[30] = output[stride*30] + output[stride*29];
  step[31] = output[stride*31] + output[stride*28];

  // Stage 6
  output[stride*0] = step[0];
  output[stride*1] = step[1];
  output[stride*2] = step[2];
  output[stride*3] = step[3];
  output[stride*4] = step[4]*C28 + step[7]*C4;
  output[stride*5] = step[5]*C12 + step[6]*C20;
  output[stride*6] = step[6]*C12 + step[5]*-C20;
  output[stride*7] = step[7]*C28 + step[4]*-C4;
  output[stride*8] = step[8] + step[9];
  output[stride*9] = -step[9] + step[8];
  output[stride*10] = -step[10] + step[11];
  output[stride*11] = step[11] + step[10];
  output[stride*12] = step[12] + step[13];
  output[stride*13] = -step[13] + step[12];
  output[stride*14] = -step[14] + step[15];
  output[stride*15] = step[15] + step[14];

  output[stride*16] = step[16];
  output[stride*17] = step[17]*-C4 + step[30]*C28;
  output[stride*18] = step[18]*-C28 + step[29]*-C4;
  output[stride*19] = step[19];
  output[stride*20] = step[20];
  output[stride*21] = step[21]*-C20 + step[26]*C12;
  output[stride*22] = step[22]*-C12 + step[25]*-C20;
  output[stride*23] = step[23];
  output[stride*24] = step[24];
  output[stride*25] = step[25]*C12 + step[22]*-C20;
  output[stride*26] = step[26]*C20 + step[21]*C12;
  output[stride*27] = step[27];
  output[stride*28] = step[28];
  output[stride*29] = step[29]*C28 + step[18]*-C4;
  output[stride*30] = step[30]*C4 + step[17]*C28;
  output[stride*31] = step[31];

  // Stage 7
  step[0] = output[stride*0];
  step[1] = output[stride*1];
  step[2] = output[stride*2];
  step[3] = output[stride*3];
  step[4] = output[stride*4];
  step[5] = output[stride*5];
  step[6] = output[stride*6];
  step[7] = output[stride*7];
  step[8] = output[stride*8]*C30 + output[stride*15]*C2;
  step[9] = output[stride*9]*C14 + output[stride*14]*C18;
  step[10] = output[stride*10]*C22 + output[stride*13]*C10;
  step[11] = output[stride*11]*C6 + output[stride*12]*C26;
  step[12] = output[stride*12]*C6 + output[stride*11]*-C26;
  step[13] = output[stride*13]*C22 + output[stride*10]*-C10;
  step[14] = output[stride*14]*C14 + output[stride*9]*-C18;
  step[15] = output[stride*15]*C30 + output[stride*8]*-C2;

  step[16] = output[stride*16] + output[stride*17];
  step[17] = -output[stride*17] + output[stride*16];
  step[18] = -output[stride*18] + output[stride*19];
  step[19] = output[stride*19] + output[stride*18];
  step[20] = output[stride*20] + output[stride*21];
  step[21] = -output[stride*21] + output[stride*20];
  step[22] = -output[stride*22] + output[stride*23];
  step[23] = output[stride*23] + output[stride*22];
  step[24] = output[stride*24] + output[stride*25];
  step[25] = -output[stride*25] + output[stride*24];
  step[26] = -output[stride*26] + output[stride*27];
  step[27] = output[stride*27] + output[stride*26];
  step[28] = output[stride*28] + output[stride*29];
  step[29] = -output[stride*29] + output[stride*28];
  step[30] = -output[stride*30] + output[stride*31];
  step[31] = output[stride*31] + output[stride*30];

  // Final stage --- outputs indices are bit-reversed.
  output[stride*0] = step[0];
  output[stride*16] = step[1];
  output[stride*8] = step[2];
  output[stride*24] = step[3];
  output[stride*4] = step[4];
  output[stride*20] = step[5];
  output[stride*12] = step[6];
  output[stride*28] = step[7];
  output[stride*2] = step[8];
  output[stride*18] = step[9];
  output[stride*10] = step[10];
  output[stride*26] = step[11];
  output[stride*6] = step[12];
  output[stride*22] = step[13];
  output[stride*14] = step[14];
  output[stride*30] = step[15];

  output[stride*1] = step[16]*C31 + step[31]*C1;
  output[stride*17] = step[17]*C15 + step[30]*C17;
  output[stride*9] = step[18]*C23 + step[29]*C9;
  output[stride*25] = step[19]*C7 + step[28]*C25;
  output[stride*5] = step[20]*C27 + step[27]*C5;
  output[stride*21] = step[21]*C11 + step[26]*C21;
  output[stride*13] = step[22]*C19 + step[25]*C13;
  output[stride*29] = step[23]*C3 + step[24]*C29;
  output[stride*3] = step[24]*C3 + step[23]*-C29;
  output[stride*19] = step[25]*C19 + step[22]*-C13;
  output[stride*11] = step[26]*C11 + step[21]*-C21;
  output[stride*27] = step[27]*C27 + step[20]*-C5;
  output[stride*7] = step[28]*C7 + step[19]*-C25;
  output[stride*23] = step[29]*C23 + step[18]*-C9;
  output[stride*15] = step[30]*C15 + step[17]*-C17;
  output[stride*31] = step[31]*C31 + step[16]*-C1;
}

static void dct64_1d(double *input, double *output, int stride) {
  double step1[64], step2[64];
  int i;
  static const double C[64] = {
    1.00000000000000000000,  // cos(0 * pi / 128)
    0.99969881869620424997,  // cos(1 * pi / 128)
    0.99879545620517240501,  // cos(2 * pi / 128)
    0.99729045667869020697,  // cos(3 * pi / 128)
    0.99518472667219692873,  // cos(4 * pi / 128)
    0.99247953459870996706,  // cos(5 * pi / 128)
    0.98917650996478101444,  // cos(6 * pi / 128)
    0.98527764238894122162,  // cos(7 * pi / 128)
    0.98078528040323043058,  // cos(8 * pi / 128)
    0.97570213003852857003,  // cos(9 * pi / 128)
    0.97003125319454397424,  // cos(10 * pi / 128)
    0.96377606579543984022,  // cos(11 * pi / 128)
    0.95694033573220882438,  // cos(12 * pi / 128)
    0.94952818059303667475,  // cos(13 * pi / 128)
    0.94154406518302080631,  // cos(14 * pi / 128)
    0.93299279883473895669,  // cos(15 * pi / 128)
    0.92387953251128673848,  // cos(16 * pi / 128)
    0.91420975570353069095,  // cos(17 * pi / 128)
    0.90398929312344333820,  // cos(18 * pi / 128)
    0.89322430119551532446,  // cos(19 * pi / 128)
    0.88192126434835504956,  // cos(20 * pi / 128)
    0.87008699110871146054,  // cos(21 * pi / 128)
    0.85772861000027211809,  // cos(22 * pi / 128)
    0.84485356524970711689,  // cos(23 * pi / 128)
    0.83146961230254523567,  // cos(24 * pi / 128)
    0.81758481315158371139,  // cos(25 * pi / 128)
    0.80320753148064494287,  // cos(26 * pi / 128)
    0.78834642762660633863,  // cos(27 * pi / 128)
    0.77301045336273699338,  // cos(28 * pi / 128)
    0.75720884650648456748,  // cos(29 * pi / 128)
    0.74095112535495921691,  // cos(30 * pi / 128)
    0.72424708295146700276,  // cos(31 * pi / 128)
    0.70710678118654757274,  // cos(32 * pi / 128)
    0.68954054473706694051,  // cos(33 * pi / 128)
    0.67155895484701844111,  // cos(34 * pi / 128)
    0.65317284295377686654,  // cos(35 * pi / 128)
    0.63439328416364559882,  // cos(36 * pi / 128)
    0.61523159058062693028,  // cos(37 * pi / 128)
    0.59569930449243346793,  // cos(38 * pi / 128)
    0.57580819141784544968,  // cos(39 * pi / 128)
    0.55557023301960228867,  // cos(40 * pi / 128)
    0.53499761988709737537,  // cos(41 * pi / 128)
    0.51410274419322177231,  // cos(42 * pi / 128)
    0.49289819222978414892,  // cos(43 * pi / 128)
    0.47139673682599780857,  // cos(44 * pi / 128)
    0.44961132965460659516,  // cos(45 * pi / 128)
    0.42755509343028219593,  // cos(46 * pi / 128)
    0.40524131400498980549,  // cos(47 * pi / 128)
    0.38268343236508983729,  // cos(48 * pi / 128)
    0.35989503653498827740,  // cos(49 * pi / 128)
    0.33688985339222005111,  // cos(50 * pi / 128)
    0.31368174039889151761,  // cos(51 * pi / 128)
    0.29028467725446227554,  // cos(52 * pi / 128)
    0.26671275747489842090,  // cos(53 * pi / 128)
    0.24298017990326398197,  // cos(54 * pi / 128)
    0.21910124015686976984,  // cos(55 * pi / 128)
    0.19509032201612830359,  // cos(56 * pi / 128)
    0.17096188876030135595,  // cos(57 * pi / 128)
    0.14673047445536174793,  // cos(58 * pi / 128)
    0.12241067519921627893,  // cos(59 * pi / 128)
    0.09801714032956077016,  // cos(60 * pi / 128)
    0.07356456359966745406,  // cos(61 * pi / 128)
    0.04906767432741813290,  // cos(62 * pi / 128)
    0.02454122852291226731,  // cos(63 * pi / 128)
  };

  for (i = 0; i < 32; ++i) {
    step1[i] = input[stride * i] + input[stride * (63 - i)];
    step1[32 + i] = (input[stride * i] -
                     input[stride * (63 - i)]) * C[i * 2 + 1];
  }

  dct32_1d(step1, step2, 1);
  dct32_1d(step1 + 32, step2 + 32, 1);

  for (i = 0; i < 64; i += 2) {
    output[stride*i] = step2[i / 2];
  }
  output[stride * 1] = 2 * step2[32] * C[32];
  for (i = 3; i < 64; i += 2) {
    output[stride * i] = 2 * step2[32 + i / 2] - output[stride * (i - 2)];
  }
}

void vp9_fdct64x64_c(const int16_t *input, tran_low_t *out, int stride) {
  // vp9_clear_system_state();  // Make it simd safe : __asm emms;
  {
    int i, j;
    double output[4096];
    // First transform columns
    for (i = 0; i < 64; i++) {
      double temp_in[64], temp_out[64];
      for (j = 0; j < 64; j++)
        temp_in[j] = input[j * stride + i];
      dct64_1d(temp_in, temp_out, 1);
      for (j = 0; j < 64; j++)
        output[j * 64 + i] = temp_out[j];
    }
    // Then transform rows
    for (i = 0; i < 64; ++i) {
      double temp_in[64], temp_out[64];
      for (j = 0; j < 64; ++j)
        temp_in[j] = output[j + i * 64];
      dct64_1d(temp_in, temp_out, 1);
      for (j = 0; j < 64; ++j)
        output[j + i * 64] = temp_out[j];
    }
    // Scale by some magic number
    for (i = 0; i < 4096; i++) {
      out[i] = (tran_low_t)round(output[i] / 16);
    }
  }
  // vp9_clear_system_state();  // Make it simd safe : __asm emms;
}

void vp9_fdct64x64_1_c(const int16_t *input, tran_low_t *output, int stride) {
  int r, c;
  tran_low_t sum = 0;
  for (r = 0; r < 64; ++r)
    for (c = 0; c < 64; ++c)
      sum += input[r * stride + c];

  output[0] = sum >> 5;
  output[1] = 0;
}
#endif

#if CONFIG_VP9_HIGHBITDEPTH
void vp9_highbd_fdct4x4_c(const int16_t *input, tran_low_t *output,
                          int stride) {
  vp9_fdct4x4_c(input, output, stride);
}

void vp9_highbd_fht4x4_c(const int16_t *input, tran_low_t *output,
                         int stride, int tx_type) {
  vp9_fht4x4_c(input, output, stride, tx_type);
}

void vp9_highbd_fdct8x8_1_c(const int16_t *input, tran_low_t *final_output,
                            int stride) {
  vp9_fdct8x8_1_c(input, final_output, stride);
}

void vp9_highbd_fdct8x8_c(const int16_t *input, tran_low_t *final_output,
                          int stride) {
  vp9_fdct8x8_c(input, final_output, stride);
}

void vp9_highbd_fdct16x16_1_c(const int16_t *input, tran_low_t *output,
                              int stride) {
  vp9_fdct16x16_1_c(input, output, stride);
}

void vp9_highbd_fdct16x16_c(const int16_t *input, tran_low_t *output,
                            int stride) {
  vp9_fdct16x16_c(input, output, stride);
}

void vp9_highbd_fht8x8_c(const int16_t *input, tran_low_t *output,
                         int stride, int tx_type) {
  vp9_fht8x8_c(input, output, stride, tx_type);
}

void vp9_highbd_fwht4x4_c(const int16_t *input, tran_low_t *output,
                          int stride) {
  vp9_fwht4x4_c(input, output, stride);
}

void vp9_highbd_fht16x16_c(const int16_t *input, tran_low_t *output,
                           int stride, int tx_type) {
  vp9_fht16x16_c(input, output, stride, tx_type);
}

void vp9_highbd_fdct32x32_1_c(const int16_t *input, tran_low_t *out,
                              int stride) {
  vp9_fdct32x32_1_c(input, out, stride);
}

void vp9_highbd_fdct32x32_c(const int16_t *input, tran_low_t *out, int stride) {
  vp9_fdct32x32_c(input, out, stride);
}

void vp9_highbd_fdct32x32_rd_c(const int16_t *input, tran_low_t *out,
                               int stride) {
  vp9_fdct32x32_rd_c(input, out, stride);
}

#if CONFIG_TX64X64
void vp9_highbd_fdct64x64_1_c(const int16_t *input, tran_low_t *out,
                              int stride) {
  vp9_fdct64x64_1_c(input, out, stride);
}

void vp9_highbd_fdct64x64_c(const int16_t *input, tran_low_t *out, int stride) {
  vp9_fdct64x64_c(input, out, stride);
}
#endif  // CONFIG_TX64X64
#endif  // CONFIG_VP9_HIGHBITDEPTH
