/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_dsp/fwd_txfm.h"

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
      out[0] = (tran_low_t)fdct_round_shift(temp1);
      out[2] = (tran_low_t)fdct_round_shift(temp2);
      temp1 = step[2] * cospi_24_64 + step[3] * cospi_8_64;
      temp2 = -step[2] * cospi_8_64 + step[3] * cospi_24_64;
      out[1] = (tran_low_t)fdct_round_shift(temp1);
      out[3] = (tran_low_t)fdct_round_shift(temp2);
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

void vp9_fdct8x8_c(const int16_t *input, tran_low_t *final_output, int stride) {
  int i, j;
  tran_low_t intermediate[64];
  int pass;
  tran_low_t *output = intermediate;
  const tran_low_t *in = NULL;

  // Transform columns
  for (pass = 0; pass < 2; ++pass) {
    tran_high_t s0, s1, s2, s3, s4, s5, s6, s7;  // canbe16
    tran_high_t t0, t1, t2, t3;                  // needs32
    tran_high_t x0, x1, x2, x3;                  // canbe16

    int i;
    for (i = 0; i < 8; i++) {
      // stage 1
      if (pass == 0) {
        s0 = (input[0 * stride] + input[7 * stride]) * 4;
        s1 = (input[1 * stride] + input[6 * stride]) * 4;
        s2 = (input[2 * stride] + input[5 * stride]) * 4;
        s3 = (input[3 * stride] + input[4 * stride]) * 4;
        s4 = (input[3 * stride] - input[4 * stride]) * 4;
        s5 = (input[2 * stride] - input[5 * stride]) * 4;
        s6 = (input[1 * stride] - input[6 * stride]) * 4;
        s7 = (input[0 * stride] - input[7 * stride]) * 4;
        ++input;
      } else {
        s0 = in[0 * 8] + in[7 * 8];
        s1 = in[1 * 8] + in[6 * 8];
        s2 = in[2 * 8] + in[5 * 8];
        s3 = in[3 * 8] + in[4 * 8];
        s4 = in[3 * 8] - in[4 * 8];
        s5 = in[2 * 8] - in[5 * 8];
        s6 = in[1 * 8] - in[6 * 8];
        s7 = in[0 * 8] - in[7 * 8];
        ++in;
      }

      // fdct4(step, step);
      x0 = s0 + s3;
      x1 = s1 + s2;
      x2 = s1 - s2;
      x3 = s0 - s3;
      t0 = (x0 + x1) * cospi_16_64;
      t1 = (x0 - x1) * cospi_16_64;
      t2 =  x2 * cospi_24_64 + x3 *  cospi_8_64;
      t3 = -x2 * cospi_8_64  + x3 * cospi_24_64;
      output[0] = (tran_low_t)fdct_round_shift(t0);
      output[2] = (tran_low_t)fdct_round_shift(t2);
      output[4] = (tran_low_t)fdct_round_shift(t1);
      output[6] = (tran_low_t)fdct_round_shift(t3);

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
      output[1] = (tran_low_t)fdct_round_shift(t0);
      output[3] = (tran_low_t)fdct_round_shift(t2);
      output[5] = (tran_low_t)fdct_round_shift(t1);
      output[7] = (tran_low_t)fdct_round_shift(t3);
      output += 8;
    }
    in  = intermediate;
    output = final_output;
  }

  // Rows
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      final_output[j + i * 8] /= 2;
  }
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
        out[0] = (tran_low_t)fdct_round_shift(t0);
        out[4] = (tran_low_t)fdct_round_shift(t2);
        out[8] = (tran_low_t)fdct_round_shift(t1);
        out[12] = (tran_low_t)fdct_round_shift(t3);

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
        out[2] = (tran_low_t)fdct_round_shift(t0);
        out[6] = (tran_low_t)fdct_round_shift(t2);
        out[10] = (tran_low_t)fdct_round_shift(t1);
        out[14] = (tran_low_t)fdct_round_shift(t3);
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
        out[1] = (tran_low_t)fdct_round_shift(temp1);
        out[9] = (tran_low_t)fdct_round_shift(temp2);
        temp1 = step1[2] * cospi_22_64 + step1[5] * cospi_10_64;
        temp2 = step1[3] *  cospi_6_64 + step1[4] * cospi_26_64;
        out[5] = (tran_low_t)fdct_round_shift(temp1);
        out[13] = (tran_low_t)fdct_round_shift(temp2);
        temp1 = step1[3] * -cospi_26_64 + step1[4] *  cospi_6_64;
        temp2 = step1[2] * -cospi_10_64 + step1[5] * cospi_22_64;
        out[3] = (tran_low_t)fdct_round_shift(temp1);
        out[11] = (tran_low_t)fdct_round_shift(temp2);
        temp1 = step1[1] * -cospi_18_64 + step1[6] * cospi_14_64;
        temp2 = step1[0] *  -cospi_2_64 + step1[7] * cospi_30_64;
        out[7] = (tran_low_t)fdct_round_shift(temp1);
        out[15] = (tran_low_t)fdct_round_shift(temp2);
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

#if CONFIG_VP9_HIGHBITDEPTH
void vp9_highbd_fdct4x4_c(const int16_t *input, tran_low_t *output,
                          int stride) {
  vp9_fdct4x4_c(input, output, stride);
}

void vp9_highbd_fdct8x8_c(const int16_t *input, tran_low_t *final_output,
                          int stride) {
  vp9_fdct8x8_c(input, final_output, stride);
}

void vp9_highbd_fdct16x16_c(const int16_t *input, tran_low_t *output,
                            int stride) {
  vp9_fdct16x16_c(input, output, stride);
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
