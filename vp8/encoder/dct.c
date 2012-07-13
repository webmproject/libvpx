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
#include "vpx_ports/config.h"
#include "vp8/common/idct.h"

#if CONFIG_INT_8X8FDCT

static const int xC1S7 = 16069;
static const int xC2S6 = 15137;
static const int xC3S5 = 13623;
static const int xC4S4 = 11585;
static const int xC5S3 =  9102;
static const int xC6S2 =  6270;
static const int xC7S1 =  3196;

#define SHIFT_BITS 14
#define DOROUND(X) X += (1<<(SHIFT_BITS-1));

#define FINAL_SHIFT 3
#define FINAL_ROUNDING (1<<(FINAL_SHIFT -1))
#define IN_SHIFT (FINAL_SHIFT+1)


void vp8_short_fdct8x8_c(short *InputData, short *OutputData, int pitch) {
  int loop;
  int short_pitch = pitch >> 1;
  int is07, is12, is34, is56;
  int is0734, is1256;
  int id07, id12, id34, id56;
  int irot_input_x, irot_input_y;
  int icommon_product1;      // Re-used product  (c4s4 * (s12 - s56))
  int icommon_product2;      // Re-used product  (c4s4 * (d12 + d56))
  int temp1, temp2;          // intermediate variable for computation

  int  InterData[64];
  int  *ip = InterData;
  short *op = OutputData;

  for (loop = 0; loop < 8; loop++) {
    // Pre calculate some common sums and differences.
    is07 = (InputData[0] + InputData[7]) << IN_SHIFT;
    is12 = (InputData[1] + InputData[2]) << IN_SHIFT;
    is34 = (InputData[3] + InputData[4]) << IN_SHIFT;
    is56 = (InputData[5] + InputData[6]) << IN_SHIFT;
    id07 = (InputData[0] - InputData[7]) << IN_SHIFT;
    id12 = (InputData[1] - InputData[2]) << IN_SHIFT;
    id34 = (InputData[3] - InputData[4]) << IN_SHIFT;
    id56 = (InputData[5] - InputData[6]) << IN_SHIFT;

    is0734 = is07 + is34;
    is1256 = is12 + is56;

    // Pre-Calculate some common product terms.
    icommon_product1 = xC4S4 * (is12 - is56);
    DOROUND(icommon_product1)
    icommon_product1 >>= SHIFT_BITS;

    icommon_product2 = xC4S4 * (id12 + id56);
    DOROUND(icommon_product2)
    icommon_product2 >>= SHIFT_BITS;


    ip[0] = (xC4S4 * (is0734 + is1256));
    DOROUND(ip[0]);
    ip[0] >>= SHIFT_BITS;

    ip[4] = (xC4S4 * (is0734 - is1256));
    DOROUND(ip[4]);
    ip[4] >>= SHIFT_BITS;

    // Define inputs to rotation for outputs 2 and 6
    irot_input_x = id12 - id56;
    irot_input_y = is07 - is34;

    // Apply rotation for outputs 2 and 6.
    temp1 = xC6S2 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC2S6 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    ip[2] = temp1 + temp2;

    temp1 = xC6S2 * irot_input_y;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC2S6 * irot_input_x;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    ip[6] = temp1 - temp2;

    // Define inputs to rotation for outputs 1 and 7
    irot_input_x = icommon_product1 + id07;
    irot_input_y = -(id34 + icommon_product2);

    // Apply rotation for outputs 1 and 7.
    temp1 = xC1S7 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC7S1 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    ip[1] = temp1 - temp2;

    temp1 = xC7S1 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC1S7 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    ip[7] = temp1 + temp2;

    // Define inputs to rotation for outputs 3 and 5
    irot_input_x = id07 - icommon_product1;
    irot_input_y = id34 - icommon_product2;

    // Apply rotation for outputs 3 and 5.
    temp1 = xC3S5 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC5S3 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    ip[3] = temp1 - temp2;


    temp1 = xC5S3 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC3S5 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    ip[5] = temp1 + temp2;

    // Increment data pointer for next row
    InputData += short_pitch;
    ip += 8;
  }

  // Performed DCT on rows, now transform the columns
  ip = InterData;
  for (loop = 0; loop < 8; loop++) {
    // Pre calculate some common sums and differences.
    is07 = ip[0 * 8] + ip[7 * 8];
    is12 = ip[1 * 8] + ip[2 * 8];
    is34 = ip[3 * 8] + ip[4 * 8];
    is56 = ip[5 * 8] + ip[6 * 8];

    id07 = ip[0 * 8] - ip[7 * 8];
    id12 = ip[1 * 8] - ip[2 * 8];
    id34 = ip[3 * 8] - ip[4 * 8];
    id56 = ip[5 * 8] - ip[6 * 8];

    is0734 = is07 + is34;
    is1256 = is12 + is56;

    // Pre-Calculate some common product terms
    icommon_product1 = xC4S4 * (is12 - is56);
    icommon_product2 = xC4S4 * (id12 + id56);
    DOROUND(icommon_product1)
    DOROUND(icommon_product2)
    icommon_product1 >>= SHIFT_BITS;
    icommon_product2 >>= SHIFT_BITS;


    temp1 = xC4S4 * (is0734 + is1256);
    temp2 = xC4S4 * (is0734 - is1256);
    DOROUND(temp1);
    DOROUND(temp2);
    temp1 >>= SHIFT_BITS;

    temp2 >>= SHIFT_BITS;
    op[0 * 8] = (temp1 + FINAL_ROUNDING) >> FINAL_SHIFT;
    op[4 * 8] = (temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;

    // Define inputs to rotation for outputs 2 and 6
    irot_input_x = id12 - id56;
    irot_input_y = is07 - is34;

    // Apply rotation for outputs 2 and 6.
    temp1 = xC6S2 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC2S6 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    op[2 * 8] = (temp1 + temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;

    temp1 = xC6S2 * irot_input_y;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC2S6 * irot_input_x;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    op[6 * 8] = (temp1 - temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;

    // Define inputs to rotation for outputs 1 and 7
    irot_input_x = icommon_product1 + id07;
    irot_input_y = -(id34 + icommon_product2);

    // Apply rotation for outputs 1 and 7.
    temp1 = xC1S7 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC7S1 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    op[1 * 8] = (temp1 - temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;

    temp1 = xC7S1 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC1S7 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    op[7 * 8] = (temp1 + temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;

    // Define inputs to rotation for outputs 3 and 5
    irot_input_x = id07 - icommon_product1;
    irot_input_y = id34 - icommon_product2;

    // Apply rotation for outputs 3 and 5.
    temp1 = xC3S5 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC5S3 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    op[3 * 8] = (temp1 - temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;


    temp1 = xC5S3 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC3S5 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    op[5 * 8] = (temp1 + temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;

    // Increment data pointer for next column.
    ip++;
    op++;
  }
}
#else

void vp8_short_fdct8x8_c(short *block, short *coefs, int pitch) {
  int j1, i, j, k;
  float b[8];
  float b1[8];
  float d[8][8];
  float f0 = (float) .7071068;
  float f1 = (float) .4903926;
  float f2 = (float) .4619398;
  float f3 = (float) .4157348;
  float f4 = (float) .3535534;
  float f5 = (float) .2777851;
  float f6 = (float) .1913417;
  float f7 = (float) .0975452;
  pitch = pitch / 2;
  for (i = 0, k = 0; i < 8; i++, k += pitch) {
    for (j = 0; j < 8; j++) {
      b[j] = (float)(block[k + j] << 3);
    }
    /* Horizontal transform */
    for (j = 0; j < 4; j++) {
      j1 = 7 - j;
      b1[j] = b[j] + b[j1];
      b1[j1] = b[j] - b[j1];
    }
    b[0] = b1[0] + b1[3];
    b[1] = b1[1] + b1[2];
    b[2] = b1[1] - b1[2];
    b[3] = b1[0] - b1[3];
    b[4] = b1[4];
    b[5] = (b1[6] - b1[5]) * f0;
    b[6] = (b1[6] + b1[5]) * f0;
    b[7] = b1[7];
    d[i][0] = (b[0] + b[1]) * f4;
    d[i][4] = (b[0] - b[1]) * f4;
    d[i][2] = b[2] * f6 + b[3] * f2;
    d[i][6] = b[3] * f6 - b[2] * f2;
    b1[4] = b[4] + b[5];
    b1[7] = b[7] + b[6];
    b1[5] = b[4] - b[5];
    b1[6] = b[7] - b[6];
    d[i][1] = b1[4] * f7 + b1[7] * f1;
    d[i][5] = b1[5] * f3 + b1[6] * f5;
    d[i][7] = b1[7] * f7 - b1[4] * f1;
    d[i][3] = b1[6] * f3 - b1[5] * f5;
  }
  /* Vertical transform */
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 4; j++) {
      j1 = 7 - j;
      b1[j] = d[j][i] + d[j1][i];
      b1[j1] = d[j][i] - d[j1][i];
    }
    b[0] = b1[0] + b1[3];
    b[1] = b1[1] + b1[2];
    b[2] = b1[1] - b1[2];
    b[3] = b1[0] - b1[3];
    b[4] = b1[4];
    b[5] = (b1[6] - b1[5]) * f0;
    b[6] = (b1[6] + b1[5]) * f0;
    b[7] = b1[7];
    d[0][i] = (b[0] + b[1]) * f4;
    d[4][i] = (b[0] - b[1]) * f4;
    d[2][i] = b[2] * f6 + b[3] * f2;
    d[6][i] = b[3] * f6 - b[2] * f2;
    b1[4] = b[4] + b[5];
    b1[7] = b[7] + b[6];
    b1[5] = b[4] - b[5];
    b1[6] = b[7] - b[6];
    d[1][i] = b1[4] * f7 + b1[7] * f1;
    d[5][i] = b1[5] * f3 + b1[6] * f5;
    d[7][i] = b1[7] * f7 - b1[4] * f1;
    d[3][i] = b1[6] * f3 - b1[5] * f5;
  }
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      *(coefs + j + i * 8) = (short) floor(d[i][j] + 0.5);
    }
  }
  return;
}

#endif

void vp8_short_fhaar2x2_c(short *input, short *output, int pitch) { // pitch = 8
  /* [1 1; 1 -1] orthogonal transform */
  /* use position: 0,1, 4, 8 */
  int i;
  short *ip1 = input;
  short *op1 = output;
  for (i = 0; i < 16; i++) {
    op1[i] = 0;
  }

  op1[0] = (ip1[0] + ip1[1] + ip1[4] + ip1[8] + 1) >> 1;
  op1[1] = (ip1[0] - ip1[1] + ip1[4] - ip1[8]) >> 1;
  op1[4] = (ip1[0] + ip1[1] - ip1[4] - ip1[8]) >> 1;
  op1[8] = (ip1[0] - ip1[1] - ip1[4] + ip1[8]) >> 1;

}

void vp8_short_fdct4x4_c(short *input, short *output, int pitch) {
  int i;
  int a1, b1, c1, d1;
  short *ip = input;
  short *op = output;

  for (i = 0; i < 4; i++) {
    a1 = ((ip[0] + ip[3]) << 5);
    b1 = ((ip[1] + ip[2]) << 5);
    c1 = ((ip[1] - ip[2]) << 5);
    d1 = ((ip[0] - ip[3]) << 5);

    op[0] = a1 + b1;
    op[2] = a1 - b1;

    op[1] = (c1 * 2217 + d1 * 5352 +  14500) >> 12;
    op[3] = (d1 * 2217 - c1 * 5352 +   7500) >> 12;

    ip += pitch / 2;
    op += 4;

  }
  ip = output;
  op = output;
  for (i = 0; i < 4; i++) {
    a1 = ip[0] + ip[12];
    b1 = ip[4] + ip[8];
    c1 = ip[4] - ip[8];
    d1 = ip[0] - ip[12];

    op[0]  = (a1 + b1 + 7) >> 4;
    op[8]  = (a1 - b1 + 7) >> 4;

    op[4]  = ((c1 * 2217 + d1 * 5352 +  12000) >> 16) + (d1 != 0);
    op[12] = (d1 * 2217 - c1 * 5352 +  51000) >> 16;

    ip++;
    op++;
  }
}

void vp8_short_fdct8x4_c(short *input, short *output, int pitch) {
  vp8_short_fdct4x4_c(input,   output,    pitch);
  vp8_short_fdct4x4_c(input + 4, output + 16, pitch);
}

void vp8_short_walsh4x4_c(short *input, short *output, int pitch) {
  int i;
  int a1, b1, c1, d1;
  short *ip = input;
  short *op = output;
  int pitch_short = pitch >> 1;

  for (i = 0; i < 4; i++) {
    a1 = ip[0 * pitch_short] + ip[3 * pitch_short];
    b1 = ip[1 * pitch_short] + ip[2 * pitch_short];
    c1 = ip[1 * pitch_short] - ip[2 * pitch_short];
    d1 = ip[0 * pitch_short] - ip[3 * pitch_short];

    op[0] = (a1 + b1 + 1) >> 1;
    op[4] = (c1 + d1) >> 1;
    op[8] = (a1 - b1) >> 1;
    op[12] = (d1 - c1) >> 1;

    ip++;
    op++;
  }
  ip = output;
  op = output;

  for (i = 0; i < 4; i++) {
    a1 = ip[0] + ip[3];
    b1 = ip[1] + ip[2];
    c1 = ip[1] - ip[2];
    d1 = ip[0] - ip[3];

    op[0] = (a1 + b1 + 1) >> 1;
    op[1] = (c1 + d1) >> 1;
    op[2] = (a1 - b1) >> 1;
    op[3] = (d1 - c1) >> 1;

    ip += 4;
    op += 4;
  }
}

#if CONFIG_LOSSLESS
void vp8_short_walsh4x4_lossless_c(short *input, short *output, int pitch) {
  int i;
  int a1, b1, c1, d1;
  short *ip = input;
  short *op = output;
  int pitch_short = pitch >> 1;

  for (i = 0; i < 4; i++) {
    a1 = (ip[0 * pitch_short] + ip[3 * pitch_short]) >> Y2_WHT_UPSCALE_FACTOR;
    b1 = (ip[1 * pitch_short] + ip[2 * pitch_short]) >> Y2_WHT_UPSCALE_FACTOR;
    c1 = (ip[1 * pitch_short] - ip[2 * pitch_short]) >> Y2_WHT_UPSCALE_FACTOR;
    d1 = (ip[0 * pitch_short] - ip[3 * pitch_short]) >> Y2_WHT_UPSCALE_FACTOR;

    op[0] = (a1 + b1 + 1) >> 1;
    op[4] = (c1 + d1) >> 1;
    op[8] = (a1 - b1) >> 1;
    op[12] = (d1 - c1) >> 1;

    ip++;
    op++;
  }
  ip = output;
  op = output;

  for (i = 0; i < 4; i++) {
    a1 = ip[0] + ip[3];
    b1 = ip[1] + ip[2];
    c1 = ip[1] - ip[2];
    d1 = ip[0] - ip[3];

    op[0] = ((a1 + b1 + 1) >> 1) << Y2_WHT_UPSCALE_FACTOR;
    op[1] = ((c1 + d1) >> 1) << Y2_WHT_UPSCALE_FACTOR;
    op[2] = ((a1 - b1) >> 1) << Y2_WHT_UPSCALE_FACTOR;
    op[3] = ((d1 - c1) >> 1) << Y2_WHT_UPSCALE_FACTOR;

    ip += 4;
    op += 4;
  }
}

void vp8_short_walsh4x4_x8_c(short *input, short *output, int pitch) {
  int i;
  int a1, b1, c1, d1;
  short *ip = input;
  short *op = output;
  int pitch_short = pitch >> 1;

  for (i = 0; i < 4; i++) {
    a1 = ip[0 * pitch_short] + ip[3 * pitch_short];
    b1 = ip[1 * pitch_short] + ip[2 * pitch_short];
    c1 = ip[1 * pitch_short] - ip[2 * pitch_short];
    d1 = ip[0 * pitch_short] - ip[3 * pitch_short];

    op[0] = (a1 + b1 + 1) >> 1;
    op[4] = (c1 + d1) >> 1;
    op[8] = (a1 - b1) >> 1;
    op[12] = (d1 - c1) >> 1;

    ip++;
    op++;
  }
  ip = output;
  op = output;

  for (i = 0; i < 4; i++) {
    a1 = ip[0] + ip[3];
    b1 = ip[1] + ip[2];
    c1 = ip[1] - ip[2];
    d1 = ip[0] - ip[3];

    op[0] = ((a1 + b1 + 1) >> 1) << WHT_UPSCALE_FACTOR;
    op[1] = ((c1 + d1) >> 1) << WHT_UPSCALE_FACTOR;
    op[2] = ((a1 - b1) >> 1) << WHT_UPSCALE_FACTOR;
    op[3] = ((d1 - c1) >> 1) << WHT_UPSCALE_FACTOR;

    ip += 4;
    op += 4;
  }
}

void vp8_short_walsh8x4_x8_c(short *input, short *output, int pitch) {
  vp8_short_walsh4x4_x8_c(input,   output,    pitch);
  vp8_short_walsh4x4_x8_c(input + 4, output + 16, pitch);
}
#endif
