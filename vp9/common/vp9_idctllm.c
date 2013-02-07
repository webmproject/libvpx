/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


/****************************************************************************
 * Notes:
 *
 * This implementation makes use of 16 bit fixed point verio of two multiply
 * constants:
 *         1.   sqrt(2) * cos (pi/8)
 *         2.   sqrt(2) * sin (pi/8)
 * Becuase the first constant is bigger than 1, to maintain the same 16 bit
 * fixed point precision as the second one, we use a trick of
 *         x * a = x + x*(a-1)
 * so
 *         x * sqrt(2) * cos (pi/8) = x + x * (sqrt(2) *cos(pi/8)-1).
 **************************************************************************/
#include <assert.h>
#include <math.h>
#include "./vpx_config.h"
#include "vp9/common/vp9_systemdependent.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_idct.h"



static const int16_t idct_i4[16] = {
  8192,  10703,  8192,   4433,
  8192,   4433, -8192, -10703,
  8192,  -4433, -8192,  10703,
  8192, -10703,  8192,  -4433
};

static const int16_t iadst_i4[16] = {
   3736,  9459, 10757,   7021,
   7021,  9459, -3736, -10757,
   9459,     0, -9459,   9459,
  10757, -9459,  7021,  -3736
};

static const int16_t idct_i8[64] = {
   5793,  8035,  7568,  6811,
   5793,  4551,  3135,  1598,
   5793,  6811,  3135, -1598,
  -5793, -8035, -7568, -4551,
   5793,  4551, -3135, -8035,
  -5793,  1598,  7568,  6811,
   5793,  1598, -7568, -4551,
   5793,  6811, -3135, -8035,
   5793, -1598, -7568,  4551,
   5793, -6811, -3135,  8035,
   5793, -4551, -3135,  8035,
  -5793, -1598,  7568, -6811,
   5793, -6811,  3135,  1598,
  -5793,  8035, -7568,  4551,
   5793, -8035,  7568, -6811,
   5793, -4551,  3135, -1598
};

static const int16_t iadst_i8[64] = {
   1460,  4184,  6342,  7644,
   7914,  7114,  5354,  2871,
   2871,  7114,  7644,  4184,
  -1460, -6342, -7914, -5354,
   4184,  7914,  2871, -5354,
  -7644, -1460,  6342,  7114,
   5354,  6342, -4184, -7114,
   2871,  7644, -1460, -7914,
   6342,  2871, -7914,  1460,
   7114, -5354, -4184,  7644,
   7114, -1460, -5354,  7914,
  -4184, -2871,  7644, -6342,
   7644, -5354,  1460,  2871,
  -6342,  7914, -7114,  4184,
   7914, -7644,  7114, -6342,
   5354, -4184,  2871, -1460
};



static const int16_t idct_i16[256] = {
   4096,  5765,  5681,  5543,  5352,  5109,  4816,  4478,
   4096,  3675,  3218,  2731,  2217,  1682,  1130,   568,
   4096,  5543,  4816,  3675,  2217,   568, -1130, -2731,
  -4096, -5109, -5681, -5765, -5352, -4478, -3218, -1682,
   4096,  5109,  3218,   568, -2217, -4478, -5681, -5543,
  -4096, -1682,  1130,  3675,  5352,  5765,  4816,  2731,
   4096,  4478,  1130, -2731, -5352, -5543, -3218,   568,
   4096,  5765,  4816,  1682, -2217, -5109, -5681, -3675,
   4096,  3675, -1130, -5109, -5352, -1682,  3218,  5765,
   4096,  -568, -4816, -5543, -2217,  2731,  5681,  4478,
   4096,  2731, -3218, -5765, -2217,  3675,  5681,  1682,
  -4096, -5543, -1130,  4478,  5352,   568, -4816, -5109,
   4096,  1682, -4816, -4478,  2217,  5765,  1130, -5109,
  -4096,  2731,  5681,   568, -5352, -3675,  3218,  5543,
   4096,   568, -5681, -1682,  5352,  2731, -4816, -3675,
   4096,  4478, -3218, -5109,  2217,  5543, -1130, -5765,
   4096,  -568, -5681,  1682,  5352, -2731, -4816,  3675,
   4096, -4478, -3218,  5109,  2217, -5543, -1130,  5765,
   4096, -1682, -4816,  4478,  2217, -5765,  1130,  5109,
  -4096, -2731,  5681,  -568, -5352,  3675,  3218, -5543,
   4096, -2731, -3218,  5765, -2217, -3675,  5681, -1682,
  -4096,  5543, -1130, -4478,  5352,  -568, -4816,  5109,
   4096, -3675, -1130,  5109, -5352,  1682,  3218, -5765,
   4096,   568, -4816,  5543, -2217, -2731,  5681, -4478,
   4096, -4478,  1130,  2731, -5352,  5543, -3218,  -568,
   4096, -5765,  4816, -1682, -2217,  5109, -5681,  3675,
   4096, -5109,  3218,  -568, -2217,  4478, -5681,  5543,
  -4096,  1682,  1130, -3675,  5352, -5765,  4816, -2731,
   4096, -5543,  4816, -3675,  2217,  -568, -1130,  2731,
  -4096,  5109, -5681,  5765, -5352,  4478, -3218,  1682,
   4096, -5765,  5681, -5543,  5352, -5109,  4816, -4478,
   4096, -3675,  3218, -2731,  2217, -1682,  1130,  -568
};

#if CONFIG_INTHT
static const int16_t iadst_i16[256] = {
   284,   850,  1407,  1951,  2476,  2977,  3450,  3889,
  4291,  4652,  4967,  5235,  5453,  5618,  5729,  5784,
   850,  2476,  3889,  4967,  5618,  5784,  5453,  4652,
  3450,  1951,   284, -1407, -2977, -4291, -5235, -5729,
  1407,  3889,  5453,  5729,  4652,  2476,  -284, -2977,
 -4967, -5784, -5235, -3450,  -850,  1951,  4291,  5618,
  1951,  4967,  5729,  3889,   284, -3450, -5618, -5235,
 -2476,  1407,  4652,  5784,  4291,   850, -2977, -5453,
  2476,  5618,  4652,   284, -4291, -5729, -2977,  1951,
  5453,  4967,   850, -3889, -5784, -3450,  1407,  5235,
  2977,  5784,  2476, -3450, -5729, -1951,  3889,  5618,
  1407, -4291, -5453,  -850,  4652,  5235,   284, -4967,
  3450,  5453,  -284, -5618, -2977,  3889,  5235,  -850,
 -5729, -2476,  4291,  4967, -1407, -5784, -1951,  4652,
  3889,  4652, -2977, -5235,  1951,  5618,  -850, -5784,
  -284,  5729,  1407, -5453, -2476,  4967,  3450, -4291,
  4291,  3450, -4967, -2476,  5453,  1407, -5729,  -284,
  5784,  -850, -5618,  1951,  5235, -2977, -4652,  3889,
  4652,  1951, -5784,  1407,  4967, -4291, -2476,  5729,
  -850, -5235,  3889,  2977, -5618,   284,  5453, -3450,
  4967,   284, -5235,  4652,   850, -5453,  4291,  1407,
 -5618,  3889,  1951, -5729,  3450,  2476, -5784,  2977,
  5235, -1407, -3450,  5784, -3889,  -850,  4967, -5453,
  1951,  2977, -5729,  4291,   284, -4652,  5618, -2476,
  5453, -2977,  -850,  4291, -5784,  4652, -1407, -2476,
  5235, -5618,  3450,   284, -3889,  5729, -4967,  1951,
  5618, -4291,  1951,   850, -3450,  5235, -5784,  4967,
 -2977,   284,  2476, -4652,  5729, -5453,  3889, -1407,
  5729, -5235,  4291, -2977,  1407,   284, -1951,  3450,
 -4652,  5453, -5784,  5618, -4967,  3889, -2476,   850,
  5784, -5729,  5618, -5453,  5235, -4967,  4652, -4291,
  3889, -3450,  2977, -2476,  1951, -1407,   850,  -284
};
#else
static const int16_t iadst_i16[256] = {
    542,  1607,  2614,  3526,  4311,  4940,  5390,  5646,
   5698,  5543,  5189,  4646,  3936,  3084,  2120,  1080,
   1080,  3084,  4646,  5543,  5646,  4940,  3526,  1607,
   -542, -2614, -4311, -5390, -5698, -5189, -3936, -2120,
   1607,  4311,  5646,  5189,  3084,     0, -3084, -5189,
  -5646, -4311, -1607,  1607,  4311,  5646,  5189,  3084,
   2120,  5189,  5390,  2614, -1607, -4940, -5543, -3084,
   1080,  4646,  5646,  3526, -542,  -4311, -5698, -3936,
   2614,  5646,  3936, -1080, -5189, -4940,  -542,  4311,
   5543,  2120, -3084, -5698, -3526,  1607,  5390,  4646,
   3084,  5646,  1607, -4311, -5189,     0,  5189,  4311,
  -1607, -5646, -3084,  3084,  5646,  1607, -4311, -5189,
   3526,  5189, -1080, -5698, -1607,  4940,  3936, -3084,
  -5390,   542,  5646,  2120, -4646, -4311,  2614,  5543,
   3936,  4311, -3526, -4646,  3084,  4940, -2614, -5189,
   2120,  5390, -1607, -5543,  1080,  5646,  -542, -5698,
   4311,  3084, -5189, -1607,  5646,     0, -5646,  1607,
   5189, -3084, -4311,  4311,  3084, -5189, -1607,  5646,
   4646,  1607, -5698,  2120,  4311, -4940, -1080,  5646,
  -2614, -3936,  5189,   542, -5543,  3084,  3526, -5390,
   4940,     0, -4940,  4940,     0, -4940,  4940,     0,
  -4940,  4940,     0, -4940,  4940,     0, -4940,  4940,
   5189, -1607, -3084,  5646, -4311,     0,  4311, -5646,
   3084,  1607, -5189,  5189, -1607, -3084,  5646, -4311,
   5390, -3084,  -542,  3936, -5646,  4940, -2120, -1607,
   4646, -5698,  4311, -1080, -2614,  5189, -5543,  3526,
   5543, -4311,  2120,   542, -3084,  4940, -5698,  5189,
  -3526,  1080,  1607, -3936,  5390, -5646,  4646, -2614,
   5646, -5189,  4311, -3084,  1607,     0, -1607,  3084,
  -4311,  5189, -5646,  5646, -5189,  4311, -3084,  1607,
   5698, -5646,  5543, -5390,  5189, -4940,  4646, -4311,
   3936, -3526,  3084, -2614,  2120, -1607,  1080,  -542
};
#endif

/* Converted the transforms to integer form. */
#define HORIZONTAL_SHIFT 14  // 16
#define HORIZONTAL_ROUNDING ((1 << (HORIZONTAL_SHIFT - 1)) - 1)
#define VERTICAL_SHIFT 17  // 15
#define VERTICAL_ROUNDING ((1 << (VERTICAL_SHIFT - 1)) - 1)
void vp9_ihtllm_c(const int16_t *input, int16_t *output, int pitch,
                      TX_TYPE tx_type, int tx_dim, uint16_t eobs) {
  int i, j, k;
  int nz_dim;
  int16_t imbuf[256];

  const int16_t *ip = input;
  int16_t *op = output;
  int16_t *im = &imbuf[0];

  /* pointers to vertical and horizontal transforms. */
  const int16_t *ptv = NULL, *pth = NULL;
  int shortpitch = pitch >> 1;

  switch (tx_type) {
    case ADST_ADST :
      ptv = pth = (tx_dim == 4) ? &iadst_i4[0]
                                  : ((tx_dim == 8) ? &iadst_i8[0]
                                                     : &iadst_i16[0]);
      break;
    case ADST_DCT  :
      ptv = (tx_dim == 4) ? &iadst_i4[0]
                            : ((tx_dim == 8) ? &iadst_i8[0] : &iadst_i16[0]);
      pth = (tx_dim == 4) ? &idct_i4[0]
                            : ((tx_dim == 8) ? &idct_i8[0] : &idct_i16[0]);
      break;
    case  DCT_ADST :
      ptv = (tx_dim == 4) ? &idct_i4[0]
                            : ((tx_dim == 8) ? &idct_i8[0] : &idct_i16[0]);
      pth = (tx_dim == 4) ? &iadst_i4[0]
                            : ((tx_dim == 8) ? &iadst_i8[0] : &iadst_i16[0]);
      break;
    case  DCT_DCT :
      ptv = pth = (tx_dim == 4) ? &idct_i4[0]
                                  : ((tx_dim == 8) ? &idct_i8[0]
                                                     : &idct_i16[0]);
      break;
    default:
      assert(0);
      break;
  }

  nz_dim = tx_dim;
  if(tx_dim > 4) {
    if(eobs < 36) {
      vpx_memset(im, 0, 512);
      nz_dim = 8;
      if(eobs < 3) {
        nz_dim = 2;
      } else if(eobs < 10) {
        nz_dim = 4;
      }
    }
  }

  /* 2-D inverse transform X = M1*Z*Transposed_M2 is calculated in 2 steps
   * from right to left:
   * 1. horizontal transform: Y= Z*Transposed_M2
   * 2. vertical transform: X = M1*Y
   * In SIMD, doing this way could eliminate the transpose needed if it is
   * calculated from left to right.
   */
  /* Horizontal transformation */
  for (j = 0; j < tx_dim; j++) {
    for (i = 0; i < nz_dim; i++) {
      int temp = 0;

      for (k = 0; k < nz_dim; k++) {
        temp += ip[k] * pth[k];
      }

      /* Calculate im and store it in its transposed position. */
      im[i] = (int16_t)((temp + HORIZONTAL_ROUNDING) >> HORIZONTAL_SHIFT);
      ip += tx_dim;
    }
    im += tx_dim;
    pth += tx_dim;
    ip = input;
  }

  /* Vertical transformation */
  im = &imbuf[0];

  for (i = 0; i < tx_dim; i++) {
    for (j = 0; j < tx_dim; j++) {
      int temp = 0;

      for (k = 0; k < nz_dim; k++) {
        temp += ptv[k] * im[k];
      }

      op[j] = (int16_t)((temp + VERTICAL_ROUNDING) >> VERTICAL_SHIFT);
      im += tx_dim;
    }
    im = &imbuf[0];
    ptv += tx_dim;
    op += shortpitch;
  }
}

void vp9_short_inv_walsh4x4_c(int16_t *input, int16_t *output) {
  int i;
  int a1, b1, c1, d1;
  int16_t *ip = input;
  int16_t *op = output;

  for (i = 0; i < 4; i++) {
    a1 = ((ip[0] + ip[3]));
    b1 = ((ip[1] + ip[2]));
    c1 = ((ip[1] - ip[2]));
    d1 = ((ip[0] - ip[3]));

    op[0] = (a1 + b1 + 1) >> 1;
    op[1] = (c1 + d1) >> 1;
    op[2] = (a1 - b1) >> 1;
    op[3] = (d1 - c1) >> 1;

    ip += 4;
    op += 4;
  }

  ip = output;
  op = output;
  for (i = 0; i < 4; i++) {
    a1 = ip[0] + ip[12];
    b1 = ip[4] + ip[8];
    c1 = ip[4] - ip[8];
    d1 = ip[0] - ip[12];
    op[0] = (a1 + b1 + 1) >> 1;
    op[4] = (c1 + d1) >> 1;
    op[8] = (a1 - b1) >> 1;
    op[12] = (d1 - c1) >> 1;
    ip++;
    op++;
  }
}

void vp9_short_inv_walsh4x4_1_c(int16_t *in, int16_t *out) {
  int i;
  int16_t tmp[4];
  int16_t *ip = in;
  int16_t *op = tmp;

  op[0] = (ip[0] + 1) >> 1;
  op[1] = op[2] = op[3] = (ip[0] >> 1);

  ip = tmp;
  op = out;
  for (i = 0; i < 4; i++) {
    op[0] = (ip[0] + 1) >> 1;
    op[4] = op[8] = op[12] = (ip[0] >> 1);
    ip++;
    op++;
  }
}

#if CONFIG_LOSSLESS
void vp9_short_inv_walsh4x4_lossless_c(int16_t *input, int16_t *output) {
  int i;
  int a1, b1, c1, d1;
  int16_t *ip = input;
  int16_t *op = output;

  for (i = 0; i < 4; i++) {
    a1 = ((ip[0] + ip[3])) >> Y2_WHT_UPSCALE_FACTOR;
    b1 = ((ip[1] + ip[2])) >> Y2_WHT_UPSCALE_FACTOR;
    c1 = ((ip[1] - ip[2])) >> Y2_WHT_UPSCALE_FACTOR;
    d1 = ((ip[0] - ip[3])) >> Y2_WHT_UPSCALE_FACTOR;

    op[0] = (a1 + b1 + 1) >> 1;
    op[1] = (c1 + d1) >> 1;
    op[2] = (a1 - b1) >> 1;
    op[3] = (d1 - c1) >> 1;

    ip += 4;
    op += 4;
  }

  ip = output;
  op = output;
  for (i = 0; i < 4; i++) {
    a1 = ip[0] + ip[12];
    b1 = ip[4] + ip[8];
    c1 = ip[4] - ip[8];
    d1 = ip[0] - ip[12];


    op[0] = ((a1 + b1 + 1) >> 1) << Y2_WHT_UPSCALE_FACTOR;
    op[4] = ((c1 + d1) >> 1) << Y2_WHT_UPSCALE_FACTOR;
    op[8] = ((a1 - b1) >> 1) << Y2_WHT_UPSCALE_FACTOR;
    op[12] = ((d1 - c1) >> 1) << Y2_WHT_UPSCALE_FACTOR;

    ip++;
    op++;
  }
}

void vp9_short_inv_walsh4x4_1_lossless_c(int16_t *in, int16_t *out) {
  int i;
  int16_t tmp[4];
  int16_t *ip = in;
  int16_t *op = tmp;

  op[0] = ((ip[0] >> Y2_WHT_UPSCALE_FACTOR) + 1) >> 1;
  op[1] = op[2] = op[3] = ((ip[0] >> Y2_WHT_UPSCALE_FACTOR) >> 1);

  ip = tmp;
  op = out;
  for (i = 0; i < 4; i++) {
    op[0] = ((ip[0] + 1) >> 1) << Y2_WHT_UPSCALE_FACTOR;
    op[4] = op[8] = op[12] = ((ip[0] >> 1)) << Y2_WHT_UPSCALE_FACTOR;
    ip++;
    op++;
  }
}

void vp9_short_inv_walsh4x4_x8_c(int16_t *input, int16_t *output, int pitch) {
  int i;
  int a1, b1, c1, d1;
  int16_t *ip = input;
  int16_t *op = output;
  int shortpitch = pitch >> 1;

  for (i = 0; i < 4; i++) {
    a1 = ((ip[0] + ip[3])) >> WHT_UPSCALE_FACTOR;
    b1 = ((ip[1] + ip[2])) >> WHT_UPSCALE_FACTOR;
    c1 = ((ip[1] - ip[2])) >> WHT_UPSCALE_FACTOR;
    d1 = ((ip[0] - ip[3])) >> WHT_UPSCALE_FACTOR;

    op[0] = (a1 + b1 + 1) >> 1;
    op[1] = (c1 + d1) >> 1;
    op[2] = (a1 - b1) >> 1;
    op[3] = (d1 - c1) >> 1;

    ip += 4;
    op += shortpitch;
  }

  ip = output;
  op = output;
  for (i = 0; i < 4; i++) {
    a1 = ip[shortpitch * 0] + ip[shortpitch * 3];
    b1 = ip[shortpitch * 1] + ip[shortpitch * 2];
    c1 = ip[shortpitch * 1] - ip[shortpitch * 2];
    d1 = ip[shortpitch * 0] - ip[shortpitch * 3];


    op[shortpitch * 0] = (a1 + b1 + 1) >> 1;
    op[shortpitch * 1] = (c1 + d1) >> 1;
    op[shortpitch * 2] = (a1 - b1) >> 1;
    op[shortpitch * 3] = (d1 - c1) >> 1;

    ip++;
    op++;
  }
}

void vp9_short_inv_walsh4x4_1_x8_c(int16_t *in, int16_t *out, int pitch) {
  int i;
  int16_t tmp[4];
  int16_t *ip = in;
  int16_t *op = tmp;
  int shortpitch = pitch >> 1;

  op[0] = ((ip[0] >> WHT_UPSCALE_FACTOR) + 1) >> 1;
  op[1] = op[2] = op[3] = ((ip[0] >> WHT_UPSCALE_FACTOR) >> 1);


  ip = tmp;
  op = out;
  for (i = 0; i < 4; i++) {
    op[shortpitch * 0] = (ip[0] + 1) >> 1;
    op[shortpitch * 1] = op[shortpitch * 2] = op[shortpitch * 3] = ip[0] >> 1;
    ip++;
    op++;
  }
}

void vp9_dc_only_inv_walsh_add_c(short input_dc, uint8_t *pred_ptr,
                                 uint8_t *dst_ptr,
                                 int pitch, int stride) {
  int r, c;
  short tmp[16];
  vp9_short_inv_walsh4x4_1_x8_c(&input_dc, tmp, 4 << 1);

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      dst_ptr[c] = clip_pixel(tmp[r * 4 + c] + pred_ptr[c]);
    }

    dst_ptr += stride;
    pred_ptr += pitch;
  }
}
#endif


void idct4_1d(int16_t *input, int16_t *output) {
  int16_t step[4];
  int temp1, temp2;
  // stage 1
  temp1 = (input[0] + input[2]) * cospi_16_64;
  temp2 = (input[0] - input[2]) * cospi_16_64;
  step[0] = dct_const_round_shift(temp1);
  step[1] = dct_const_round_shift(temp2);
  temp1 = input[1] * cospi_24_64 - input[3] * cospi_8_64;
  temp2 = input[1] * cospi_8_64 + input[3] * cospi_24_64;
  step[2] = dct_const_round_shift(temp1);
  step[3] = dct_const_round_shift(temp2);

  // stage 2
  output[0] = step[0] + step[3];
  output[1] = step[1] + step[2];
  output[2] = step[1] - step[2];
  output[3] = step[0] - step[3];
}

void vp9_short_idct4x4llm_c(int16_t *input, int16_t *output, int pitch) {
  int16_t out[4 * 4];
  int16_t *outptr = &out[0];
  const int short_pitch = pitch >> 1;
  int i, j;
  int16_t temp_in[4], temp_out[4];
  // First transform rows
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      temp_in[j] = input[j];
    idct4_1d(temp_in, outptr);
    input += 4;
    outptr += 4;
  }
  // Then transform columns
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      temp_in[j] = out[j * 4 + i];
    idct4_1d(temp_in, temp_out);
    for (j = 0; j < 4; ++j)
      output[j * short_pitch + i] = (temp_out[j] + 8) >> 4;
  }
}

void vp9_short_idct4x4llm_1_c(int16_t *input, int16_t *output, int pitch) {
  int i;
  int a1;
  int16_t *op = output;
  int shortpitch = pitch >> 1;
  int tmp;
  int16_t out;
  tmp = input[0] * cospi_16_64;
  out = dct_const_round_shift(tmp);
  tmp = out * cospi_16_64;
  out = dct_const_round_shift(tmp);
  a1 = (out + 8) >> 4;

  for (i = 0; i < 4; i++) {
    op[0] = a1;
    op[1] = a1;
    op[2] = a1;
    op[3] = a1;
    op += shortpitch;
  }
}

void vp9_dc_only_idct_add_c(int input_dc, uint8_t *pred_ptr,
                            uint8_t *dst_ptr, int pitch, int stride) {
  int a1;
  int r, c;
  int tmp;
  int16_t out;
  tmp = input_dc * cospi_16_64;
  out = dct_const_round_shift(tmp);
  tmp = out * cospi_16_64;
  out = dct_const_round_shift(tmp);
  a1 = (out + 8) >> 4;

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      dst_ptr[c] = clip_pixel(a1 + pred_ptr[c]);
    }
    dst_ptr += stride;
    pred_ptr += pitch;
  }
}

void idct8_1d(int16_t *input, int16_t *output) {
  int16_t step1[8], step2[8];
  int temp1, temp2;
  // stage 1
  step1[0] = input[0];
  step1[2] = input[4];
  step1[1] = input[2];
  step1[3] = input[6];
  temp1 = input[1] * cospi_28_64 - input[7] * cospi_4_64;
  temp2 = input[1] * cospi_4_64 + input[7] * cospi_28_64;
  step1[4] = dct_const_round_shift(temp1);
  step1[7] = dct_const_round_shift(temp2);
  temp1 = input[5] * cospi_12_64 - input[3] * cospi_20_64;
  temp2 = input[5] * cospi_20_64 + input[3] * cospi_12_64;
  step1[5] = dct_const_round_shift(temp1);
  step1[6] = dct_const_round_shift(temp2);

  // stage 2 & stage 3 - even half
  idct4_1d(step1, step1);

  // stage 2 - odd half
  step2[4] = step1[4] + step1[5];
  step2[5] = step1[4] - step1[5];
  step2[6] = -step1[6] + step1[7];
  step2[7] = step1[6] + step1[7];

  // stage 3 -odd half
  step1[4] = step2[4];
  temp1 = (step2[6] - step2[5]) * cospi_16_64;
  temp2 = (step2[5] + step2[6]) * cospi_16_64;
  step1[5] = dct_const_round_shift(temp1);
  step1[6] = dct_const_round_shift(temp2);
  step1[7] = step2[7];

  // stage 4
  output[0] = step1[0] + step1[7];
  output[1] = step1[1] + step1[6];
  output[2] = step1[2] + step1[5];
  output[3] = step1[3] + step1[4];
  output[4] = step1[3] - step1[4];
  output[5] = step1[2] - step1[5];
  output[6] = step1[1] - step1[6];
  output[7] = step1[0] - step1[7];
}

void vp9_short_idct8x8_c(int16_t *input, int16_t *output, int pitch) {
  int16_t out[8 * 8];
  int16_t *outptr = &out[0];
  const int short_pitch = pitch >> 1;
  int i, j;
  int16_t temp_in[8], temp_out[8];

  // First transform rows
  for (i = 0; i < 8; ++i) {
    idct8_1d(input, outptr);
    input += 8;
    outptr += 8;
  }

  // Then transform columns
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = out[j * 8 + i];
    idct8_1d(temp_in, temp_out);
    for (j = 0; j < 8; ++j)
        output[j * short_pitch + i] = (temp_out[j] + 16) >> 5;
    }
}

#if CONFIG_INTHT
static void iadst8_1d(int16_t *input, int16_t *output) {
  int x0, x1, x2, x3, x4, x5, x6, x7;
  int s0, s1, s2, s3, s4, s5, s6, s7;

  x0 = input[7];
  x1 = input[0];
  x2 = input[5];
  x3 = input[2];
  x4 = input[3];
  x5 = input[4];
  x6 = input[1];
  x7 = input[6];

  if (!(x0 | x1 | x2 | x3 | x4 | x5 | x6 | x7)) {
    output[0] = output[1] = output[2] = output[3] = output[4]
                    = output[5] = output[6] = output[7] = 0;
    return;
  }

  // stage 1
  s0 = cospi_2_64  * x0 + cospi_30_64 * x1;
  s1 = cospi_30_64 * x0 - cospi_2_64  * x1;
  s2 = cospi_10_64 * x2 + cospi_22_64 * x3;
  s3 = cospi_22_64 * x2 - cospi_10_64 * x3;
  s4 = cospi_18_64 * x4 + cospi_14_64 * x5;
  s5 = cospi_14_64 * x4 - cospi_18_64 * x5;
  s6 = cospi_26_64 * x6 + cospi_6_64  * x7;
  s7 = cospi_6_64  * x6 - cospi_26_64 * x7;

  x0 = dct_const_round_shift(s0 + s4);
  x1 = dct_const_round_shift(s1 + s5);
  x2 = dct_const_round_shift(s2 + s6);
  x3 = dct_const_round_shift(s3 + s7);
  x4 = dct_const_round_shift(s0 - s4);
  x5 = dct_const_round_shift(s1 - s5);
  x6 = dct_const_round_shift(s2 - s6);
  x7 = dct_const_round_shift(s3 - s7);

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
  x4 = dct_const_round_shift(s4 + s6);
  x5 = dct_const_round_shift(s5 + s7);
  x6 = dct_const_round_shift(s4 - s6);
  x7 = dct_const_round_shift(s5 - s7);

  // stage 3
  s2 = cospi_16_64 * (x2 + x3);
  s3 = cospi_16_64 * (x2 - x3);
  s6 = cospi_16_64 * (x6 + x7);
  s7 = cospi_16_64 * (x6 - x7);

  x2 = dct_const_round_shift(s2);
  x3 = dct_const_round_shift(s3);
  x6 = dct_const_round_shift(s6);
  x7 = dct_const_round_shift(s7);

  output[0] =   x0;
  output[1] = - x4;
  output[2] =   x6;
  output[3] = - x2;
  output[4] =   x3;
  output[5] = - x7;
  output[6] =   x5;
  output[7] = - x1;

  return;
}

void vp9_short_iht8x8_c(int16_t *input, int16_t *output,
                        TX_TYPE tx_type, int pitch) {
  int16_t out[8 * 8];
  int16_t *outptr = &out[0];
  const int short_pitch = pitch >> 1;
  int i, j;
  int16_t temp_in[8], temp_out[8];

  void (*invr)(int16_t*, int16_t*);
  void (*invc)(int16_t*, int16_t*);

  switch (tx_type) {
    case ADST_ADST:
      invc = &iadst8_1d;
      invr = &iadst8_1d;
      break;
    case ADST_DCT:
      invc = &iadst8_1d;
      invr = &idct8_1d;
      break;
    case DCT_ADST:
      invc = &idct8_1d;
      invr = &iadst8_1d;
      break;
    case DCT_DCT:
      invc = &idct8_1d;
      invr = &idct8_1d;
      break;
    default:
      assert(0);
  }

  // inverse transform row vectors
  for (i = 0; i < 8; ++i) {
    invr(input, outptr);
    input += 8;
    outptr += 8;
  }

  // inverse transform column vectors
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = out[j * 8 + i];
    invc(temp_in, temp_out);
    for (j = 0; j < 8; ++j)
      output[j * short_pitch + i] = (temp_out[j] + 16) >> 5;
  }
}
#endif


void vp9_short_idct10_8x8_c(int16_t *input, int16_t *output, int pitch) {
  int16_t out[8 * 8];
  int16_t *outptr = &out[0];
  const int short_pitch = pitch >> 1;
  int i, j;
  int16_t temp_in[8], temp_out[8];

  vpx_memset(out, 0, sizeof(out));
  // First transform rows
  // only first 4 row has non-zero coefs
  for (i = 0; i < 4; ++i) {
    idct8_1d(input, outptr);
    input += 8;
    outptr += 8;
  }

  // Then transform columns
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = out[j * 8 + i];
    idct8_1d(temp_in, temp_out);
    for (j = 0; j < 8; ++j)
        output[j * short_pitch + i] = (temp_out[j] + 16) >> 5;
    }
}

void vp9_short_idct1_8x8_c(int16_t *input, int16_t *output) {
  int tmp;
  int16_t out;
  tmp = input[0] * cospi_16_64;
  out = dct_const_round_shift(tmp);
  tmp = out * cospi_16_64;
  out = dct_const_round_shift(tmp);
  *output = (out + 16) >> 5;
}

void vp9_short_ihaar2x2_c(int16_t *input, int16_t *output, int pitch) {
  int i;
  int16_t *ip = input;  // 0, 1, 4, 8
  int16_t *op = output;
  for (i = 0; i < 16; i++) {
    op[i] = 0;
  }

  op[0] = (ip[0] + ip[1] + ip[4] + ip[8] + 1) >> 1;
  op[1] = (ip[0] - ip[1] + ip[4] - ip[8]) >> 1;
  op[4] = (ip[0] + ip[1] - ip[4] - ip[8]) >> 1;
  op[8] = (ip[0] - ip[1] - ip[4] + ip[8]) >> 1;
}

void idct16_1d(int16_t *input, int16_t *output) {
  int16_t step1[16], step2[16];
  int temp1, temp2;

  // stage 1
  step1[0] = input[0/2];
  step1[1] = input[16/2];
  step1[2] = input[8/2];
  step1[3] = input[24/2];
  step1[4] = input[4/2];
  step1[5] = input[20/2];
  step1[6] = input[12/2];
  step1[7] = input[28/2];
  step1[8] = input[2/2];
  step1[9] = input[18/2];
  step1[10] = input[10/2];
  step1[11] = input[26/2];
  step1[12] = input[6/2];
  step1[13] = input[22/2];
  step1[14] = input[14/2];
  step1[15] = input[30/2];

  // stage 2
  step2[0] = step1[0];
  step2[1] = step1[1];
  step2[2] = step1[2];
  step2[3] = step1[3];
  step2[4] = step1[4];
  step2[5] = step1[5];
  step2[6] = step1[6];
  step2[7] = step1[7];

  temp1 = step1[8] * cospi_30_64 - step1[15] * cospi_2_64;
  temp2 = step1[8] * cospi_2_64 + step1[15] * cospi_30_64;
  step2[8] = dct_const_round_shift(temp1);
  step2[15] = dct_const_round_shift(temp2);

  temp1 = step1[9] * cospi_14_64 - step1[14] * cospi_18_64;
  temp2 = step1[9] * cospi_18_64 + step1[14] * cospi_14_64;
  step2[9] = dct_const_round_shift(temp1);
  step2[14] = dct_const_round_shift(temp2);

  temp1 = step1[10] * cospi_22_64 - step1[13] * cospi_10_64;
  temp2 = step1[10] * cospi_10_64 + step1[13] * cospi_22_64;
  step2[10] = dct_const_round_shift(temp1);
  step2[13] = dct_const_round_shift(temp2);

  temp1 = step1[11] * cospi_6_64 - step1[12] * cospi_26_64;
  temp2 = step1[11] * cospi_26_64 + step1[12] * cospi_6_64;
  step2[11] = dct_const_round_shift(temp1);
  step2[12] = dct_const_round_shift(temp2);

  // stage 3
  step1[0] = step2[0];
  step1[1] = step2[1];
  step1[2] = step2[2];
  step1[3] = step2[3];

  temp1 = step2[4] * cospi_28_64 - step2[7] * cospi_4_64;
  temp2 = step2[4] * cospi_4_64 + step2[7] * cospi_28_64;
  step1[4] = dct_const_round_shift(temp1);
  step1[7] = dct_const_round_shift(temp2);
  temp1 = step2[5] * cospi_12_64 - step2[6] * cospi_20_64;
  temp2 = step2[5] * cospi_20_64 + step2[6] * cospi_12_64;
  step1[5] = dct_const_round_shift(temp1);
  step1[6] = dct_const_round_shift(temp2);

  step1[8] = step2[8] + step2[9];
  step1[9] = step2[8] - step2[9];
  step1[10] = -step2[10] + step2[11];
  step1[11] = step2[10] + step2[11];
  step1[12] = step2[12] + step2[13];
  step1[13] = step2[12] - step2[13];
  step1[14] = -step2[14] + step2[15];
  step1[15] = step2[14] + step2[15];

  temp1 = (step1[0] + step1[1]) * cospi_16_64;
  temp2 = (step1[0] - step1[1]) * cospi_16_64;
  step2[0] = dct_const_round_shift(temp1);
  step2[1] = dct_const_round_shift(temp2);
  temp1 = step1[2] * cospi_24_64 - step1[3] * cospi_8_64;
  temp2 = step1[2] * cospi_8_64 + step1[3] * cospi_24_64;
  step2[2] = dct_const_round_shift(temp1);
  step2[3] = dct_const_round_shift(temp2);
  step2[4] = step1[4] + step1[5];
  step2[5] = step1[4] - step1[5];
  step2[6] = -step1[6] + step1[7];
  step2[7] = step1[6] + step1[7];

  step2[8] = step1[8];
  step2[15] = step1[15];
  temp1 = -step1[9] * cospi_8_64 + step1[14] * cospi_24_64;
  temp2 = step1[9] * cospi_24_64 + step1[14] * cospi_8_64;
  step2[9] = dct_const_round_shift(temp1);
  step2[14] = dct_const_round_shift(temp2);
  temp1 = -step1[10] * cospi_24_64 - step1[13] * cospi_8_64;
  temp2 = -step1[10] * cospi_8_64 + step1[13] * cospi_24_64;
  step2[10] = dct_const_round_shift(temp1);
  step2[13] = dct_const_round_shift(temp2);
  step2[11] = step1[11];
  step2[12] = step1[12];

  // stage 5
  step1[0] = step2[0] + step2[3];
  step1[1] = step2[1] + step2[2];
  step1[2] = step2[1] - step2[2];
  step1[3] = step2[0] - step2[3];
  step1[4] = step2[4];
  temp1 = (step2[6] - step2[5]) * cospi_16_64;
  temp2 = (step2[5] + step2[6]) * cospi_16_64;
  step1[5] = dct_const_round_shift(temp1);
  step1[6] = dct_const_round_shift(temp2);
  step1[7] = step2[7];

  step1[8] = step2[8] + step2[11];
  step1[9] = step2[9] + step2[10];
  step1[10] = step2[9] - step2[10];
  step1[11] = step2[8] - step2[11];
  step1[12] = -step2[12] + step2[15];
  step1[13] = -step2[13] + step2[14];
  step1[14] = step2[13] + step2[14];
  step1[15] = step2[12] + step2[15];

  // stage 6
  step2[0] = step1[0] + step1[7];
  step2[1] = step1[1] + step1[6];
  step2[2] = step1[2] + step1[5];
  step2[3] = step1[3] + step1[4];
  step2[4] = step1[3] - step1[4];
  step2[5] = step1[2] - step1[5];
  step2[6] = step1[1] - step1[6];
  step2[7] = step1[0] - step1[7];
  step2[8] = step1[8];
  step2[9] = step1[9];
  temp1 = (-step1[10] + step1[13]) * cospi_16_64;
  temp2 = (step1[10] + step1[13]) * cospi_16_64;
  step2[10] = dct_const_round_shift(temp1);
  step2[13] = dct_const_round_shift(temp2);
  temp1 = (-step1[11] + step1[12]) * cospi_16_64;
  temp2 = (step1[11] + step1[12]) * cospi_16_64;
  step2[11] = dct_const_round_shift(temp1);
  step2[12] = dct_const_round_shift(temp2);
  step2[14] = step1[14];
  step2[15] = step1[15];

  // stage 7
  output[0] = step2[0] + step2[15];
  output[1] = step2[1] + step2[14];
  output[2] = step2[2] + step2[13];
  output[3] = step2[3] + step2[12];
  output[4] = step2[4] + step2[11];
  output[5] = step2[5] + step2[10];
  output[6] = step2[6] + step2[9];
  output[7] = step2[7] + step2[8];
  output[8] = step2[7] - step2[8];
  output[9] = step2[6] - step2[9];
  output[10] = step2[5] - step2[10];
  output[11] = step2[4] - step2[11];
  output[12] = step2[3] - step2[12];
  output[13] = step2[2] - step2[13];
  output[14] = step2[1] - step2[14];
  output[15] = step2[0] - step2[15];
}

void vp9_short_idct16x16_c(int16_t *input, int16_t *output, int pitch) {
  int16_t out[16 * 16];
  int16_t *outptr = &out[0];
  const int short_pitch = pitch >> 1;
  int i, j;
  int16_t temp_in[16], temp_out[16];

  // First transform rows
  for (i = 0; i < 16; ++i) {
    idct16_1d(input, outptr);
    input += short_pitch;
    outptr += 16;
  }

  // Then transform columns
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j)
      temp_in[j] = out[j * 16 + i];
    idct16_1d(temp_in, temp_out);
    for (j = 0; j < 16; ++j)
        output[j * 16 + i] = (temp_out[j] + 32) >> 6;
    }
}

void vp9_short_idct10_16x16_c(int16_t *input, int16_t *output, int pitch) {
    int16_t out[16 * 16];
    int16_t *outptr = &out[0];
    const int short_pitch = pitch >> 1;
    int i, j;
    int16_t temp_in[16], temp_out[16];

    /* First transform rows. Since all non-zero dct coefficients are in
     * upper-left 4x4 area, we only need to calculate first 4 rows here.
     */
    vpx_memset(out, 0, sizeof(out));
    for (i = 0; i < 4; ++i) {
      idct16_1d(input, outptr);
      input += short_pitch;
      outptr += 16;
    }

    // Then transform columns
    for (i = 0; i < 16; ++i) {
      for (j = 0; j < 16; ++j)
        temp_in[j] = out[j*16 + i];
      idct16_1d(temp_in, temp_out);
      for (j = 0; j < 16; ++j)
        output[j*16 + i] = (temp_out[j] + 32) >> 6;
    }
}


void vp9_short_idct1_16x16_c(int16_t *input, int16_t *output) {
  int tmp;
  int16_t out;
  tmp = input[0] * cospi_16_64;
  out = dct_const_round_shift(tmp);
  tmp = out * cospi_16_64;
  out = dct_const_round_shift(tmp);
  *output = (out + 32) >> 6;
}


#if !CONFIG_DWTDCTHYBRID
void idct32_1d(int16_t *input, int16_t *output) {
  int16_t step1[32], step2[32];
  int temp1, temp2;

  // stage 1
  step1[0] = input[0];
  step1[1] = input[16];
  step1[2] = input[8];
  step1[3] = input[24];
  step1[4] = input[4];
  step1[5] = input[20];
  step1[6] = input[12];
  step1[7] = input[28];
  step1[8] = input[2];
  step1[9] = input[18];
  step1[10] = input[10];
  step1[11] = input[26];
  step1[12] = input[6];
  step1[13] = input[22];
  step1[14] = input[14];
  step1[15] = input[30];

  temp1 = input[1] * cospi_31_64 - input[31] * cospi_1_64;
  temp2 = input[1] * cospi_1_64 + input[31] * cospi_31_64;
  step1[16] = dct_const_round_shift(temp1);
  step1[31] = dct_const_round_shift(temp2);

  temp1 = input[17] * cospi_15_64 - input[15] * cospi_17_64;
  temp2 = input[17] * cospi_17_64 + input[15] * cospi_15_64;
  step1[17] = dct_const_round_shift(temp1);
  step1[30] = dct_const_round_shift(temp2);

  temp1 = input[9] * cospi_23_64 - input[23] * cospi_9_64;
  temp2 = input[9] * cospi_9_64 + input[23] * cospi_23_64;
  step1[18] = dct_const_round_shift(temp1);
  step1[29] = dct_const_round_shift(temp2);

  temp1 = input[25] * cospi_7_64 - input[7] * cospi_25_64;
  temp2 = input[25] * cospi_25_64 + input[7] * cospi_7_64;
  step1[19] = dct_const_round_shift(temp1);
  step1[28] = dct_const_round_shift(temp2);

  temp1 = input[5] * cospi_27_64 - input[27] * cospi_5_64;
  temp2 = input[5] * cospi_5_64 + input[27] * cospi_27_64;
  step1[20] = dct_const_round_shift(temp1);
  step1[27] = dct_const_round_shift(temp2);

  temp1 = input[21] * cospi_11_64 - input[11] * cospi_21_64;
  temp2 = input[21] * cospi_21_64 + input[11] * cospi_11_64;
  step1[21] = dct_const_round_shift(temp1);
  step1[26] = dct_const_round_shift(temp2);

  temp1 = input[13] * cospi_19_64 - input[19] * cospi_13_64;
  temp2 = input[13] * cospi_13_64 + input[19] * cospi_19_64;
  step1[22] = dct_const_round_shift(temp1);
  step1[25] = dct_const_round_shift(temp2);

  temp1 = input[29] * cospi_3_64 - input[3] * cospi_29_64;
  temp2 = input[29] * cospi_29_64 + input[3] * cospi_3_64;
  step1[23] = dct_const_round_shift(temp1);
  step1[24] = dct_const_round_shift(temp2);

  // stage 2
  step2[0] = step1[0];
  step2[1] = step1[1];
  step2[2] = step1[2];
  step2[3] = step1[3];
  step2[4] = step1[4];
  step2[5] = step1[5];
  step2[6] = step1[6];
  step2[7] = step1[7];

  temp1 = step1[8] * cospi_30_64 - step1[15] * cospi_2_64;
  temp2 = step1[8] * cospi_2_64 + step1[15] * cospi_30_64;
  step2[8] = dct_const_round_shift(temp1);
  step2[15] = dct_const_round_shift(temp2);

  temp1 = step1[9] * cospi_14_64 - step1[14] * cospi_18_64;
  temp2 = step1[9] * cospi_18_64 + step1[14] * cospi_14_64;
  step2[9] = dct_const_round_shift(temp1);
  step2[14] = dct_const_round_shift(temp2);

  temp1 = step1[10] * cospi_22_64 - step1[13] * cospi_10_64;
  temp2 = step1[10] * cospi_10_64 + step1[13] * cospi_22_64;
  step2[10] = dct_const_round_shift(temp1);
  step2[13] = dct_const_round_shift(temp2);

  temp1 = step1[11] * cospi_6_64 - step1[12] * cospi_26_64;
  temp2 = step1[11] * cospi_26_64 + step1[12] * cospi_6_64;
  step2[11] = dct_const_round_shift(temp1);
  step2[12] = dct_const_round_shift(temp2);

  step2[16] = step1[16] + step1[17];
  step2[17] = step1[16] - step1[17];
  step2[18] = -step1[18] + step1[19];
  step2[19] = step1[18] + step1[19];
  step2[20] = step1[20] + step1[21];
  step2[21] = step1[20] - step1[21];
  step2[22] = -step1[22] + step1[23];
  step2[23] = step1[22] + step1[23];
  step2[24] = step1[24] + step1[25];
  step2[25] = step1[24] - step1[25];
  step2[26] = -step1[26] + step1[27];
  step2[27] = step1[26] + step1[27];
  step2[28] = step1[28] + step1[29];
  step2[29] = step1[28] - step1[29];
  step2[30] = -step1[30] + step1[31];
  step2[31] = step1[30] + step1[31];

  // stage 3
  step1[0] = step2[0];
  step1[1] = step2[1];
  step1[2] = step2[2];
  step1[3] = step2[3];

  temp1 = step2[4] * cospi_28_64 - step2[7] * cospi_4_64;
  temp2 = step2[4] * cospi_4_64 + step2[7] * cospi_28_64;
  step1[4] = dct_const_round_shift(temp1);
  step1[7] = dct_const_round_shift(temp2);
  temp1 = step2[5] * cospi_12_64 - step2[6] * cospi_20_64;
  temp2 = step2[5] * cospi_20_64 + step2[6] * cospi_12_64;
  step1[5] = dct_const_round_shift(temp1);
  step1[6] = dct_const_round_shift(temp2);

  step1[8] = step2[8] + step2[9];
  step1[9] = step2[8] - step2[9];
  step1[10] = -step2[10] + step2[11];
  step1[11] = step2[10] + step2[11];
  step1[12] = step2[12] + step2[13];
  step1[13] = step2[12] - step2[13];
  step1[14] = -step2[14] + step2[15];
  step1[15] = step2[14] + step2[15];

  step1[16] = step2[16];
  step1[31] = step2[31];
  temp1 = -step2[17] * cospi_4_64 + step2[30] * cospi_28_64;
  temp2 = step2[17] * cospi_28_64 + step2[30] * cospi_4_64;
  step1[17] = dct_const_round_shift(temp1);
  step1[30] = dct_const_round_shift(temp2);
  temp1 = -step2[18] * cospi_28_64 - step2[29] * cospi_4_64;
  temp2 = -step2[18] * cospi_4_64 + step2[29] * cospi_28_64;
  step1[18] = dct_const_round_shift(temp1);
  step1[29] = dct_const_round_shift(temp2);
  step1[19] = step2[19];
  step1[20] = step2[20];
  temp1 = -step2[21] * cospi_20_64 + step2[26] * cospi_12_64;
  temp2 = step2[21] * cospi_12_64 + step2[26] * cospi_20_64;
  step1[21] = dct_const_round_shift(temp1);
  step1[26] = dct_const_round_shift(temp2);
  temp1 = -step2[22] * cospi_12_64 - step2[25] * cospi_20_64;
  temp2 = -step2[22] * cospi_20_64 + step2[25] * cospi_12_64;
  step1[22] = dct_const_round_shift(temp1);
  step1[25] = dct_const_round_shift(temp2);
  step1[23] = step2[23];
  step1[24] = step2[24];
  step1[27] = step2[27];
  step1[28] = step2[28];

  // stage 4
  temp1 = (step1[0] + step1[1]) * cospi_16_64;
  temp2 = (step1[0] - step1[1]) * cospi_16_64;
  step2[0] = dct_const_round_shift(temp1);
  step2[1] = dct_const_round_shift(temp2);
  temp1 = step1[2] * cospi_24_64 - step1[3] * cospi_8_64;
  temp2 = step1[2] * cospi_8_64 + step1[3] * cospi_24_64;
  step2[2] = dct_const_round_shift(temp1);
  step2[3] = dct_const_round_shift(temp2);
  step2[4] = step1[4] + step1[5];
  step2[5] = step1[4] - step1[5];
  step2[6] = -step1[6] + step1[7];
  step2[7] = step1[6] + step1[7];

  step2[8] = step1[8];
  step2[15] = step1[15];
  temp1 = -step1[9] * cospi_8_64 + step1[14] * cospi_24_64;
  temp2 = step1[9] * cospi_24_64 + step1[14] * cospi_8_64;
  step2[9] = dct_const_round_shift(temp1);
  step2[14] = dct_const_round_shift(temp2);
  temp1 = -step1[10] * cospi_24_64 - step1[13] * cospi_8_64;
  temp2 = -step1[10] * cospi_8_64 + step1[13] * cospi_24_64;
  step2[10] = dct_const_round_shift(temp1);
  step2[13] = dct_const_round_shift(temp2);
  step2[11] = step1[11];
  step2[12] = step1[12];

  step2[16] = step1[16] + step1[19];
  step2[17] = step1[17] + step1[18];
  step2[18] = step1[17] - step1[18];
  step2[19] = step1[16] - step1[19];
  step2[20] = -step1[20] + step1[23];
  step2[21] = -step1[21] + step1[22];
  step2[22] = step1[21] + step1[22];
  step2[23] = step1[20] + step1[23];

  step2[24] = step1[24] + step1[27];
  step2[25] = step1[25] + step1[26];
  step2[26] = step1[25] - step1[26];
  step2[27] = step1[24] - step1[27];
  step2[28] = -step1[28] + step1[31];
  step2[29] = -step1[29] + step1[30];
  step2[30] = step1[29] + step1[30];
  step2[31] = step1[28] + step1[31];

  // stage 5
  step1[0] = step2[0] + step2[3];
  step1[1] = step2[1] + step2[2];
  step1[2] = step2[1] - step2[2];
  step1[3] = step2[0] - step2[3];
  step1[4] = step2[4];
  temp1 = (step2[6] - step2[5]) * cospi_16_64;
  temp2 = (step2[5] + step2[6]) * cospi_16_64;
  step1[5] = dct_const_round_shift(temp1);
  step1[6] = dct_const_round_shift(temp2);
  step1[7] = step2[7];

  step1[8] = step2[8] + step2[11];
  step1[9] = step2[9] + step2[10];
  step1[10] = step2[9] - step2[10];
  step1[11] = step2[8] - step2[11];
  step1[12] = -step2[12] + step2[15];
  step1[13] = -step2[13] + step2[14];
  step1[14] = step2[13] + step2[14];
  step1[15] = step2[12] + step2[15];

  step1[16] = step2[16];
  step1[17] = step2[17];
  temp1 = -step2[18] * cospi_8_64 + step2[29] * cospi_24_64;
  temp2 = step2[18] * cospi_24_64 + step2[29] * cospi_8_64;
  step1[18] = dct_const_round_shift(temp1);
  step1[29] = dct_const_round_shift(temp2);
  temp1 = -step2[19] * cospi_8_64 + step2[28] * cospi_24_64;
  temp2 = step2[19] * cospi_24_64 + step2[28] * cospi_8_64;
  step1[19] = dct_const_round_shift(temp1);
  step1[28] = dct_const_round_shift(temp2);
  temp1 = -step2[20] * cospi_24_64 - step2[27] * cospi_8_64;
  temp2 = -step2[20] * cospi_8_64 + step2[27] * cospi_24_64;
  step1[20] = dct_const_round_shift(temp1);
  step1[27] = dct_const_round_shift(temp2);
  temp1 = -step2[21] * cospi_24_64 - step2[26] * cospi_8_64;
  temp2 = -step2[21] * cospi_8_64 + step2[26] * cospi_24_64;
  step1[21] = dct_const_round_shift(temp1);
  step1[26] = dct_const_round_shift(temp2);
  step1[22] = step2[22];
  step1[23] = step2[23];
  step1[24] = step2[24];
  step1[25] = step2[25];
  step1[30] = step2[30];
  step1[31] = step2[31];

  // stage 6
  step2[0] = step1[0] + step1[7];
  step2[1] = step1[1] + step1[6];
  step2[2] = step1[2] + step1[5];
  step2[3] = step1[3] + step1[4];
  step2[4] = step1[3] - step1[4];
  step2[5] = step1[2] - step1[5];
  step2[6] = step1[1] - step1[6];
  step2[7] = step1[0] - step1[7];
  step2[8] = step1[8];
  step2[9] = step1[9];
  temp1 = (-step1[10] + step1[13]) * cospi_16_64;
  temp2 = (step1[10] + step1[13]) * cospi_16_64;
  step2[10] = dct_const_round_shift(temp1);
  step2[13] = dct_const_round_shift(temp2);
  temp1 = (-step1[11] + step1[12]) * cospi_16_64;
  temp2 = (step1[11] + step1[12]) * cospi_16_64;
  step2[11] = dct_const_round_shift(temp1);
  step2[12] = dct_const_round_shift(temp2);
  step2[14] = step1[14];
  step2[15] = step1[15];

  step2[16] = step1[16] + step1[23];
  step2[17] = step1[17] + step1[22];
  step2[18] = step1[18] + step1[21];
  step2[19] = step1[19] + step1[20];
  step2[20] = step1[19] - step1[20];
  step2[21] = step1[18] - step1[21];
  step2[22] = step1[17] - step1[22];
  step2[23] = step1[16] - step1[23];

  step2[24] = -step1[24] + step1[31];
  step2[25] = -step1[25] + step1[30];
  step2[26] = -step1[26] + step1[29];
  step2[27] = -step1[27] + step1[28];
  step2[28] = step1[27] + step1[28];
  step2[29] = step1[26] + step1[29];
  step2[30] = step1[25] + step1[30];
  step2[31] = step1[24] + step1[31];

  // stage 7
  step1[0] = step2[0] + step2[15];
  step1[1] = step2[1] + step2[14];
  step1[2] = step2[2] + step2[13];
  step1[3] = step2[3] + step2[12];
  step1[4] = step2[4] + step2[11];
  step1[5] = step2[5] + step2[10];
  step1[6] = step2[6] + step2[9];
  step1[7] = step2[7] + step2[8];
  step1[8] = step2[7] - step2[8];
  step1[9] = step2[6] - step2[9];
  step1[10] = step2[5] - step2[10];
  step1[11] = step2[4] - step2[11];
  step1[12] = step2[3] - step2[12];
  step1[13] = step2[2] - step2[13];
  step1[14] = step2[1] - step2[14];
  step1[15] = step2[0] - step2[15];

  step1[16] = step2[16];
  step1[17] = step2[17];
  step1[18] = step2[18];
  step1[19] = step2[19];
  temp1 = (-step2[20] + step2[27]) * cospi_16_64;
  temp2 = (step2[20] + step2[27]) * cospi_16_64;
  step1[20] = dct_const_round_shift(temp1);
  step1[27] = dct_const_round_shift(temp2);
  temp1 = (-step2[21] + step2[26]) * cospi_16_64;
  temp2 = (step2[21] + step2[26]) * cospi_16_64;
  step1[21] = dct_const_round_shift(temp1);
  step1[26] = dct_const_round_shift(temp2);
  temp1 = (-step2[22] + step2[25]) * cospi_16_64;
  temp2 = (step2[22] + step2[25]) * cospi_16_64;
  step1[22] = dct_const_round_shift(temp1);
  step1[25] = dct_const_round_shift(temp2);
  temp1 = (-step2[23] + step2[24]) * cospi_16_64;
  temp2 = (step2[23] + step2[24]) * cospi_16_64;
  step1[23] = dct_const_round_shift(temp1);
  step1[24] = dct_const_round_shift(temp2);
  step1[28] = step2[28];
  step1[29] = step2[29];
  step1[30] = step2[30];
  step1[31] = step2[31];

  // final stage
  output[0] = step1[0] + step1[31];
  output[1] = step1[1] + step1[30];
  output[2] = step1[2] + step1[29];
  output[3] = step1[3] + step1[28];
  output[4] = step1[4] + step1[27];
  output[5] = step1[5] + step1[26];
  output[6] = step1[6] + step1[25];
  output[7] = step1[7] + step1[24];
  output[8] = step1[8] + step1[23];
  output[9] = step1[9] + step1[22];
  output[10] = step1[10] + step1[21];
  output[11] = step1[11] + step1[20];
  output[12] = step1[12] + step1[19];
  output[13] = step1[13] + step1[18];
  output[14] = step1[14] + step1[17];
  output[15] = step1[15] + step1[16];
  output[16] = step1[15] - step1[16];
  output[17] = step1[14] - step1[17];
  output[18] = step1[13] - step1[18];
  output[19] = step1[12] - step1[19];
  output[20] = step1[11] - step1[20];
  output[21] = step1[10] - step1[21];
  output[22] = step1[9] - step1[22];
  output[23] = step1[8] - step1[23];
  output[24] = step1[7] - step1[24];
  output[25] = step1[6] - step1[25];
  output[26] = step1[5] - step1[26];
  output[27] = step1[4] - step1[27];
  output[28] = step1[3] - step1[28];
  output[29] = step1[2] - step1[29];
  output[30] = step1[1] - step1[30];
  output[31] = step1[0] - step1[31];
}


void vp9_short_idct32x32_c(int16_t *input, int16_t *output, int pitch) {
  int16_t out[32 * 32];
  int16_t *outptr = &out[0];
  const int short_pitch = pitch >> 1;
  int i, j;
  int16_t temp_in[32], temp_out[32];

  // First transform rows
  for (i = 0; i < 32; ++i) {
    idct32_1d(input, outptr);
    input += short_pitch;
    outptr += 32;
  }
  // Then transform columns
  for (i = 0; i < 32; ++i) {
    for (j = 0; j < 32; ++j)
      temp_in[j] = out[j * 32 + i];
    idct32_1d(temp_in, temp_out);
    for (j = 0; j < 32; ++j)
      output[j * 32 + i] = (temp_out[j] + 32) >> 6;
  }
}

void vp9_short_idct1_32x32_c(int16_t *input, int16_t *output) {
  int tmp;
  int16_t out;
  tmp = input[0] * cospi_16_64;
  out = dct_const_round_shift(tmp);
  tmp = out * cospi_16_64;
  out = dct_const_round_shift(tmp);
  *output = (out + 32) >> 6;
}

#else  // !CONFIG_DWTDCTHYBRID

#if DWT_TYPE == 53

// Note: block length must be even for this implementation
static void synthesis_53_row(int length, int16_t *lowpass, int16_t *highpass,
                             int16_t *x) {
  int16_t r, *a, *b;
  int n;

  n = length >> 1;
  b = highpass;
  a = lowpass;
  r = *highpass;
  while (n--) {
    *a++ -= (r + (*b) + 1) >> 1;
    r = *b++;
  }

  n = length >> 1;
  b = highpass;
  a = lowpass;
  while (--n) {
    *x++ = ((r = *a++) + 1) >> 1;
    *x++ = *b++ + ((r + (*a) + 2) >> 2);
  }
  *x++ = ((r = *a) + 1) >> 1;
  *x++ = *b + ((r + 1) >> 1);
}

static void synthesis_53_col(int length, int16_t *lowpass, int16_t *highpass,
                             int16_t *x) {
  int16_t r, *a, *b;
  int n;

  n = length >> 1;
  b = highpass;
  a = lowpass;
  r = *highpass;
  while (n--) {
    *a++ -= (r + (*b) + 1) >> 1;
    r = *b++;
  }

  n = length >> 1;
  b = highpass;
  a = lowpass;
  while (--n) {
    r = *a++;
    *x++ = r;
    *x++ = ((*b++) << 1) + ((r + (*a) + 1) >> 1);
  }
  *x++ = *a;
  *x++ = ((*b) << 1) + *a;
}

static void dyadic_synthesize_53(int levels, int width, int height, int16_t *c,
                                 int pitch_c, int16_t *x, int pitch_x) {
  int th[16], tw[16], lv, i, j, nh, nw, hh = height, hw = width;
  short buffer[2 * DWT_MAX_LENGTH];

  th[0] = hh;
  tw[0] = hw;
  for (i = 1; i <= levels; i++) {
    th[i] = (th[i - 1] + 1) >> 1;
    tw[i] = (tw[i - 1] + 1) >> 1;
  }
  for (lv = levels - 1; lv >= 0; lv--) {
    nh = th[lv];
    nw = tw[lv];
    hh = th[lv + 1];
    hw = tw[lv + 1];
    if ((nh < 2) || (nw < 2)) continue;
    for (j = 0; j < nw; j++) {
      for (i = 0; i < nh; i++)
        buffer[i] = c[i * pitch_c + j];
      synthesis_53_col(nh, buffer, buffer + hh, buffer + nh);
      for (i = 0; i < nh; i++)
        c[i * pitch_c + j] = buffer[i + nh];
    }
    for (i = 0; i < nh; i++) {
      memcpy(buffer, &c[i * pitch_c], nw * sizeof(*buffer));
      synthesis_53_row(nw, buffer, buffer + hw, &c[i * pitch_c]);
    }
  }
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      x[i * pitch_x + j] = c[i * pitch_c + j] >= 0 ?
          ((c[i * pitch_c + j] + DWT_PRECISION_RND) >> DWT_PRECISION_BITS) :
          -((-c[i * pitch_c + j] + DWT_PRECISION_RND) >> DWT_PRECISION_BITS);
    }
  }
}

#elif DWT_TYPE == 26

// Note: block length must be even for this implementation
static void synthesis_26_row(int length, int16_t *lowpass, int16_t *highpass,
                             int16_t *x) {
  int16_t r, s, *a, *b;
  int i, n = length >> 1;

  if (n >= 4) {
    a = lowpass;
    b = highpass;
    r = *lowpass;
    while (--n) {
      *b++ += (r - a[1] + 4) >> 3;
      r = *a++;
    }
    *b += (r - *a + 4) >> 3;
  }
  a = lowpass;
  b = highpass;
  for (i = length >> 1; i; i--) {
    s = *b++;
    r = *a++;
    *x++ = (r + s + 1) >> 1;
    *x++ = (r - s + 1) >> 1;
  }
}

static void synthesis_26_col(int length, int16_t *lowpass, int16_t *highpass,
                             int16_t *x) {
  int16_t r, s, *a, *b;
  int i, n = length >> 1;

  if (n >= 4) {
    a = lowpass;
    b = highpass;
    r = *lowpass;
    while (--n) {
      *b++ += (r - a[1] + 4) >> 3;
      r = *a++;
    }
    *b += (r - *a + 4) >> 3;
  }
  a = lowpass;
  b = highpass;
  for (i = length >> 1; i; i--) {
    s = *b++;
    r = *a++;
    *x++ = r + s;
    *x++ = r - s;
  }
}

static void dyadic_synthesize_26(int levels, int width, int height, int16_t *c,
                                 int pitch_c, int16_t *x, int pitch_x) {
  int th[16], tw[16], lv, i, j, nh, nw, hh = height, hw = width;
  int16_t buffer[2 * DWT_MAX_LENGTH];

  th[0] = hh;
  tw[0] = hw;
  for (i = 1; i <= levels; i++) {
    th[i] = (th[i - 1] + 1) >> 1;
    tw[i] = (tw[i - 1] + 1) >> 1;
  }
  for (lv = levels - 1; lv >= 0; lv--) {
    nh = th[lv];
    nw = tw[lv];
    hh = th[lv + 1];
    hw = tw[lv + 1];
    if ((nh < 2) || (nw < 2)) continue;
    for (j = 0; j < nw; j++) {
      for (i = 0; i < nh; i++)
        buffer[i] = c[i * pitch_c + j];
      synthesis_26_col(nh, buffer, buffer + hh, buffer + nh);
      for (i = 0; i < nh; i++)
        c[i * pitch_c + j] = buffer[i + nh];
    }
    for (i = 0; i < nh; i++) {
      memcpy(buffer, &c[i * pitch_c], nw * sizeof(*buffer));
      synthesis_26_row(nw, buffer, buffer + hw, &c[i * pitch_c]);
    }
  }
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      x[i * pitch_x + j] = c[i * pitch_c + j] >= 0 ?
          ((c[i * pitch_c + j] + DWT_PRECISION_RND) >> DWT_PRECISION_BITS) :
          -((-c[i * pitch_c + j] + DWT_PRECISION_RND) >> DWT_PRECISION_BITS);
    }
  }
}

#elif DWT_TYPE == 97

static void synthesis_97(int length, double *lowpass, double *highpass,
                         double *x) {
  static const double a_predict1 = -1.586134342;
  static const double a_update1 = -0.05298011854;
  static const double a_predict2 = 0.8829110762;
  static const double a_update2 = 0.4435068522;
  static const double s_low = 1.149604398;
  static const double s_high = 1/1.149604398;
  static const double inv_s_low = 1 / s_low;
  static const double inv_s_high = 1 / s_high;
  int i;
  double y[DWT_MAX_LENGTH];
  // Undo pack and scale
  for (i = 0; i < length / 2; i++) {
    y[i * 2] = lowpass[i] * inv_s_low;
    y[i * 2 + 1] = highpass[i] * inv_s_high;
  }
  memcpy(x, y, sizeof(*y) * length);
  // Undo update 2
  for (i = 2; i < length; i += 2) {
    x[i] -= a_update2 * (x[i-1] + x[i+1]);
  }
  x[0] -= 2 * a_update2 * x[1];
  // Undo predict 2
  for (i = 1; i < length - 2; i += 2) {
    x[i] -= a_predict2 * (x[i - 1] + x[i + 1]);
  }
  x[length - 1] -= 2 * a_predict2 * x[length - 2];
  // Undo update 1
  for (i = 2; i < length; i += 2) {
    x[i] -= a_update1 * (x[i - 1] + x[i + 1]);
  }
  x[0] -= 2 * a_update1 * x[1];
  // Undo predict 1
  for (i = 1; i < length - 2; i += 2) {
    x[i] -= a_predict1 * (x[i - 1] + x[i + 1]);
  }
  x[length - 1] -= 2 * a_predict1 * x[length - 2];
}

static void dyadic_synthesize_97(int levels, int width, int height, int16_t *c,
                                 int pitch_c, int16_t *x, int pitch_x) {
  int th[16], tw[16], lv, i, j, nh, nw, hh = height, hw = width;
  double buffer[2 * DWT_MAX_LENGTH];
  double y[DWT_MAX_LENGTH * DWT_MAX_LENGTH];

  th[0] = hh;
  tw[0] = hw;
  for (i = 1; i <= levels; i++) {
    th[i] = (th[i - 1] + 1) >> 1;
    tw[i] = (tw[i - 1] + 1) >> 1;
  }
  for (lv = levels - 1; lv >= 0; lv--) {
    nh = th[lv];
    nw = tw[lv];
    hh = th[lv + 1];
    hw = tw[lv + 1];
    if ((nh < 2) || (nw < 2)) continue;
    for (j = 0; j < nw; j++) {
      for (i = 0; i < nh; i++)
        buffer[i] = c[i * pitch_c + j];
      synthesis_97(nh, buffer, buffer + hh, buffer + nh);
      for (i = 0; i < nh; i++)
        y[i * DWT_MAX_LENGTH + j] = buffer[i + nh];
    }
    for (i = 0; i < nh; i++) {
      memcpy(buffer, &y[i * DWT_MAX_LENGTH], nw * sizeof(*buffer));
      synthesis_97(nw, buffer, buffer + hw, &y[i * DWT_MAX_LENGTH]);
    }
  }
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++)
      x[i * pitch_x + j] = round(y[i * DWT_MAX_LENGTH + j] /
                                 (1 << DWT_PRECISION_BITS));
}

#endif  // DWT_TYPE

// TODO(debargha): Implement scaling differently so as not to have to use the
// floating point 16x16 dct
static void butterfly_16x16_idct_1d_f(double input[16], double output[16]) {
  static const double C1 = 0.995184726672197;
  static const double C2 = 0.98078528040323;
  static const double C3 = 0.956940335732209;
  static const double C4 = 0.923879532511287;
  static const double C5 = 0.881921264348355;
  static const double C6 = 0.831469612302545;
  static const double C7 = 0.773010453362737;
  static const double C8 = 0.707106781186548;
  static const double C9 = 0.634393284163646;
  static const double C10 = 0.555570233019602;
  static const double C11 = 0.471396736825998;
  static const double C12 = 0.38268343236509;
  static const double C13 = 0.290284677254462;
  static const double C14 = 0.195090322016128;
  static const double C15 = 0.098017140329561;

  vp9_clear_system_state();  // Make it simd safe : __asm emms;
  {
    double step[16];
    double intermediate[16];
    double temp1, temp2;


    // step 1 and 2
    step[ 0] = input[0] + input[8];
    step[ 1] = input[0] - input[8];

    temp1 = input[4]*C12;
    temp2 = input[12]*C4;

    temp1 -= temp2;
    temp1 *= C8;

    step[ 2] = 2*(temp1);

    temp1 = input[4]*C4;
    temp2 = input[12]*C12;
    temp1 += temp2;
    temp1 = (temp1);
    temp1 *= C8;
    step[ 3] = 2*(temp1);

    temp1 = input[2]*C8;
    temp1 = 2*(temp1);
    temp2 = input[6] + input[10];

    step[ 4] = temp1 + temp2;
    step[ 5] = temp1 - temp2;

    temp1 = input[14]*C8;
    temp1 = 2*(temp1);
    temp2 = input[6] - input[10];

    step[ 6] = temp2 - temp1;
    step[ 7] = temp2 + temp1;

    // for odd input
    temp1 = input[3]*C12;
    temp2 = input[13]*C4;
    temp1 += temp2;
    temp1 = (temp1);
    temp1 *= C8;
    intermediate[ 8] = 2*(temp1);

    temp1 = input[3]*C4;
    temp2 = input[13]*C12;
    temp2 -= temp1;
    temp2 = (temp2);
    temp2 *= C8;
    intermediate[ 9] = 2*(temp2);

    intermediate[10] = 2*(input[9]*C8);
    intermediate[11] = input[15] - input[1];
    intermediate[12] = input[15] + input[1];
    intermediate[13] = 2*((input[7]*C8));

    temp1 = input[11]*C12;
    temp2 = input[5]*C4;
    temp2 -= temp1;
    temp2 = (temp2);
    temp2 *= C8;
    intermediate[14] = 2*(temp2);

    temp1 = input[11]*C4;
    temp2 = input[5]*C12;
    temp1 += temp2;
    temp1 = (temp1);
    temp1 *= C8;
    intermediate[15] = 2*(temp1);

    step[ 8] = intermediate[ 8] + intermediate[14];
    step[ 9] = intermediate[ 9] + intermediate[15];
    step[10] = intermediate[10] + intermediate[11];
    step[11] = intermediate[10] - intermediate[11];
    step[12] = intermediate[12] + intermediate[13];
    step[13] = intermediate[12] - intermediate[13];
    step[14] = intermediate[ 8] - intermediate[14];
    step[15] = intermediate[ 9] - intermediate[15];

    // step 3
    output[0] = step[ 0] + step[ 3];
    output[1] = step[ 1] + step[ 2];
    output[2] = step[ 1] - step[ 2];
    output[3] = step[ 0] - step[ 3];

    temp1 = step[ 4]*C14;
    temp2 = step[ 7]*C2;
    temp1 -= temp2;
    output[4] =  (temp1);

    temp1 = step[ 4]*C2;
    temp2 = step[ 7]*C14;
    temp1 += temp2;
    output[7] =  (temp1);

    temp1 = step[ 5]*C10;
    temp2 = step[ 6]*C6;
    temp1 -= temp2;
    output[5] =  (temp1);

    temp1 = step[ 5]*C6;
    temp2 = step[ 6]*C10;
    temp1 += temp2;
    output[6] =  (temp1);

    output[8] = step[ 8] + step[11];
    output[9] = step[ 9] + step[10];
    output[10] = step[ 9] - step[10];
    output[11] = step[ 8] - step[11];
    output[12] = step[12] + step[15];
    output[13] = step[13] + step[14];
    output[14] = step[13] - step[14];
    output[15] = step[12] - step[15];

    // output 4
    step[ 0] = output[0] + output[7];
    step[ 1] = output[1] + output[6];
    step[ 2] = output[2] + output[5];
    step[ 3] = output[3] + output[4];
    step[ 4] = output[3] - output[4];
    step[ 5] = output[2] - output[5];
    step[ 6] = output[1] - output[6];
    step[ 7] = output[0] - output[7];

    temp1 = output[8]*C7;
    temp2 = output[15]*C9;
    temp1 -= temp2;
    step[ 8] = (temp1);

    temp1 = output[9]*C11;
    temp2 = output[14]*C5;
    temp1 += temp2;
    step[ 9] = (temp1);

    temp1 = output[10]*C3;
    temp2 = output[13]*C13;
    temp1 -= temp2;
    step[10] = (temp1);

    temp1 = output[11]*C15;
    temp2 = output[12]*C1;
    temp1 += temp2;
    step[11] = (temp1);

    temp1 = output[11]*C1;
    temp2 = output[12]*C15;
    temp2 -= temp1;
    step[12] = (temp2);

    temp1 = output[10]*C13;
    temp2 = output[13]*C3;
    temp1 += temp2;
    step[13] = (temp1);

    temp1 = output[9]*C5;
    temp2 = output[14]*C11;
    temp2 -= temp1;
    step[14] = (temp2);

    temp1 = output[8]*C9;
    temp2 = output[15]*C7;
    temp1 += temp2;
    step[15] = (temp1);

    // step 5
    output[0] = (step[0] + step[15]);
    output[1] = (step[1] + step[14]);
    output[2] = (step[2] + step[13]);
    output[3] = (step[3] + step[12]);
    output[4] = (step[4] + step[11]);
    output[5] = (step[5] + step[10]);
    output[6] = (step[6] + step[ 9]);
    output[7] = (step[7] + step[ 8]);

    output[15] = (step[0] - step[15]);
    output[14] = (step[1] - step[14]);
    output[13] = (step[2] - step[13]);
    output[12] = (step[3] - step[12]);
    output[11] = (step[4] - step[11]);
    output[10] = (step[5] - step[10]);
    output[9] = (step[6] - step[ 9]);
    output[8] = (step[7] - step[ 8]);
  }
  vp9_clear_system_state();  // Make it simd safe : __asm emms;
}

static void vp9_short_idct16x16_c_f(int16_t *input, int16_t *output, int pitch,
                                    int scale) {
  vp9_clear_system_state();  // Make it simd safe : __asm emms;
  {
    double out[16*16], out2[16*16];
    const int short_pitch = pitch >> 1;
    int i, j;
      // First transform rows
    for (i = 0; i < 16; ++i) {
      double temp_in[16], temp_out[16];
      for (j = 0; j < 16; ++j)
        temp_in[j] = input[j + i*short_pitch];
      butterfly_16x16_idct_1d_f(temp_in, temp_out);
      for (j = 0; j < 16; ++j)
        out[j + i*16] = temp_out[j];
    }
    // Then transform columns
    for (i = 0; i < 16; ++i) {
      double temp_in[16], temp_out[16];
      for (j = 0; j < 16; ++j)
        temp_in[j] = out[j*16 + i];
      butterfly_16x16_idct_1d_f(temp_in, temp_out);
      for (j = 0; j < 16; ++j)
        out2[j*16 + i] = temp_out[j];
    }
    for (i = 0; i < 16*16; ++i)
      output[i] = round(out2[i] / (128 >> scale));
  }
  vp9_clear_system_state();  // Make it simd safe : __asm emms;
}

static void idct8_1d_f(double *x) {
  int i, j;
  double t[8];
  static const double idctmat[64] = {
    0.35355339059327,  0.49039264020162,  0.46193976625564,  0.41573480615127,
    0.35355339059327,   0.2777851165098,  0.19134171618254, 0.097545161008064,
    0.35355339059327,  0.41573480615127,  0.19134171618254, -0.097545161008064,
    -0.35355339059327, -0.49039264020161, -0.46193976625564,  -0.2777851165098,
    0.35355339059327,   0.2777851165098, -0.19134171618254, -0.49039264020162,
    -0.35355339059327, 0.097545161008064,  0.46193976625564,  0.41573480615127,
    0.35355339059327, 0.097545161008063, -0.46193976625564,  -0.2777851165098,
    0.35355339059327,  0.41573480615127, -0.19134171618254, -0.49039264020162,
    0.35355339059327, -0.097545161008063, -0.46193976625564,   0.2777851165098,
    0.35355339059327, -0.41573480615127, -0.19134171618255,  0.49039264020162,
    0.35355339059327,  -0.2777851165098, -0.19134171618254,  0.49039264020161,
    -0.35355339059327, -0.097545161008064,  0.46193976625564, -0.41573480615127,
    0.35355339059327, -0.41573480615127,  0.19134171618254, 0.097545161008065,
    -0.35355339059327,  0.49039264020162, -0.46193976625564,   0.2777851165098,
    0.35355339059327, -0.49039264020162,  0.46193976625564, -0.41573480615127,
    0.35355339059327,  -0.2777851165098,  0.19134171618255, -0.097545161008064
  };
  for (i = 0; i < 8; ++i) {
    t[i] = 0;
    for (j = 0; j < 8; ++j)
      t[i] += idctmat[i * 8 + j] * x[j];
  }
  for (i = 0; i < 8; ++i) {
    x[i] = t[i];
  }
}

static void vp9_short_idct8x8_c_f(int16_t *coefs, int16_t *block, int pitch,
                                  int scale) {
  double X[8 * 8], Y[8];
  int i, j;
  int shortpitch = pitch >> 1;

  vp9_clear_system_state();  // Make it simd safe : __asm emms;
  {
    for (i = 0; i < 8; i++) {
      for (j = 0; j < 8; j++) {
        X[i * 8 + j] = (double)coefs[i * shortpitch + j];
      }
    }
    for (i = 0; i < 8; i++)
      idct8_1d_f(X + 8 * i);
    for (i = 0; i < 8; i++) {
      for (j = 0; j < 8; ++j)
        Y[j] = X[i + 8 * j];
      idct8_1d_f(Y);
      for (j = 0; j < 8; ++j)
        X[i + 8 * j] = Y[j];
    }
    for (i = 0; i < 8; i++) {
      for (j = 0; j < 8; j++) {
        block[i * 8 + j] = (int16_t)round(X[i * 8 + j] / (8 >> scale));
      }
    }
  }
  vp9_clear_system_state();  // Make it simd safe : __asm emms;
}

#define multiply_bits(d, n) ((n) < 0 ? (d) >> (n) : (d) << (n))

#if DWTDCT_TYPE == DWTDCT16X16_LEAN

void vp9_short_idct32x32_c(int16_t *input, int16_t *output, int pitch) {
  // assume output is a 32x32 buffer
  // Temporary buffer to hold a 16x16 block for 16x16 inverse dct
  int16_t buffer[16 * 16];
  // Temporary buffer to hold a 32x32 block for inverse 32x32 dwt
  int16_t buffer2[32 * 32];
  // Note: pitch is in bytes, short_pitch is in short units
  const int short_pitch = pitch >> 1;
  int i, j;

  // TODO(debargha): Implement more efficiently by adding output pitch
  // argument to the idct16x16 function
  vp9_short_idct16x16_c_f(input, buffer, pitch,
                          1 + DWT_PRECISION_BITS);
  for (i = 0; i < 16; ++i) {
    vpx_memcpy(buffer2 + i * 32, buffer + i * 16, sizeof(*buffer2) * 16);
  }
  for (i = 0; i < 16; ++i) {
    for (j = 16; j < 32; ++j) {
      buffer2[i * 32 + j] =
          multiply_bits(input[i * short_pitch + j], DWT_PRECISION_BITS - 2);
    }
  }
  for (i = 16; i < 32; ++i) {
    for (j = 0; j < 32; ++j) {
      buffer2[i * 32 + j] =
          multiply_bits(input[i * short_pitch + j], DWT_PRECISION_BITS - 2);
    }
  }
#if DWT_TYPE == 26
  dyadic_synthesize_26(1, 32, 32, buffer2, 32, output, 32);
#elif DWT_TYPE == 97
  dyadic_synthesize_97(1, 32, 32, buffer2, 32, output, 32);
#elif DWT_TYPE == 53
  dyadic_synthesize_53(1, 32, 32, buffer2, 32, output, 32);
#endif
}

#elif DWTDCT_TYPE == DWTDCT16X16

void vp9_short_idct32x32_c(int16_t *input, int16_t *output, int pitch) {
  // assume output is a 32x32 buffer
  // Temporary buffer to hold a 16x16 block for 16x16 inverse dct
  int16_t buffer[16 * 16];
  // Temporary buffer to hold a 32x32 block for inverse 32x32 dwt
  int16_t buffer2[32 * 32];
  // Note: pitch is in bytes, short_pitch is in short units
  const int short_pitch = pitch >> 1;
  int i, j;

  // TODO(debargha): Implement more efficiently by adding output pitch
  // argument to the idct16x16 function
  vp9_short_idct16x16_c_f(input, buffer, pitch,
                          1 + DWT_PRECISION_BITS);
  for (i = 0; i < 16; ++i) {
    vpx_memcpy(buffer2 + i * 32, buffer + i * 16, sizeof(*buffer2) * 16);
  }
  vp9_short_idct16x16_c_f(input + 16, buffer, pitch,
                          1 + DWT_PRECISION_BITS);
  for (i = 0; i < 16; ++i) {
    vpx_memcpy(buffer2 + i * 32 + 16, buffer + i * 16, sizeof(*buffer2) * 16);
  }
  vp9_short_idct16x16_c_f(input + 16 * short_pitch, buffer, pitch,
                          1 + DWT_PRECISION_BITS);
  for (i = 0; i < 16; ++i) {
    vpx_memcpy(buffer2 + i * 32 + 16 * 32, buffer + i * 16,
               sizeof(*buffer2) * 16);
  }
  vp9_short_idct16x16_c_f(input + 16 * short_pitch + 16, buffer, pitch,
                          1 + DWT_PRECISION_BITS);
  for (i = 0; i < 16; ++i) {
    vpx_memcpy(buffer2 + i * 32 + 16 * 33, buffer + i * 16,
               sizeof(*buffer2) * 16);
  }
#if DWT_TYPE == 26
  dyadic_synthesize_26(1, 32, 32, buffer2, 32, output, 32);
#elif DWT_TYPE == 97
  dyadic_synthesize_97(1, 32, 32, buffer2, 32, output, 32);
#elif DWT_TYPE == 53
  dyadic_synthesize_53(1, 32, 32, buffer2, 32, output, 32);
#endif
}

#elif DWTDCT_TYPE == DWTDCT8X8

void vp9_short_idct32x32_c(int16_t *input, int16_t *output, int pitch) {
  // assume output is a 32x32 buffer
  // Temporary buffer to hold a 16x16 block for 16x16 inverse dct
  int16_t buffer[8 * 8];
  // Temporary buffer to hold a 32x32 block for inverse 32x32 dwt
  int16_t buffer2[32 * 32];
  // Note: pitch is in bytes, short_pitch is in short units
  const int short_pitch = pitch >> 1;
  int i, j;

  // TODO(debargha): Implement more efficiently by adding output pitch
  // argument to the idct16x16 function
  vp9_short_idct8x8_c_f(input, buffer, pitch,
                        1 + DWT_PRECISION_BITS);
  for (i = 0; i < 8; ++i) {
    vpx_memcpy(buffer2 + i * 32, buffer + i * 8, sizeof(*buffer2) * 8);
  }
  vp9_short_idct8x8_c_f(input + 8, buffer, pitch,
                        1 + DWT_PRECISION_BITS);
  for (i = 0; i < 8; ++i) {
    vpx_memcpy(buffer2 + i * 32 + 8, buffer + i * 8, sizeof(*buffer2) * 8);
  }
  vp9_short_idct8x8_c_f(input + 8 * short_pitch, buffer, pitch,
                        1 + DWT_PRECISION_BITS);
  for (i = 0; i < 8; ++i) {
    vpx_memcpy(buffer2 + i * 32 + 8 * 32, buffer + i * 8,
               sizeof(*buffer2) * 8);
  }
  vp9_short_idct8x8_c_f(input + 8 * short_pitch + 8, buffer, pitch,
                        1 + DWT_PRECISION_BITS);
  for (i = 0; i < 8; ++i) {
    vpx_memcpy(buffer2 + i * 32 + 8 * 33, buffer + i * 8,
               sizeof(*buffer2) * 8);
  }
  for (i = 0; i < 16; ++i) {
    for (j = 16; j < 32; ++j) {
      buffer2[i * 32 + j] =
          multiply_bits(input[i * short_pitch + j], DWT_PRECISION_BITS - 2);
    }
  }
  for (i = 16; i < 32; ++i) {
    for (j = 0; j < 32; ++j) {
      buffer2[i * 32 + j] =
          multiply_bits(input[i * short_pitch + j], DWT_PRECISION_BITS - 2);
    }
  }
#if DWT_TYPE == 26
  dyadic_synthesize_26(2, 32, 32, buffer2, 32, output, 32);
#elif DWT_TYPE == 97
  dyadic_synthesize_97(2, 32, 32, buffer2, 32, output, 32);
#elif DWT_TYPE == 53
  dyadic_synthesize_53(2, 32, 32, buffer2, 32, output, 32);
#endif
}

#endif

#if CONFIG_TX64X64
void vp9_short_idct64x64_c(int16_t *input, int16_t *output, int pitch) {
  // assume output is a 64x64 buffer
  // Temporary buffer to hold a 16x16 block for 16x16 inverse dct
  int16_t buffer[16 * 16];
  // Temporary buffer to hold a 32x32 block for inverse 32x32 dwt
  int16_t buffer2[64 * 64];
  // Note: pitch is in bytes, short_pitch is in short units
  const int short_pitch = pitch >> 1;
  int i, j;

  // TODO(debargha): Implement more efficiently by adding output pitch
  // argument to the idct16x16 function
  vp9_short_idct16x16_c_f(input, buffer, pitch,
                          2 + DWT_PRECISION_BITS);
  for (i = 0; i < 16; ++i) {
    vpx_memcpy(buffer2 + i * 64, buffer + i * 16, sizeof(*buffer2) * 16);
  }
#if DWTDCT_TYPE == DWTDCT16X16_LEAN
  for (i = 0; i < 16; ++i) {
    for (j = 16; j < 64; ++j) {
      buffer2[i * 64 + j] =
          multiply_bits(input[i * short_pitch + j], DWT_PRECISION_BITS - 1);
    }
  }
  for (i = 16; i < 64; ++i) {
    for (j = 0; j < 64; ++j) {
      buffer2[i * 64 + j] =
          multiply_bits(input[i * short_pitch + j], DWT_PRECISION_BITS - 1);
    }
  }
#elif DWTDCT_TYPE == DWTDCT16X16
  vp9_short_idct16x16_c_f(input + 16, buffer, pitch,
                          2 + DWT_PRECISION_BITS);
  for (i = 0; i < 16; ++i) {
    vpx_memcpy(buffer2 + i * 64 + 16, buffer + i * 16, sizeof(*buffer2) * 16);
  }
  vp9_short_idct16x16_c_f(input + 16 * short_pitch, buffer, pitch,
                          2 + DWT_PRECISION_BITS);
  for (i = 0; i < 16; ++i) {
    vpx_memcpy(buffer2 + i * 64 + 16 * 64, buffer + i * 16,
               sizeof(*buffer2) * 16);
  }
  vp9_short_idct16x16_c_f(input + 16 * short_pitch + 16, buffer, pitch,
                          2 + DWT_PRECISION_BITS);
  for (i = 0; i < 16; ++i) {
    vpx_memcpy(buffer2 + i * 64 + 16 * 65, buffer + i * 16,
               sizeof(*buffer2) * 16);
  }

  // Copying and scaling highest bands into buffer2
  for (i = 0; i < 32; ++i) {
    for (j = 32; j < 64; ++j) {
      buffer2[i * 64 + j] =
          multiply_bits(input[i * short_pitch + j], DWT_PRECISION_BITS - 1);
    }
  }
  for (i = 32; i < 64; ++i) {
    for (j = 0; j < 64; ++j) {
      buffer2[i * 64 + j] =
          multiply_bits(input[i * short_pitch + j], DWT_PRECISION_BITS - 1);
    }
  }
#endif  // DWTDCT_TYPE

#if DWT_TYPE == 26
  dyadic_synthesize_26(2, 64, 64, buffer2, 64, output, 64);
#elif DWT_TYPE == 97
  dyadic_synthesize_97(2, 64, 64, buffer2, 64, output, 64);
#elif DWT_TYPE == 53
  dyadic_synthesize_53(2, 64, 64, buffer2, 64, output, 64);
#endif
}
#endif  // CONFIG_TX64X64
#endif  // !CONFIG_DWTDCTHYBRID
