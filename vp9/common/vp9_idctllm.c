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
#include "vpx_ports/config.h"
#include "vp9/common/vp9_systemdependent.h"

#include "vp9/common/vp9_blockd.h"

static const int cospi8sqrt2minus1 = 20091;
static const int sinpi8sqrt2      = 35468;
static const int rounding = 0;

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


/* Converted the transforms to integer form. */
#define VERTICAL_SHIFT 14  // 16
#define VERTICAL_ROUNDING ((1 << (VERTICAL_SHIFT - 1)) - 1)
#define HORIZONTAL_SHIFT 17  // 15
#define HORIZONTAL_ROUNDING ((1 << (HORIZONTAL_SHIFT - 1)) - 1)
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

  /* vertical transformation */
  for (j = 0; j < tx_dim; j++) {
    for (i = 0; i < nz_dim; i++) {
      int temp = 0;

      for (k = 0; k < nz_dim; k++) {
        temp += ptv[k] * ip[(k * tx_dim)];
      }

      im[i] = (int16_t)((temp + VERTICAL_ROUNDING) >> VERTICAL_SHIFT);
      ip++;
    }
    im += tx_dim;  // 16
    ptv += tx_dim;
    ip = input;
  }

  /* horizontal transformation */
  im = &imbuf[0];

  for (j = 0; j < tx_dim; j++) {
    const int16_t *pthc = pth;

    for (i = 0; i < tx_dim; i++) {
      int temp = 0;

      for (k = 0; k < nz_dim; k++) {
        temp += im[k] * pthc[k];
      }

      op[i] = (int16_t)((temp + HORIZONTAL_ROUNDING) >> HORIZONTAL_SHIFT);
      pthc += tx_dim;
    }

    im += tx_dim;  // 16
    op += shortpitch;
  }
}

void vp9_short_idct4x4llm_c(short *input, short *output, int pitch) {
  int i;
  int a1, b1, c1, d1;

  short *ip = input;
  short *op = output;
  int temp1, temp2;
  int shortpitch = pitch >> 1;

  for (i = 0; i < 4; i++) {
    a1 = ip[0] + ip[8];
    b1 = ip[0] - ip[8];

    temp1 = (ip[4] * sinpi8sqrt2 + rounding) >> 16;
    temp2 = ip[12] + ((ip[12] * cospi8sqrt2minus1 + rounding) >> 16);
    c1 = temp1 - temp2;

    temp1 = ip[4] + ((ip[4] * cospi8sqrt2minus1 + rounding) >> 16);
    temp2 = (ip[12] * sinpi8sqrt2 + rounding) >> 16;
    d1 = temp1 + temp2;

    op[shortpitch * 0] = a1 + d1;
    op[shortpitch * 3] = a1 - d1;

    op[shortpitch * 1] = b1 + c1;
    op[shortpitch * 2] = b1 - c1;

    ip++;
    op++;
  }

  ip = output;
  op = output;

  for (i = 0; i < 4; i++) {
    a1 = ip[0] + ip[2];
    b1 = ip[0] - ip[2];

    temp1 = (ip[1] * sinpi8sqrt2 + rounding) >> 16;
    temp2 = ip[3] + ((ip[3] * cospi8sqrt2minus1 + rounding) >> 16);
    c1 = temp1 - temp2;

    temp1 = ip[1] + ((ip[1] * cospi8sqrt2minus1 + rounding) >> 16);
    temp2 = (ip[3] * sinpi8sqrt2 + rounding) >> 16;
    d1 = temp1 + temp2;

    op[0] = (a1 + d1 + 16) >> 5;
    op[3] = (a1 - d1 + 16) >> 5;

    op[1] = (b1 + c1 + 16) >> 5;
    op[2] = (b1 - c1 + 16) >> 5;

    ip += shortpitch;
    op += shortpitch;
  }
}

void vp9_short_idct4x4llm_1_c(short *input, short *output, int pitch) {
  int i;
  int a1;
  short *op = output;
  int shortpitch = pitch >> 1;
  a1 = ((input[0] + 16) >> 5);
  for (i = 0; i < 4; i++) {
    op[0] = a1;
    op[1] = a1;
    op[2] = a1;
    op[3] = a1;
    op += shortpitch;
  }
}

void vp9_dc_only_idct_add_c(short input_dc, unsigned char *pred_ptr,
                            unsigned char *dst_ptr, int pitch, int stride) {
  int a1 = ((input_dc + 16) >> 5);
  int r, c;

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      int a = a1 + pred_ptr[c];

      if (a < 0)
        a = 0;

      if (a > 255)
        a = 255;

      dst_ptr[c] = (unsigned char) a;
    }

    dst_ptr += stride;
    pred_ptr += pitch;
  }
}

void vp9_short_inv_walsh4x4_c(short *input, short *output) {
  int i;
  int a1, b1, c1, d1;
  short *ip = input;
  short *op = output;

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

void vp9_short_inv_walsh4x4_1_c(short *in, short *out) {
  int i;
  short tmp[4];
  short *ip = in;
  short *op = tmp;

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
void vp9_short_inv_walsh4x4_lossless_c(short *input, short *output) {
  int i;
  int a1, b1, c1, d1;
  short *ip = input;
  short *op = output;

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

void vp9_short_inv_walsh4x4_1_lossless_c(short *in, short *out) {
  int i;
  short tmp[4];
  short *ip = in;
  short *op = tmp;

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

void vp9_short_inv_walsh4x4_x8_c(short *input, short *output, int pitch) {
  int i;
  int a1, b1, c1, d1;
  short *ip = input;
  short *op = output;
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

void vp9_short_inv_walsh4x4_1_x8_c(short *in, short *out, int pitch) {
  int i;
  short tmp[4];
  short *ip = in;
  short *op = tmp;
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

void vp9_dc_only_inv_walsh_add_c(short input_dc, unsigned char *pred_ptr,
                                 unsigned char *dst_ptr,
                                 int pitch, int stride) {
  int r, c;
  short tmp[16];
  vp9_short_inv_walsh4x4_1_x8_c(&input_dc, tmp, 4 << 1);

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      int a = tmp[r * 4 + c] + pred_ptr[c];
      if (a < 0)
        a = 0;

      if (a > 255)
        a = 255;

      dst_ptr[c] = (unsigned char) a;
    }

    dst_ptr += stride;
    pred_ptr += pitch;
  }
}
#endif

void vp9_dc_only_idct_add_8x8_c(short input_dc,
                                unsigned char *pred_ptr,
                                unsigned char *dst_ptr,
                                int pitch, int stride) {
  int a1 = ((input_dc + 16) >> 5);
  int r, c, b;
  unsigned char *orig_pred = pred_ptr;
  unsigned char *orig_dst = dst_ptr;
  for (b = 0; b < 4; b++) {
    for (r = 0; r < 4; r++) {
      for (c = 0; c < 4; c++) {
        int a = a1 + pred_ptr[c];

        if (a < 0)
          a = 0;

        if (a > 255)
          a = 255;

        dst_ptr[c] = (unsigned char) a;
      }

      dst_ptr += stride;
      pred_ptr += pitch;
    }
    dst_ptr = orig_dst + (b + 1) % 2 * 4 + (b + 1) / 2 * 4 * stride;
    pred_ptr = orig_pred + (b + 1) % 2 * 4 + (b + 1) / 2 * 4 * pitch;
  }
}

#define W1 2841                 /* 2048*sqrt(2)*cos(1*pi/16) */
#define W2 2676                 /* 2048*sqrt(2)*cos(2*pi/16) */
#define W3 2408                 /* 2048*sqrt(2)*cos(3*pi/16) */
#define W5 1609                 /* 2048*sqrt(2)*cos(5*pi/16) */
#define W6 1108                 /* 2048*sqrt(2)*cos(6*pi/16) */
#define W7 565                  /* 2048*sqrt(2)*cos(7*pi/16) */

/* row (horizontal) IDCT
 *
 * 7                       pi         1 dst[k] = sum c[l] * src[l] * cos( -- *
 * ( k + - ) * l ) l=0                      8          2
 *
 * where: c[0]    = 128 c[1..7] = 128*sqrt(2) */

static void idctrow(int *blk) {
  int x0, x1, x2, x3, x4, x5, x6, x7, x8;
  /* shortcut */
  if (!((x1 = blk[4] << 11) | (x2 = blk[6]) | (x3 = blk[2]) |
        (x4 = blk[1]) | (x5 = blk[7]) | (x6 = blk[5]) | (x7 = blk[3]))) {
    blk[0] = blk[1] = blk[2] = blk[3] = blk[4]
                                        = blk[5] = blk[6] = blk[7] = blk[0] << 3;
    return;
  }

  x0 = (blk[0] << 11) + 128;    /* for proper rounding in the fourth stage */
  /* first stage */
  x8 = W7 * (x4 + x5);
  x4 = x8 + (W1 - W7) * x4;
  x5 = x8 - (W1 + W7) * x5;
  x8 = W3 * (x6 + x7);
  x6 = x8 - (W3 - W5) * x6;
  x7 = x8 - (W3 + W5) * x7;

  /* second stage */
  x8 = x0 + x1;
  x0 -= x1;
  x1 = W6 * (x3 + x2);
  x2 = x1 - (W2 + W6) * x2;
  x3 = x1 + (W2 - W6) * x3;
  x1 = x4 + x6;
  x4 -= x6;
  x6 = x5 + x7;
  x5 -= x7;

  /* third stage */
  x7 = x8 + x3;
  x8 -= x3;
  x3 = x0 + x2;
  x0 -= x2;
  x2 = (181 * (x4 + x5) + 128) >> 8;
  x4 = (181 * (x4 - x5) + 128) >> 8;

  /* fourth stage */
  blk[0] = (x7 + x1) >> 8;
  blk[1] = (x3 + x2) >> 8;
  blk[2] = (x0 + x4) >> 8;
  blk[3] = (x8 + x6) >> 8;
  blk[4] = (x8 - x6) >> 8;
  blk[5] = (x0 - x4) >> 8;
  blk[6] = (x3 - x2) >> 8;
  blk[7] = (x7 - x1) >> 8;
}

/* column (vertical) IDCT
 *
 * 7                         pi         1 dst[8*k] = sum c[l] * src[8*l] *
 * cos( -- * ( k + - ) * l ) l=0                        8          2
 *
 * where: c[0]    = 1/1024 c[1..7] = (1/1024)*sqrt(2) */
static void idctcol(int *blk) {
  int x0, x1, x2, x3, x4, x5, x6, x7, x8;

  /* shortcut */
  if (!((x1 = (blk[8 * 4] << 8)) | (x2 = blk[8 * 6]) | (x3 = blk[8 * 2]) |
        (x4 = blk[8 * 1]) | (x5 = blk[8 * 7]) | (x6 = blk[8 * 5]) |
        (x7 = blk[8 * 3]))) {
    blk[8 * 0] = blk[8 * 1] = blk[8 * 2] = blk[8 * 3]
                                           = blk[8 * 4] = blk[8 * 5] = blk[8 * 6]
                                                                       = blk[8 * 7] = ((blk[8 * 0] + 32) >> 6);
    return;
  }

  x0 = (blk[8 * 0] << 8) + 16384;

  /* first stage */
  x8 = W7 * (x4 + x5) + 4;
  x4 = (x8 + (W1 - W7) * x4) >> 3;
  x5 = (x8 - (W1 + W7) * x5) >> 3;
  x8 = W3 * (x6 + x7) + 4;
  x6 = (x8 - (W3 - W5) * x6) >> 3;
  x7 = (x8 - (W3 + W5) * x7) >> 3;

  /* second stage */
  x8 = x0 + x1;
  x0 -= x1;
  x1 = W6 * (x3 + x2) + 4;
  x2 = (x1 - (W2 + W6) * x2) >> 3;
  x3 = (x1 + (W2 - W6) * x3) >> 3;
  x1 = x4 + x6;
  x4 -= x6;
  x6 = x5 + x7;
  x5 -= x7;

  /* third stage */
  x7 = x8 + x3;
  x8 -= x3;
  x3 = x0 + x2;
  x0 -= x2;
  x2 = (181 * (x4 + x5) + 128) >> 8;
  x4 = (181 * (x4 - x5) + 128) >> 8;

  /* fourth stage */
  blk[8 * 0] = (x7 + x1) >> 14;
  blk[8 * 1] = (x3 + x2) >> 14;
  blk[8 * 2] = (x0 + x4) >> 14;
  blk[8 * 3] = (x8 + x6) >> 14;
  blk[8 * 4] = (x8 - x6) >> 14;
  blk[8 * 5] = (x0 - x4) >> 14;
  blk[8 * 6] = (x3 - x2) >> 14;
  blk[8 * 7] = (x7 - x1) >> 14;
}

#define TX_DIM 8
void vp9_short_idct8x8_c(short *coefs, short *block, int pitch) {
  int X[TX_DIM * TX_DIM];
  int i, j;
  int shortpitch = pitch >> 1;

  for (i = 0; i < TX_DIM; i++) {
    for (j = 0; j < TX_DIM; j++) {
      X[i * TX_DIM + j] = (int)(coefs[i * TX_DIM + j] + 1
                                + (coefs[i * TX_DIM + j] < 0)) >> 2;
    }
  }
  for (i = 0; i < 8; i++)
    idctrow(X + 8 * i);

  for (i = 0; i < 8; i++)
    idctcol(X + i);

  for (i = 0; i < TX_DIM; i++) {
    for (j = 0; j < TX_DIM; j++) {
      block[i * shortpitch + j]  = X[i * TX_DIM + j] >> 1;
    }
  }
}

/* Row IDCT when only first 4 coefficients are non-zero. */
static void idctrow10(int *blk) {
  int x0, x1, x2, x3, x4, x5, x6, x7, x8;

  /* shortcut */
  if (!((x1 = blk[4] << 11) | (x2 = blk[6]) | (x3 = blk[2]) |
        (x4 = blk[1]) | (x5 = blk[7]) | (x6 = blk[5]) | (x7 = blk[3]))) {
    blk[0] = blk[1] = blk[2] = blk[3] = blk[4]
           = blk[5] = blk[6] = blk[7] = blk[0] << 3;
    return;
  }

  x0 = (blk[0] << 11) + 128;    /* for proper rounding in the fourth stage */
  /* first stage */
  x5 = W7 * x4;
  x4 = W1 * x4;
  x6 = W3 * x7;
  x7 = -W5 * x7;

  /* second stage */
  x2 = W6 * x3;
  x3 = W2 * x3;
  x1 = x4 + x6;
  x4 -= x6;
  x6 = x5 + x7;
  x5 -= x7;

  /* third stage */
  x7 = x0 + x3;
  x8 = x0 - x3;
  x3 = x0 + x2;
  x0 -= x2;
  x2 = (181 * (x4 + x5) + 128) >> 8;
  x4 = (181 * (x4 - x5) + 128) >> 8;

  /* fourth stage */
  blk[0] = (x7 + x1) >> 8;
  blk[1] = (x3 + x2) >> 8;
  blk[2] = (x0 + x4) >> 8;
  blk[3] = (x8 + x6) >> 8;
  blk[4] = (x8 - x6) >> 8;
  blk[5] = (x0 - x4) >> 8;
  blk[6] = (x3 - x2) >> 8;
  blk[7] = (x7 - x1) >> 8;
}

/* Column (vertical) IDCT when only first 4 coefficients are non-zero. */
static void idctcol10(int *blk) {
  int x0, x1, x2, x3, x4, x5, x6, x7, x8;

  /* shortcut */
  if (!((x1 = (blk[8 * 4] << 8)) | (x2 = blk[8 * 6]) | (x3 = blk[8 * 2]) |
        (x4 = blk[8 * 1]) | (x5 = blk[8 * 7]) | (x6 = blk[8 * 5]) |
        (x7 = blk[8 * 3]))) {
    blk[8 * 0] = blk[8 * 1] = blk[8 * 2] = blk[8 * 3]
        = blk[8 * 4] = blk[8 * 5] = blk[8 * 6]
        = blk[8 * 7] = ((blk[8 * 0] + 32) >> 6);
    return;
  }

  x0 = (blk[8 * 0] << 8) + 16384;

  /* first stage */
  x5 = (W7 * x4 + 4) >> 3;
  x4 = (W1 * x4 + 4) >> 3;
  x6 = (W3 * x7 + 4) >> 3;
  x7 = (-W5 * x7 + 4) >> 3;

  /* second stage */
  x2 = (W6 * x3 + 4) >> 3;
  x3 = (W2 * x3 + 4) >> 3;
  x1 = x4 + x6;
  x4 -= x6;
  x6 = x5 + x7;
  x5 -= x7;

  /* third stage */
  x7 = x0 + x3;
  x8 = x0 - x3;
  x3 = x0 + x2;
  x0 -= x2;
  x2 = (181 * (x4 + x5) + 128) >> 8;
  x4 = (181 * (x4 - x5) + 128) >> 8;

  /* fourth stage */
  blk[8 * 0] = (x7 + x1) >> 14;
  blk[8 * 1] = (x3 + x2) >> 14;
  blk[8 * 2] = (x0 + x4) >> 14;
  blk[8 * 3] = (x8 + x6) >> 14;
  blk[8 * 4] = (x8 - x6) >> 14;
  blk[8 * 5] = (x0 - x4) >> 14;
  blk[8 * 6] = (x3 - x2) >> 14;
  blk[8 * 7] = (x7 - x1) >> 14;
}

void vp9_short_idct10_8x8_c(short *coefs, short *block, int pitch) {
  int X[TX_DIM * TX_DIM];
  int i, j;
  int shortpitch = pitch >> 1;

  for (i = 0; i < TX_DIM; i++) {
    for (j = 0; j < TX_DIM; j++) {
      X[i * TX_DIM + j] = (int)(coefs[i * TX_DIM + j] + 1
                                + (coefs[i * TX_DIM + j] < 0)) >> 2;
    }
  }

  /* Do first 4 row idct only since non-zero dct coefficients are all in
   *  upper-left 4x4 area. */
  for (i = 0; i < 4; i++)
    idctrow10(X + 8 * i);

  for (i = 0; i < 8; i++)
    idctcol10(X + i);

  for (i = 0; i < TX_DIM; i++) {
    for (j = 0; j < TX_DIM; j++) {
      block[i * shortpitch + j]  = X[i * TX_DIM + j] >> 1;
    }
  }
}

void vp9_short_ihaar2x2_c(short *input, short *output, int pitch) {
  int i;
  short *ip = input; // 0,1, 4, 8
  short *op = output;
  for (i = 0; i < 16; i++) {
    op[i] = 0;
  }

  op[0] = (ip[0] + ip[1] + ip[4] + ip[8] + 1) >> 1;
  op[1] = (ip[0] - ip[1] + ip[4] - ip[8]) >> 1;
  op[4] = (ip[0] + ip[1] - ip[4] - ip[8]) >> 1;
  op[8] = (ip[0] - ip[1] - ip[4] + ip[8]) >> 1;
}


#if 0
// Keep a really bad float version as reference for now.
void vp9_short_idct16x16_c(short *input, short *output, int pitch) {

  vp9_clear_system_state(); // Make it simd safe : __asm emms;
  {
    double x;
    const int short_pitch = pitch >> 1;
    int i, j, k, l;
    for (l = 0; l < 16; ++l) {
      for (k = 0; k < 16; ++k) {
        double s = 0;
        for (i = 0; i < 16; ++i) {
          for (j = 0; j < 16; ++j) {
            x=cos(PI*j*(l+0.5)/16.0)*cos(PI*i*(k+0.5)/16.0)*input[i*16+j]/32;
            if (i != 0)
              x *= sqrt(2.0);
            if (j != 0)
              x *= sqrt(2.0);
            s += x;
          }
        }
        output[k*short_pitch+l] = (short)round(s);
      }
    }
  }
  vp9_clear_system_state(); // Make it simd safe : __asm emms;
}
#endif

#define TEST_INT_16x16_IDCT 1
#if !TEST_INT_16x16_IDCT
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


static void butterfly_16x16_idct_1d(double input[16], double output[16]) {

  vp9_clear_system_state(); // Make it simd safe : __asm emms;
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
  vp9_clear_system_state(); // Make it simd safe : __asm emms;
}

// Remove once an int version of iDCT is written
#if 0
void reference_16x16_idct_1d(double input[16], double output[16]) {

  vp9_clear_system_state(); // Make it simd safe : __asm emms;
  {
    const double kPi = 3.141592653589793238462643383279502884;
    const double kSqrt2 = 1.414213562373095048801688724209698;
    for (int k = 0; k < 16; k++) {
      output[k] = 0.0;
      for (int n = 0; n < 16; n++) {
        output[k] += input[n]*cos(kPi*(2*k+1)*n/32.0);
        if (n == 0)
          output[k] = output[k]/kSqrt2;
      }
    }
  }
  vp9_clear_system_state(); // Make it simd safe : __asm emms;
}
#endif

void vp9_short_idct16x16_c(short *input, short *output, int pitch) {

  vp9_clear_system_state(); // Make it simd safe : __asm emms;
  {
    double out[16*16], out2[16*16];
    const int short_pitch = pitch >> 1;
    int i, j;
      // First transform rows
    for (i = 0; i < 16; ++i) {
      double temp_in[16], temp_out[16];
      for (j = 0; j < 16; ++j)
        temp_in[j] = input[j + i*short_pitch];
      butterfly_16x16_idct_1d(temp_in, temp_out);
      for (j = 0; j < 16; ++j)
        out[j + i*16] = temp_out[j];
    }
    // Then transform columns
    for (i = 0; i < 16; ++i) {
      double temp_in[16], temp_out[16];
      for (j = 0; j < 16; ++j)
        temp_in[j] = out[j*16 + i];
      butterfly_16x16_idct_1d(temp_in, temp_out);
      for (j = 0; j < 16; ++j)
        out2[j*16 + i] = temp_out[j];
    }
    for (i = 0; i < 16*16; ++i)
      output[i] = round(out2[i]/128);
  }
  vp9_clear_system_state(); // Make it simd safe : __asm emms;
}

#else
static const int16_t C1 = 16305;
static const int16_t C2 = 16069;
static const int16_t C3 = 15679;
static const int16_t C4 = 15137;
static const int16_t C5 = 14449;
static const int16_t C6 = 13623;
static const int16_t C7 = 12665;
static const int16_t C8 = 11585;
static const int16_t C9 = 10394;
static const int16_t C10 = 9102;
static const int16_t C11 = 7723;
static const int16_t C12 = 6270;
static const int16_t C13 = 4756;
static const int16_t C14 = 3196;
static const int16_t C15 = 1606;

#define INITIAL_SHIFT 2
#define INITIAL_ROUNDING (1 << (INITIAL_SHIFT - 1))
#define RIGHT_SHIFT 14
#define RIGHT_ROUNDING (1 << (RIGHT_SHIFT - 1))

static void butterfly_16x16_idct_1d(int16_t input[16], int16_t output[16],
                                    int last_shift_bits) {
    int16_t step[16];
    int intermediate[16];
    int temp1, temp2;

    int step1_shift = RIGHT_SHIFT + INITIAL_SHIFT;
    int step1_rounding = 1 << (step1_shift - 1);
    int last_rounding = 0;

    if (last_shift_bits > 0)
      last_rounding = 1 << (last_shift_bits - 1);

    // step 1 and 2
    step[ 0] = (input[0] + input[8] + INITIAL_ROUNDING) >> INITIAL_SHIFT;
    step[ 1] = (input[0] - input[8] + INITIAL_ROUNDING) >> INITIAL_SHIFT;

    temp1 = input[4] * C12;
    temp2 = input[12] * C4;
    temp1 = (temp1 - temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;
    temp1  *= C8;
    step[ 2] = (2 * (temp1) + step1_rounding) >> step1_shift;

    temp1 = input[4] * C4;
    temp2 = input[12] * C12;
    temp1 = (temp1 + temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;
    temp1 *= C8;
    step[ 3] = (2 * (temp1) + step1_rounding) >> step1_shift;

    temp1 = input[2] * C8;
    temp1 = (2 * (temp1) +   RIGHT_ROUNDING) >> RIGHT_SHIFT;
    temp2 = input[6] + input[10];
    step[ 4] = (temp1 + temp2 + INITIAL_ROUNDING) >> INITIAL_SHIFT;
    step[ 5] = (temp1 - temp2 + INITIAL_ROUNDING) >> INITIAL_SHIFT;

    temp1 = input[14] * C8;
    temp1 = (2 * (temp1) +   RIGHT_ROUNDING) >> RIGHT_SHIFT;
    temp2 = input[6] - input[10];
    step[ 6] = (temp2 - temp1 + INITIAL_ROUNDING) >> INITIAL_SHIFT;
    step[ 7] = (temp2 + temp1 + INITIAL_ROUNDING) >> INITIAL_SHIFT;

    // for odd input
    temp1 = input[3] * C12;
    temp2 = input[13] * C4;
    temp1 = (temp1 + temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;
    temp1 *= C8;
    intermediate[ 8] = (2 * (temp1) +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = input[3] * C4;
    temp2 = input[13] * C12;
    temp2 = (temp2 - temp1 + RIGHT_ROUNDING) >> RIGHT_SHIFT;
    temp2 *= C8;
    intermediate[ 9] = (2 * (temp2) +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    intermediate[10] = (2 * (input[9] * C8) + RIGHT_ROUNDING) >> RIGHT_SHIFT;
    intermediate[11] = input[15] - input[1];
    intermediate[12] = input[15] + input[1];
    intermediate[13] = (2 * (input[7] * C8) + RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = input[11] * C12;
    temp2 = input[5] * C4;
    temp2 = (temp2 - temp1 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;
    temp2 *= C8;
    intermediate[14] = (2 * (temp2) +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = input[11] * C4;
    temp2 = input[5] * C12;
    temp1 = (temp1 + temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;
    temp1 *= C8;
    intermediate[15] = (2 * (temp1) +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    step[ 8] = (intermediate[ 8] + intermediate[14] + INITIAL_ROUNDING)
        >> INITIAL_SHIFT;
    step[ 9] = (intermediate[ 9] + intermediate[15] + INITIAL_ROUNDING)
        >> INITIAL_SHIFT;
    step[10] = (intermediate[10] + intermediate[11] + INITIAL_ROUNDING)
        >> INITIAL_SHIFT;
    step[11] = (intermediate[10] - intermediate[11] + INITIAL_ROUNDING)
        >> INITIAL_SHIFT;
    step[12] = (intermediate[12] + intermediate[13] + INITIAL_ROUNDING)
        >> INITIAL_SHIFT;
    step[13] = (intermediate[12] - intermediate[13] + INITIAL_ROUNDING)
        >> INITIAL_SHIFT;
    step[14] = (intermediate[ 8] - intermediate[14] + INITIAL_ROUNDING)
        >> INITIAL_SHIFT;
    step[15] = (intermediate[ 9] - intermediate[15] + INITIAL_ROUNDING)
        >> INITIAL_SHIFT;

    // step 3
    output[0] = step[ 0] + step[ 3];
    output[1] = step[ 1] + step[ 2];
    output[2] = step[ 1] - step[ 2];
    output[3] = step[ 0] - step[ 3];

    temp1 = step[ 4] * C14;
    temp2 = step[ 7] * C2;
    output[4] =  (temp1 - temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[ 4] * C2;
    temp2 = step[ 7] * C14;
    output[7] =  (temp1 + temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[ 5] * C10;
    temp2 = step[ 6] * C6;
    output[5] =  (temp1 - temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[ 5] * C6;
    temp2 = step[ 6] * C10;
    output[6] =  (temp1 + temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

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

    temp1 = output[8] * C7;
    temp2 = output[15] * C9;
    step[ 8] = (temp1 - temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[9] * C11;
    temp2 = output[14] * C5;
    step[ 9] = (temp1 + temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[10] * C3;
    temp2 = output[13] * C13;
    step[10] = (temp1 - temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[11] * C15;
    temp2 = output[12] * C1;
    step[11] = (temp1 + temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[11] * C1;
    temp2 = output[12] * C15;
    step[12] = (temp2 - temp1 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[10] * C13;
    temp2 = output[13] * C3;
    step[13] = (temp1 + temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[9] * C5;
    temp2 = output[14] * C11;
    step[14] = (temp2 - temp1 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[8] * C9;
    temp2 = output[15] * C7;
    step[15] = (temp1 + temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    // step 5
    output[0] = (step[0] + step[15] + last_rounding) >> last_shift_bits;
    output[1] = (step[1] + step[14] + last_rounding) >> last_shift_bits;
    output[2] = (step[2] + step[13] + last_rounding) >> last_shift_bits;
    output[3] = (step[3] + step[12] + last_rounding) >> last_shift_bits;
    output[4] = (step[4] + step[11] + last_rounding) >> last_shift_bits;
    output[5] = (step[5] + step[10] + last_rounding) >> last_shift_bits;
    output[6] = (step[6] + step[ 9] + last_rounding) >> last_shift_bits;
    output[7] = (step[7] + step[ 8] + last_rounding) >> last_shift_bits;

    output[15] = (step[0] - step[15] + last_rounding) >> last_shift_bits;
    output[14] = (step[1] - step[14] + last_rounding) >> last_shift_bits;
    output[13] = (step[2] - step[13] + last_rounding) >> last_shift_bits;
    output[12] = (step[3] - step[12] + last_rounding) >> last_shift_bits;
    output[11] = (step[4] - step[11] + last_rounding) >> last_shift_bits;
    output[10] = (step[5] - step[10] + last_rounding) >> last_shift_bits;
    output[9] = (step[6] - step[ 9] + last_rounding) >> last_shift_bits;
    output[8] = (step[7] - step[ 8] + last_rounding) >> last_shift_bits;
}

void vp9_short_idct16x16_c(int16_t *input, int16_t *output, int pitch) {
    int16_t out[16 * 16];
    int16_t *outptr = &out[0];
    const int short_pitch = pitch >> 1;
    int i, j;
    int16_t temp_in[16], temp_out[16];

    // First transform rows
    for (i = 0; i < 16; ++i) {
      butterfly_16x16_idct_1d(input, outptr, 0);
      input += short_pitch;
      outptr += 16;
    }

    // Then transform columns
    for (i = 0; i < 16; ++i) {
      for (j = 0; j < 16; ++j)
        temp_in[j] = out[j * 16 + i];
      butterfly_16x16_idct_1d(temp_in, temp_out, 3);
      for (j = 0; j < 16; ++j)
        output[j * 16 + i] = temp_out[j];
    }
}

/* The following function is called when we know the maximum number of non-zero
 * dct coefficients is less or equal 10.
 */
static void butterfly_16x16_idct10_1d(int16_t input[16], int16_t output[16],
                                      int last_shift_bits) {
    int16_t step[16] = {0};
    int intermediate[16] = {0};
    int temp1, temp2;
    int last_rounding = 0;

    if (last_shift_bits > 0)
      last_rounding = 1 << (last_shift_bits - 1);

    // step 1 and 2
    step[ 0] = (input[0] + INITIAL_ROUNDING) >> INITIAL_SHIFT;
    step[ 1] = (input[0] + INITIAL_ROUNDING) >> INITIAL_SHIFT;

    temp1 = (2 * (input[2] * C8) + RIGHT_ROUNDING) >> RIGHT_SHIFT;
    step[ 4] = (temp1 + INITIAL_ROUNDING) >> INITIAL_SHIFT;
    step[ 5] = (temp1 + INITIAL_ROUNDING) >> INITIAL_SHIFT;

    // for odd input
    temp1 = (input[3] * C12 + RIGHT_ROUNDING) >> RIGHT_SHIFT;
    temp1 *= C8;
    intermediate[ 8] = (2 * (temp1) + RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = (-input[3] * C4 + RIGHT_ROUNDING) >> RIGHT_SHIFT;
    temp1 *= C8;
    intermediate[ 9] = (2 * (temp1) + RIGHT_ROUNDING) >> RIGHT_SHIFT;

    step[ 8] = (intermediate[ 8] + INITIAL_ROUNDING) >> INITIAL_SHIFT;
    step[ 9] = (intermediate[ 9] + INITIAL_ROUNDING) >> INITIAL_SHIFT;
    step[10] = (-input[1] + INITIAL_ROUNDING) >> INITIAL_SHIFT;
    step[11] = (input[1] + INITIAL_ROUNDING) >> INITIAL_SHIFT;
    step[12] = (input[1] + INITIAL_ROUNDING) >> INITIAL_SHIFT;
    step[13] = (input[1] + INITIAL_ROUNDING) >> INITIAL_SHIFT;
    step[14] = (intermediate[ 8] + INITIAL_ROUNDING) >> INITIAL_SHIFT;
    step[15] = (intermediate[ 9] + INITIAL_ROUNDING) >> INITIAL_SHIFT;

    // step 3
    output[0] = step[ 0];
    output[1] = step[ 1];
    output[2] = step[ 1];
    output[3] = step[ 0];

    temp1 = step[ 4] * C14;
    output[4] =  (temp1 + RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[ 4] * C2;
    output[7] =  (temp1 + RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[ 5] * C10;
    output[5] =  (temp1 + RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[ 5] * C6;
    output[6] =  (temp1 + RIGHT_ROUNDING) >> RIGHT_SHIFT;

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

    temp1 = output[8] * C7;
    temp2 = output[15] * C9;
    step[ 8] = (temp1 - temp2 + RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[9] * C11;
    temp2 = output[14] * C5;
    step[ 9] = (temp1 + temp2 + RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[10] * C3;
    temp2 = output[13] * C13;
    step[10] = (temp1 - temp2 + RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[11] * C15;
    temp2 = output[12] * C1;
    step[11] = (temp1 + temp2 + RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[11] * C1;
    temp2 = output[12] * C15;
    step[12] = (temp2 - temp1 + RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[10] * C13;
    temp2 = output[13] * C3;
    step[13] = (temp1 + temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[9] * C5;
    temp2 = output[14] * C11;
    step[14] = (temp2 - temp1 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[8] * C9;
    temp2 = output[15] * C7;
    step[15] = (temp1 + temp2 +   RIGHT_ROUNDING) >> RIGHT_SHIFT;

    // step 5
    output[0] = (step[0] + step[15] + last_rounding) >> last_shift_bits;
    output[1] = (step[1] + step[14] + last_rounding) >> last_shift_bits;
    output[2] = (step[2] + step[13] + last_rounding) >> last_shift_bits;
    output[3] = (step[3] + step[12] + last_rounding) >> last_shift_bits;
    output[4] = (step[4] + step[11] + last_rounding) >> last_shift_bits;
    output[5] = (step[5] + step[10] + last_rounding) >> last_shift_bits;
    output[6] = (step[6] + step[ 9] + last_rounding) >> last_shift_bits;
    output[7] = (step[7] + step[ 8] + last_rounding) >> last_shift_bits;

    output[15] = (step[0] - step[15] + last_rounding) >> last_shift_bits;
    output[14] = (step[1] - step[14] + last_rounding) >> last_shift_bits;
    output[13] = (step[2] - step[13] + last_rounding) >> last_shift_bits;
    output[12] = (step[3] - step[12] + last_rounding) >> last_shift_bits;
    output[11] = (step[4] - step[11] + last_rounding) >> last_shift_bits;
    output[10] = (step[5] - step[10] + last_rounding) >> last_shift_bits;
    output[9] = (step[6] - step[ 9] + last_rounding) >> last_shift_bits;
    output[8] = (step[7] - step[ 8] + last_rounding) >> last_shift_bits;
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
      butterfly_16x16_idct10_1d(input, outptr, 0);
      input += short_pitch;
      outptr += 16;
    }

    // Then transform columns
    for (i = 0; i < 16; ++i) {
      for (j = 0; j < 16; ++j)
        temp_in[j] = out[j*16 + i];
      butterfly_16x16_idct10_1d(temp_in, temp_out, 3);
      for (j = 0; j < 16; ++j)
        output[j*16 + i] = temp_out[j];
    }
}
#undef INITIAL_SHIFT
#undef INITIAL_ROUNDING
#undef RIGHT_SHIFT
#undef RIGHT_ROUNDING
#endif
