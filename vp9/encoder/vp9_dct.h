/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_DCT_H_
#define VP9_ENCODER_VP9_DCT_H_

#include "vp9/common/vp9_idct.h"

#ifdef __cplusplus
extern "C" {
#endif

void vp9_highbd_fdct4x4_c(const int16_t *input, tran_low_t *output, int stride);
void vp9_highbd_fdct8x8_c(const int16_t *input, tran_low_t *output, int stride);
void vp9_highbd_fdct16x16_c(const int16_t *input, tran_low_t *output,
                            int stride);
void vp9_highbd_fdct32x32_c(const int16_t *input, tran_low_t *out, int stride);
void vp9_highbd_fdct32x32_rd_c(const int16_t *input, tran_low_t *out,
                               int stride);

void vp9_fdct4(const tran_low_t *input, tran_low_t *output);
void vp9_fadst4(const tran_low_t *input, tran_low_t *output);
void vp9_fdct8(const tran_low_t *input, tran_low_t *output);
void vp9_fadst8(const tran_low_t *input, tran_low_t *output);
void vp9_fdct16(const tran_low_t in[16], tran_low_t out[16]);
void vp9_fadst16(const tran_low_t *input, tran_low_t *output);
void vp9_fdct32(const tran_high_t *input, tran_high_t *output, int round);

#if CONFIG_EXT_TX
void vp9_fdst4(const tran_low_t *input, tran_low_t *output);
void vp9_fdst8(const tran_low_t *input, tran_low_t *output);
void vp9_fdst16(const tran_low_t *input, tran_low_t *output);
#endif  // CONFIG_EXT_TX

static const transform_2d FHT_4[] = {
  { vp9_fdct4,  vp9_fdct4  },  // DCT_DCT  = 0
  { vp9_fadst4, vp9_fdct4  },  // ADST_DCT = 1
  { vp9_fdct4,  vp9_fadst4 },  // DCT_ADST = 2
  { vp9_fadst4, vp9_fadst4 },  // ADST_ADST = 3
#if CONFIG_EXT_TX
  { vp9_fadst4, vp9_fdct4  },  // FLIPADST_DCT = 4
  { vp9_fdct4,  vp9_fadst4 },  // DCT_FLIPADST = 5
  { vp9_fadst4, vp9_fadst4 },  // FLIPADST_FLIPADST = 6
  { vp9_fadst4, vp9_fadst4 },  // ADST_FLIPADST = 7
  { vp9_fadst4, vp9_fadst4 },  // FLIPADST_ADST = 8
  { vp9_fdst4,  vp9_fdst4  },  // DST_DST = 9
  { vp9_fdst4,  vp9_fdct4  },  // DST_DCT = 10
  { vp9_fdct4,  vp9_fdst4  },  // DCT_DST = 11
  { vp9_fdst4,  vp9_fadst4 },  // DST_ADST = 12
  { vp9_fadst4, vp9_fdst4  },  // ADST_DST = 13
  { vp9_fdst4,  vp9_fadst4 },  // DST_FLIPADST = 14
  { vp9_fadst4, vp9_fdst4  },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
};

static const transform_2d FHT_8[] = {
  { vp9_fdct8,  vp9_fdct8  },  // DCT_DCT  = 0
  { vp9_fadst8, vp9_fdct8  },  // ADST_DCT = 1
  { vp9_fdct8,  vp9_fadst8 },  // DCT_ADST = 2
  { vp9_fadst8, vp9_fadst8 },  // ADST_ADST = 3
#if CONFIG_EXT_TX
  { vp9_fadst8, vp9_fdct8  },  // FLIPADST_DCT = 4
  { vp9_fdct8,  vp9_fadst8 },  // DCT_FLIPADST = 5
  { vp9_fadst8, vp9_fadst8 },  // FLIPADST_FLIPADST = 6
  { vp9_fadst8, vp9_fadst8 },  // ADST_FLIPADST = 7
  { vp9_fadst8, vp9_fadst8 },  // FLIPADST_ADST = 8
  { vp9_fdst8,  vp9_fdst8  },  // DST_DST = 9
  { vp9_fdst8,  vp9_fdct8  },  // DST_DCT = 10
  { vp9_fdct8,  vp9_fdst8  },  // DCT_DST = 11
  { vp9_fdst8,  vp9_fadst8 },  // DST_ADST = 12
  { vp9_fadst8, vp9_fdst8  },  // ADST_DST = 13
  { vp9_fdst8,  vp9_fadst8 },  // DST_FLIPADST = 14
  { vp9_fadst8, vp9_fdst8  },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
};

static const transform_2d FHT_16[] = {
  { vp9_fdct16,  vp9_fdct16  },  // DCT_DCT  = 0
  { vp9_fadst16, vp9_fdct16  },  // ADST_DCT = 1
  { vp9_fdct16,  vp9_fadst16 },  // DCT_ADST = 2
  { vp9_fadst16, vp9_fadst16 },  // ADST_ADST = 3
#if CONFIG_EXT_TX
  { vp9_fadst16, vp9_fdct16  },  // FLIPADST_DCT = 4
  { vp9_fdct16,  vp9_fadst16 },  // DCT_FLIPADST = 5
  { vp9_fadst16, vp9_fadst16 },  // FLIPADST_FLIPADST = 6
  { vp9_fadst16, vp9_fadst16 },  // ADST_FLIPADST = 7
  { vp9_fadst16, vp9_fadst16 },  // FLIPADST_ADST = 8
  { vp9_fdst16,  vp9_fdst16  },  // DST_DST = 9
  { vp9_fdst16,  vp9_fdct16  },  // DST_DCT = 10
  { vp9_fdct16,  vp9_fdst16  },  // DCT_DST = 11
  { vp9_fdst16,  vp9_fadst16 },  // DST_ADST = 12
  { vp9_fadst16, vp9_fdst16  },  // ADST_DST = 13
  { vp9_fdst16,  vp9_fadst16 },  // DST_FLIPADST = 14
  { vp9_fadst16, vp9_fdst16  },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
};

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_DCT_H_
