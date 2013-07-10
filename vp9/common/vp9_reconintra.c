/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <stdio.h>

#include "./vpx_config.h"
#include "vp9_rtcd.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_onyxc_int.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/vpx_once.h"

const TX_TYPE mode2txfm_map[MB_MODE_COUNT] = {
    DCT_DCT,    // DC
    ADST_DCT,   // V
    DCT_ADST,   // H
    DCT_DCT,    // D45
    ADST_ADST,  // D135
    ADST_DCT,   // D117
    DCT_ADST,   // D153
    DCT_ADST,   // D27
    ADST_DCT,   // D63
    ADST_ADST,  // TM
    DCT_DCT,    // NEARESTMV
    DCT_DCT,    // NEARMV
    DCT_DCT,    // ZEROMV
    DCT_DCT     // NEWMV
};

#define intra_pred_sized(type, size) \
void vp9_##type##_predictor_##size##x##size##_c(uint8_t *pred_ptr, \
                                                ptrdiff_t stride, \
                                                uint8_t *above_row, \
                                                uint8_t *left_col) { \
  type##_predictor(pred_ptr, stride, size, above_row, left_col); \
}
#define intra_pred_allsizes(type) \
  intra_pred_sized(type, 4) \
  intra_pred_sized(type, 8) \
  intra_pred_sized(type, 16) \
  intra_pred_sized(type, 32)

static INLINE void d27_predictor(uint8_t *pred_ptr, ptrdiff_t stride, int bs,
                                 uint8_t *above_row, uint8_t *left_col) {
  int r, c;
  // first column
  for (r = 0; r < bs - 1; ++r) {
      pred_ptr[r * stride] = ROUND_POWER_OF_TWO(left_col[r] +
                                                   left_col[r + 1], 1);
  }
  pred_ptr[(bs - 1) * stride] = left_col[bs - 1];
  pred_ptr++;
  // second column
  for (r = 0; r < bs - 2; ++r) {
      pred_ptr[r * stride] = ROUND_POWER_OF_TWO(left_col[r] +
                                                   left_col[r + 1] * 2 +
                                                   left_col[r + 2], 2);
  }
  pred_ptr[(bs - 2) * stride] = ROUND_POWER_OF_TWO(left_col[bs - 2] +
                                                      left_col[bs - 1] * 3,
                                                      2);
  pred_ptr[(bs - 1) * stride] = left_col[bs - 1];
  pred_ptr++;

  // rest of last row
  for (c = 0; c < bs - 2; ++c) {
    pred_ptr[(bs - 1) * stride + c] = left_col[bs - 1];
  }

  for (r = bs - 2; r >= 0; --r) {
    for (c = 0; c < bs - 2; ++c) {
      pred_ptr[r * stride + c] = pred_ptr[(r + 1) * stride + c - 2];
    }
  }
}
intra_pred_allsizes(d27)

static INLINE void d63_predictor(uint8_t *pred_ptr, ptrdiff_t stride, int bs,
                                 uint8_t *above_row, uint8_t *left_col) {
  int r, c;
  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      if (r & 1) {
        pred_ptr[c] = ROUND_POWER_OF_TWO(above_row[r/2 + c] +
                                         above_row[r/2 + c + 1] * 2 +
                                         above_row[r/2 + c + 2], 2);
      } else {
        pred_ptr[c] = ROUND_POWER_OF_TWO(above_row[r/2 + c] +
                                         above_row[r/2+ c + 1], 1);
      }
    }
    pred_ptr += stride;
  }
}
intra_pred_allsizes(d63)

static INLINE void d45_predictor(uint8_t *pred_ptr, ptrdiff_t stride, int bs,
                                 uint8_t *above_row, uint8_t *left_col) {
  int r, c;
  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      if (r + c + 2 < bs * 2)
        pred_ptr[c] = ROUND_POWER_OF_TWO(above_row[r + c] +
                                         above_row[r + c + 1] * 2 +
                                         above_row[r + c + 2], 2);
      else
        pred_ptr[c] = above_row[bs * 2 - 1];
    }
    pred_ptr += stride;
  }
}
intra_pred_allsizes(d45)

static INLINE void d117_predictor(uint8_t *pred_ptr, ptrdiff_t stride, int bs,
                                  uint8_t *above_row, uint8_t *left_col) {
  int r, c;
  // first row
  for (c = 0; c < bs; c++)
    pred_ptr[c] = ROUND_POWER_OF_TWO(above_row[c - 1] + above_row[c], 1);
  pred_ptr += stride;

  // second row
  pred_ptr[0] = ROUND_POWER_OF_TWO(left_col[0] +
                                   above_row[-1] * 2 +
                                   above_row[0], 2);
  for (c = 1; c < bs; c++)
    pred_ptr[c] = ROUND_POWER_OF_TWO(above_row[c - 2] +
                                     above_row[c - 1] * 2 +
                                     above_row[c], 2);
  pred_ptr += stride;

  // the rest of first col
  pred_ptr[0] = ROUND_POWER_OF_TWO(above_row[-1] +
                                   left_col[0] * 2 +
                                   left_col[1], 2);
  for (r = 3; r < bs; ++r)
    pred_ptr[(r-2) * stride] = ROUND_POWER_OF_TWO(left_col[r - 3] +
                                                  left_col[r - 2] * 2 +
                                                  left_col[r - 1], 2);
  // the rest of the block
  for (r = 2; r < bs; ++r) {
    for (c = 1; c < bs; c++)
      pred_ptr[c] = pred_ptr[-2 * stride + c - 1];
    pred_ptr += stride;
  }
}
intra_pred_allsizes(d117)

static INLINE void d135_predictor(uint8_t *pred_ptr, ptrdiff_t stride, int bs,
                                  uint8_t *above_row, uint8_t *left_col) {
  int r, c;
  pred_ptr[0] = ROUND_POWER_OF_TWO(left_col[0] +
                                   above_row[-1] * 2 +
                                   above_row[0], 2);
  for (c = 1; c < bs; c++)
    pred_ptr[c] = ROUND_POWER_OF_TWO(above_row[c - 2] +
                                     above_row[c - 1] * 2 +
                                     above_row[c], 2);

  pred_ptr[stride] = ROUND_POWER_OF_TWO(above_row[-1] +
                                        left_col[0] * 2 +
                                        left_col[1], 2);
  for (r = 2; r < bs; ++r)
    pred_ptr[r * stride] = ROUND_POWER_OF_TWO(left_col[r - 2] +
                                              left_col[r - 1] * 2 +
                                              left_col[r], 2);

  pred_ptr += stride;
  for (r = 1; r < bs; ++r) {
    for (c = 1; c < bs; c++)
      pred_ptr[c] = pred_ptr[-stride + c - 1];
    pred_ptr += stride;
  }
}
intra_pred_allsizes(d135)

static INLINE void d153_predictor(uint8_t *pred_ptr, ptrdiff_t stride, int bs,
                                  uint8_t *above_row, uint8_t *left_col) {
  int r, c;
  pred_ptr[0] = ROUND_POWER_OF_TWO(above_row[-1] + left_col[0], 1);
  for (r = 1; r < bs; r++)
    pred_ptr[r * stride] =
        ROUND_POWER_OF_TWO(left_col[r - 1] + left_col[r], 1);
  pred_ptr++;

  pred_ptr[0] = ROUND_POWER_OF_TWO(left_col[0] +
                                   above_row[-1] * 2 +
                                   above_row[0], 2);
  pred_ptr[stride] = ROUND_POWER_OF_TWO(above_row[-1] +
                                        left_col[0] * 2 +
                                        left_col[1], 2);
  for (r = 2; r < bs; r++)
    pred_ptr[r * stride] = ROUND_POWER_OF_TWO(left_col[r - 2] +
                                              left_col[r - 1] * 2 +
                                              left_col[r], 2);
  pred_ptr++;

  for (c = 0; c < bs - 2; c++)
    pred_ptr[c] = ROUND_POWER_OF_TWO(above_row[c - 1] +
                                     above_row[c] * 2 +
                                     above_row[c + 1], 2);
  pred_ptr += stride;
  for (r = 1; r < bs; ++r) {
    for (c = 0; c < bs - 2; c++)
      pred_ptr[c] = pred_ptr[-stride + c - 2];
    pred_ptr += stride;
  }
}
intra_pred_allsizes(d153)

static INLINE void v_predictor(uint8_t *pred_ptr, ptrdiff_t stride, int bs,
                       uint8_t *above_row, uint8_t *left_col) {
  int r;

  for (r = 0; r < bs; r++) {
    vpx_memcpy(pred_ptr, above_row, bs);
    pred_ptr += stride;
  }
}
intra_pred_allsizes(v)

static INLINE void h_predictor(uint8_t *pred_ptr, ptrdiff_t stride, int bs,
                               uint8_t *above_row, uint8_t *left_col) {
  int r;

  for (r = 0; r < bs; r++) {
    vpx_memset(pred_ptr, left_col[r], bs);
    pred_ptr += stride;
  }
}
intra_pred_allsizes(h)

static INLINE void tm_predictor(uint8_t *pred_ptr, ptrdiff_t stride, int bs,
                                uint8_t *above_row, uint8_t *left_col) {
  int r, c;
  int ytop_left = above_row[-1];

  for (r = 0; r < bs; r++) {
    for (c = 0; c < bs; c++)
      pred_ptr[c] = clip_pixel(left_col[r] + above_row[c] - ytop_left);
    pred_ptr += stride;
  }
}
intra_pred_allsizes(tm)

static INLINE void dc_128_predictor(uint8_t *pred_ptr, ptrdiff_t stride, int bs,
                                    uint8_t *above_row, uint8_t *left_col) {
  int r;

  for (r = 0; r < bs; r++) {
    vpx_memset(pred_ptr, 128, bs);
    pred_ptr += stride;
  }
}
intra_pred_allsizes(dc_128)

static INLINE void dc_left_predictor(uint8_t *pred_ptr, ptrdiff_t stride,
                                     int bs,
                                     uint8_t *above_row, uint8_t *left_col) {
  int i, r;
  int expected_dc = 128;
  int average = 0;
  const int count = bs;

  for (i = 0; i < bs; i++)
    average += left_col[i];
  expected_dc = (average + (count >> 1)) / count;

  for (r = 0; r < bs; r++) {
    vpx_memset(pred_ptr, expected_dc, bs);
    pred_ptr += stride;
  }
}
intra_pred_allsizes(dc_left)

static INLINE void dc_top_predictor(uint8_t *pred_ptr, ptrdiff_t stride, int bs,
                                    uint8_t *above_row, uint8_t *left_col) {
  int i, r;
  int expected_dc = 128;
  int average = 0;
  const int count = bs;

  for (i = 0; i < bs; i++)
    average += above_row[i];
  expected_dc = (average + (count >> 1)) / count;

  for (r = 0; r < bs; r++) {
    vpx_memset(pred_ptr, expected_dc, bs);
    pred_ptr += stride;
  }
}
intra_pred_allsizes(dc_top)

static INLINE void dc_predictor(uint8_t *pred_ptr, ptrdiff_t stride, int bs,
                                uint8_t *above_row, uint8_t *left_col) {
  int i, r;
  int expected_dc = 128;
  int average = 0;
  const int count = 2 * bs;

  for (i = 0; i < bs; i++)
    average += above_row[i];
  for (i = 0; i < bs; i++)
    average += left_col[i];
  expected_dc = (average + (count >> 1)) / count;

  for (r = 0; r < bs; r++) {
    vpx_memset(pred_ptr, expected_dc, bs);
    pred_ptr += stride;
  }
}
intra_pred_allsizes(dc)
#undef intra_pred_allsizes

typedef void (*intra_pred_fn)(uint8_t *pred_ptr, ptrdiff_t stride,
                              uint8_t *above_row, uint8_t *left_col);

static intra_pred_fn pred[VP9_INTRA_MODES][4];
static intra_pred_fn dc_pred[2][2][4];

static void init_intra_pred_fn_ptrs(void) {
#define intra_pred_allsizes(l, type) \
  l[0] = vp9_##type##_predictor_4x4; \
  l[1] = vp9_##type##_predictor_8x8; \
  l[2] = vp9_##type##_predictor_16x16; \
  l[3] = vp9_##type##_predictor_32x32

  intra_pred_allsizes(pred[V_PRED], v);
  intra_pred_allsizes(pred[H_PRED], h);
  intra_pred_allsizes(pred[D27_PRED], d27);
  intra_pred_allsizes(pred[D45_PRED], d45);
  intra_pred_allsizes(pred[D63_PRED], d63);
  intra_pred_allsizes(pred[D117_PRED], d117);
  intra_pred_allsizes(pred[D135_PRED], d135);
  intra_pred_allsizes(pred[D153_PRED], d153);
  intra_pred_allsizes(pred[TM_PRED], tm);

  intra_pred_allsizes(dc_pred[0][0], dc_128);
  intra_pred_allsizes(dc_pred[0][1], dc_top);
  intra_pred_allsizes(dc_pred[1][0], dc_left);
  intra_pred_allsizes(dc_pred[1][1], dc);

#undef intra_pred_allsizes
}

static void build_intra_predictors(uint8_t *src, int src_stride,
                                   uint8_t *pred_ptr, int stride,
                                   MB_PREDICTION_MODE mode, TX_SIZE txsz,
                                   int up_available, int left_available,
                                   int right_available) {
  int i;
  DECLARE_ALIGNED_ARRAY(16, uint8_t, left_col, 64);
  DECLARE_ALIGNED_ARRAY(16, uint8_t, yabove_data, 128 + 16);
  uint8_t *above_row = yabove_data + 16;
  const int bs = 4 << txsz;

  // 127 127 127 .. 127 127 127 127 127 127
  // 129  A   B  ..  Y   Z
  // 129  C   D  ..  W   X
  // 129  E   F  ..  U   V
  // 129  G   H  ..  S   T   T   T   T   T
  // ..

  once(init_intra_pred_fn_ptrs);
  if (left_available) {
    for (i = 0; i < bs; i++)
      left_col[i] = src[i * src_stride - 1];
  } else {
    vpx_memset(left_col, 129, bs);
  }

  if (up_available) {
    uint8_t *above_ptr = src - src_stride;
    if (bs == 4 && right_available && left_available) {
      above_row = above_ptr;
    } else {
      vpx_memcpy(above_row, above_ptr, bs);
      if (bs == 4 && right_available)
        vpx_memcpy(above_row + bs, above_ptr + bs, bs);
      else
        vpx_memset(above_row + bs, above_row[bs - 1], bs);
      above_row[-1] = left_available ? above_ptr[-1] : 129;
    }
  } else {
    vpx_memset(above_row, 127, bs * 2);
    above_row[-1] = 127;
  }

  if (mode == DC_PRED) {
    dc_pred[left_available][up_available][txsz](pred_ptr, stride,
                                                above_row, left_col);
  } else {
    pred[mode][txsz](pred_ptr, stride, above_row, left_col);
  }
}

void vp9_predict_intra_block(MACROBLOCKD *xd,
                            int block_idx,
                            int bwl_in,
                            TX_SIZE tx_size,
                            int mode,
                            uint8_t *reference, int ref_stride,
                            uint8_t *predictor, int pre_stride) {
  const int bwl = bwl_in - tx_size;
  const int wmask = (1 << bwl) - 1;
  const int have_top = (block_idx >> bwl) || xd->up_available;
  const int have_left = (block_idx & wmask) || xd->left_available;
  const int have_right = ((block_idx & wmask) != wmask);

  assert(bwl >= 0);
  build_intra_predictors(reference, ref_stride,
                         predictor, pre_stride,
                         mode,
                         tx_size,
                         have_top, have_left,
                         have_right);
}
