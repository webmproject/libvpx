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

#if CONFIG_FILTERINTRA
static void filter_intra_predictors(uint8_t *ypred_ptr, int y_stride, int bs,
                                    uint8_t *yabove_row, uint8_t *yleft_col,
                                    int mode) {
  static const int prec_bits = 10;
  static const int round_val = 511;
  static const int taps[10][3] = {
      {438 , 660, -352},  // DC
      {1014, 565, -559},  // V
      {312, 1017, -312},  // H
      {0, 0, 0},          // D45
      {478, 483, 153},    // D135
      {699, 470, -122},   // D117
      {356, 707, 35},     // D153
      {0, 0, 0},          // D27
      {0, 0, 0},          // D63
      {877, 896, -812}    // TM
      };
  int k, r, c;
  int pred[17][17];
  int mean, ipred;
  const int c1 = taps[mode][0];
  const int c2 = taps[mode][1];
  const int c3 = taps[mode][2];

  k = 0;
  mean = 0;
  while (k < bs) {
    mean = mean + (int)yleft_col[r];
    mean = mean + (int)yabove_row[c];
    k++;
  }
  mean = (mean + bs) / (2 * bs);

  for (r = 0; r < bs; r++)
    pred[r + 1][0] = (int)yleft_col[r] - mean;

  for (c = 0; c < bs + 1; c++)
    pred[0][c] = (int)yabove_row[c - 1] - mean;

  for (r = 1; r < bs + 1; r++)
    for (c = 1; c < bs + 1; c++) {
      ipred = c1 * pred[r - 1][c] + c2 * pred[r][c - 1]
              + c3 * pred[r - 1][c - 1];
      pred[r][c] = ipred < 0 ? -((-ipred + round_val) >> prec_bits) :
                               ((ipred + round_val) >> prec_bits);
    }

  for (r = 0; r < bs; r++) {
    for (c = 0; c < bs; c++) {
      ipred = pred[r + 1][c + 1] + mean;
      ypred_ptr[c] = clip_pixel(ipred);
    }
    ypred_ptr += y_stride;
  }
}

static void build_filter_intra_predictors(uint8_t *src, int src_stride,
                                          uint8_t *pred_ptr, int stride,
                                          MB_PREDICTION_MODE mode, TX_SIZE txsz,
                                          int up_available, int left_available,
                                          int right_available) {
  int i;
  DECLARE_ALIGNED_ARRAY(16, uint8_t, left_col, 64);
  DECLARE_ALIGNED_ARRAY(16, uint8_t, yabove_data, 128 + 16);
  uint8_t *above_row = yabove_data + 16;
  const int bs = 4 << txsz;

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

  filter_intra_predictors(pred_ptr, stride, bs, above_row, left_col, mode);
}
#endif

void vp9_predict_intra_block(MACROBLOCKD *xd,
                            int block_idx,
                            int bwl_in,
                            TX_SIZE tx_size,
                            int mode,
#if CONFIG_FILTERINTRA
                            int filterbit,
#endif
                            uint8_t *reference, int ref_stride,
                            uint8_t *predictor, int pre_stride) {
  const int bwl = bwl_in - tx_size;
  const int wmask = (1 << bwl) - 1;
  const int have_top = (block_idx >> bwl) || xd->up_available;
  const int have_left = (block_idx & wmask) || xd->left_available;
  const int have_right = ((block_idx & wmask) != wmask);
#if CONFIG_FILTERINTRA
  int filterflag = is_filter_allowed(mode) && (tx_size <= TX_8X8) && filterbit;
#endif

  assert(bwl >= 0);
#if CONFIG_FILTERINTRA
  if (!filterflag) {
#endif
  build_intra_predictors(reference, ref_stride,
                         predictor, pre_stride,
                         mode,
                         tx_size,
                         have_top, have_left,
                         have_right);
#if CONFIG_FILTERINTRA
  } else {
    build_filter_intra_predictors(reference, ref_stride,
                                  predictor, pre_stride,
                                  mode,
                                  tx_size,
                                  have_top, have_left,
                                  have_right);
  }
#endif
}

#if CONFIG_INTERINTRA
// Intra predictor for the second square block in interintra prediction.
// Prediction of the first block (in pred_ptr) will be used to generate half of
// the boundary values.
static void build_intra_predictors_for_2nd_block_interintra
                                   (uint8_t *src, int src_stride,
                                   uint8_t *pred_ptr, int stride,
                                   MB_PREDICTION_MODE mode, TX_SIZE txsz,
                                   int up_available, int left_available,
                                   int right_available, int bwltbh) {
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
    for (i = 0; i < bs; i++) {
      if (bwltbh)
        left_col[i] = src[i * src_stride - 1];
      else
        left_col[i] = pred_ptr[i * stride - 1];
    }
  } else {
    vpx_memset(left_col, 129, bs);
  }

  if (up_available) {
    uint8_t *above_ptr;
    if (bwltbh)
      above_ptr = pred_ptr - stride;
    else
      above_ptr = src - src_stride;
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

#if CONFIG_MASKED_COMPOUND
#define MASK_WEIGHT_BITS_INTERINTRA 6

static int get_masked_weight_interintra(int m) {
  #define SMOOTHER_LEN_INTERINTRA  32
  static const uint8_t smoothfn[2 * SMOOTHER_LEN_INTERINTRA + 1] = {
      0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  1,  1,  1,
      1,  1,  2,  2,  3,  4,  5,  6,
      8,  9, 12, 14, 17, 21, 24, 28,
      32,
      36, 40, 43, 47, 50, 52, 55, 56,
      58, 59, 60, 61, 62, 62, 63, 63,
      63, 63, 63, 64, 64, 64, 64, 64,
      64, 64, 64, 64, 64, 64, 64, 64,
  };
  if (m < -SMOOTHER_LEN_INTERINTRA)
    return 0;
  else if (m > SMOOTHER_LEN_INTERINTRA)
    return (1 << MASK_WEIGHT_BITS_INTERINTRA);
  else
    return smoothfn[m + SMOOTHER_LEN_INTERINTRA];
}

static int get_hard_mask_interintra(int m) {
  return m > 0;
}

// Equation of line: f(x, y) = a[0]*(x - a[2]*w/4) + a[1]*(y - a[3]*h/4) = 0
// The soft mask is obtained by computing f(x, y) and then calling
// get_masked_weight(f(x, y)).
static const int mask_params_sml_interintra[1 << MASK_BITS_SML_INTERINTRA]
                                            [4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},
};

static const int mask_params_med_hgtw_interintra[1 << MASK_BITS_MED_INTERINTRA]
                                                 [4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},

  {-1,  2, 2, 1},
  { 1, -2, 2, 1},
  {-1,  2, 2, 3},
  { 1, -2, 2, 3},
  { 1,  2, 2, 1},
  {-1, -2, 2, 1},
  { 1,  2, 2, 3},
  {-1, -2, 2, 3},
};

static const int mask_params_med_hltw_interintra[1 << MASK_BITS_MED_INTERINTRA]
                                                 [4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},

  {-2,  1, 1, 2},
  { 2, -1, 1, 2},
  {-2,  1, 3, 2},
  { 2, -1, 3, 2},
  { 2,  1, 1, 2},
  {-2, -1, 1, 2},
  { 2,  1, 3, 2},
  {-2, -1, 3, 2},
};

static const int mask_params_med_heqw_interintra[1 << MASK_BITS_MED_INTERINTRA]
                                                 [4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},

  { 0,  2, 0, 1},
  { 0, -2, 0, 1},
  { 0,  2, 0, 3},
  { 0, -2, 0, 3},
  { 2,  0, 1, 0},
  {-2,  0, 1, 0},
  { 2,  0, 3, 0},
  {-2,  0, 3, 0},
};

static const int mask_params_big_hgtw_interintra[1 << MASK_BITS_BIG_INTERINTRA]
                                                 [4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},

  {-1,  2, 2, 1},
  { 1, -2, 2, 1},
  {-1,  2, 2, 3},
  { 1, -2, 2, 3},
  { 1,  2, 2, 1},
  {-1, -2, 2, 1},
  { 1,  2, 2, 3},
  {-1, -2, 2, 3},

  {-2,  1, 1, 2},
  { 2, -1, 1, 2},
  {-2,  1, 3, 2},
  { 2, -1, 3, 2},
  { 2,  1, 1, 2},
  {-2, -1, 1, 2},
  { 2,  1, 3, 2},
  {-2, -1, 3, 2},

  { 0,  2, 0, 1},
  { 0, -2, 0, 1},
  { 0,  2, 0, 2},
  { 0, -2, 0, 2},
  { 0,  2, 0, 3},
  { 0, -2, 0, 3},
  { 2,  0, 2, 0},
  {-2,  0, 2, 0},
};

static const int mask_params_big_hltw_interintra[1 << MASK_BITS_BIG_INTERINTRA]
                                                 [4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},

  {-1,  2, 2, 1},
  { 1, -2, 2, 1},
  {-1,  2, 2, 3},
  { 1, -2, 2, 3},
  { 1,  2, 2, 1},
  {-1, -2, 2, 1},
  { 1,  2, 2, 3},
  {-1, -2, 2, 3},

  {-2,  1, 1, 2},
  { 2, -1, 1, 2},
  {-2,  1, 3, 2},
  { 2, -1, 3, 2},
  { 2,  1, 1, 2},
  {-2, -1, 1, 2},
  { 2,  1, 3, 2},
  {-2, -1, 3, 2},

  { 0,  2, 0, 2},
  { 0, -2, 0, 2},
  { 2,  0, 1, 0},
  {-2,  0, 1, 0},
  { 2,  0, 2, 0},
  {-2,  0, 2, 0},
  { 2,  0, 3, 0},
  {-2,  0, 3, 0},
};

static const int mask_params_big_heqw_interintra[1 << MASK_BITS_BIG_INTERINTRA]
                                                 [4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},

  {-1,  2, 2, 1},
  { 1, -2, 2, 1},
  {-1,  2, 2, 3},
  { 1, -2, 2, 3},
  { 1,  2, 2, 1},
  {-1, -2, 2, 1},
  { 1,  2, 2, 3},
  {-1, -2, 2, 3},

  {-2,  1, 1, 2},
  { 2, -1, 1, 2},
  {-2,  1, 3, 2},
  { 2, -1, 3, 2},
  { 2,  1, 1, 2},
  {-2, -1, 1, 2},
  { 2,  1, 3, 2},
  {-2, -1, 3, 2},

  { 0,  2, 0, 1},
  { 0, -2, 0, 1},
  { 0,  2, 0, 3},
  { 0, -2, 0, 3},
  { 2,  0, 1, 0},
  {-2,  0, 1, 0},
  { 2,  0, 3, 0},
  {-2,  0, 3, 0},
};

static const int *get_mask_params_interintra(int mask_index,
                                             BLOCK_SIZE_TYPE sb_type,
                                             int h, int w) {
  const int *a;
  const int mask_bits = get_mask_bits_interintra(sb_type);

  if (mask_index == MASK_NONE_INTERINTRA)
    return NULL;

  if (mask_bits == MASK_BITS_SML_INTERINTRA) {
    a = mask_params_sml_interintra[mask_index];
  } else if (mask_bits == MASK_BITS_MED_INTERINTRA) {
    if (h > w)
      a = mask_params_med_hgtw_interintra[mask_index];
    else if (h < w)
      a = mask_params_med_hltw_interintra[mask_index];
    else
      a = mask_params_med_heqw_interintra[mask_index];
  } else if (mask_bits == MASK_BITS_BIG_INTERINTRA) {
    if (h > w)
      a = mask_params_big_hgtw_interintra[mask_index];
    else if (h < w)
      a = mask_params_big_hltw_interintra[mask_index];
    else
      a = mask_params_big_heqw_interintra[mask_index];
  } else {
    assert(0);
  }
  return a;
}

void vp9_generate_masked_weight_interintra(int mask_index,
                                           BLOCK_SIZE_TYPE sb_type,
                                           int h, int w,
                                           uint8_t *mask, int stride) {
  int i, j;
  const int *a = get_mask_params_interintra(mask_index, sb_type, h, w);
  if (!a) return;
  for (i = 0; i < h; ++i)
    for (j = 0; j < w; ++j) {
      int x = (j - (a[2] * w) / 4);
      int y = (i - (a[3] * h) / 4);
      int m = a[0] * x + a[1] * y;
      mask[i * stride + j] = get_masked_weight_interintra(m);
    }
}

void vp9_generate_hard_mask_interintra(int mask_index, BLOCK_SIZE_TYPE sb_type,
                            int h, int w, uint8_t *mask, int stride) {
  int i, j;
  const int *a = get_mask_params_interintra(mask_index, sb_type, h, w);
  if (!a) return;
  for (i = 0; i < h; ++i)
    for (j = 0; j < w; ++j) {
      int x = (j - (a[2] * w) / 4);
      int y = (i - (a[3] * h) / 4);
      int m = a[0] * x + a[1] * y;
      mask[i * stride + j] = get_hard_mask_interintra(m);
    }
}
#endif

static void combine_interintra(MB_PREDICTION_MODE mode,
#if CONFIG_MASKED_COMPOUND
                               int use_masked_interintra,
                               int mask_index,
                               BLOCK_SIZE_TYPE bsize,
#endif
                               uint8_t *interpred,
                               int interstride,
                               uint8_t *intrapred,
                               int intrastride,
                               int bw, int bh) {
  static const int scale_bits = 8;
  static const int scale_max = 256;
  static const int scale_round = 127;
  static const int weights1d[64] = {
      128, 125, 122, 119, 116, 114, 111, 109,
      107, 105, 103, 101,  99,  97,  96,  94,
       93,  91,  90,  89,  88,  86,  85,  84,
       83,  82,  81,  81,  80,  79,  78,  78,
       77,  76,  76,  75,  75,  74,  74,  73,
       73,  72,  72,  71,  71,  71,  70,  70,
       70,  70,  69,  69,  69,  69,  68,  68,
       68,  68,  68,  67,  67,  67,  67,  67,
  };

  int size = MAX(bw, bh);
  int size_scale = (size >= 64 ? 1 :
                    size == 32 ? 2 :
                    size == 16 ? 4 :
                    size == 8  ? 8 : 16);
  int i, j;

#if CONFIG_MASKED_COMPOUND
  uint8_t mask[4096];
  if (use_masked_interintra && get_mask_bits_interintra(bsize))
    vp9_generate_masked_weight_interintra(mask_index, bsize, bh, bw, mask, bw);
#endif

  switch (mode) {
    case V_PRED:
      for (i = 0; i < bh; ++i) {
        for (j = 0; j < bw; ++j) {
          int k = i * interstride + j;
          int scale = weights1d[i * size_scale];
#if CONFIG_MASKED_COMPOUND
          int m = mask[i * bw + j];
          if (use_masked_interintra && get_mask_bits_interintra(bsize))
              interpred[k] = (intrapred[i * intrastride + j] * m +
                              interpred[k] *
                              ((1 << MASK_WEIGHT_BITS_INTERINTRA) - m) +
                              (1 << (MASK_WEIGHT_BITS_INTERINTRA - 1))) >>
                              MASK_WEIGHT_BITS_INTERINTRA;
          else
#endif
          interpred[k] =
              ((scale_max - scale) * interpred[k] +
                  scale * intrapred[i * intrastride + j] + scale_round)
                  >> scale_bits;
        }
      }
     break;

    case H_PRED:
      for (i = 0; i < bh; ++i) {
        for (j = 0; j < bw; ++j) {
          int k = i * interstride + j;
          int scale = weights1d[j * size_scale];
#if CONFIG_MASKED_COMPOUND
          int m = mask[i * bw + j];
          if (use_masked_interintra && get_mask_bits_interintra(bsize))
              interpred[k] = (intrapred[i * intrastride + j] * m +
                              interpred[k] *
                              ((1 << MASK_WEIGHT_BITS_INTERINTRA) - m) +
                              (1 << (MASK_WEIGHT_BITS_INTERINTRA - 1))) >>
                              MASK_WEIGHT_BITS_INTERINTRA;
          else
#endif
          interpred[k] =
              ((scale_max - scale) * interpred[k] +
                  scale * intrapred[i * intrastride + j] + scale_round)
                  >> scale_bits;
        }
      }
     break;

    case D63_PRED:
    case D117_PRED:
      for (i = 0; i < bh; ++i) {
        for (j = 0; j < bw; ++j) {
          int k = i * interstride + j;
          int scale = (weights1d[i * size_scale] * 3 +
                       weights1d[j * size_scale]) >> 2;
#if CONFIG_MASKED_COMPOUND
          int m = mask[i * bw + j];
          if (use_masked_interintra && get_mask_bits_interintra(bsize))
              interpred[k] = (intrapred[i * intrastride + j] * m +
                              interpred[k] *
                              ((1 << MASK_WEIGHT_BITS_INTERINTRA) - m) +
                              (1 << (MASK_WEIGHT_BITS_INTERINTRA - 1))) >>
                              MASK_WEIGHT_BITS_INTERINTRA;
          else
#endif
          interpred[k] =
              ((scale_max - scale) * interpred[k] +
                  scale * intrapred[i * intrastride + j] + scale_round)
                  >> scale_bits;
        }
      }
     break;

    case D27_PRED:
    case D153_PRED:
      for (i = 0; i < bh; ++i) {
        for (j = 0; j < bw; ++j) {
          int k = i * interstride + j;
          int scale = (weights1d[j * size_scale] * 3 +
                       weights1d[i * size_scale]) >> 2;
#if CONFIG_MASKED_COMPOUND
          int m = mask[i * bw + j];
          if (use_masked_interintra && get_mask_bits_interintra(bsize))
              interpred[k] = (intrapred[i * intrastride + j] * m +
                              interpred[k] *
                              ((1 << MASK_WEIGHT_BITS_INTERINTRA) - m) +
                              (1 << (MASK_WEIGHT_BITS_INTERINTRA - 1))) >>
                              MASK_WEIGHT_BITS_INTERINTRA;
          else
#endif
          interpred[k] =
              ((scale_max - scale) * interpred[k] +
                  scale * intrapred[i * intrastride + j] + scale_round)
                  >> scale_bits;
        }
      }
     break;

    case D135_PRED:
      for (i = 0; i < bh; ++i) {
        for (j = 0; j < bw; ++j) {
          int k = i * interstride + j;
          int scale = weights1d[(i < j ? i : j) * size_scale];
#if CONFIG_MASKED_COMPOUND
          int m = mask[i * bw + j];
          if (use_masked_interintra && get_mask_bits_interintra(bsize))
              interpred[k] = (intrapred[i * intrastride + j] * m +
                              interpred[k] *
                              ((1 << MASK_WEIGHT_BITS_INTERINTRA) - m) +
                              (1 << (MASK_WEIGHT_BITS_INTERINTRA - 1))) >>
                              MASK_WEIGHT_BITS_INTERINTRA;
          else
#endif
          interpred[k] =
              ((scale_max - scale) * interpred[k] +
                  scale * intrapred[i * intrastride + j] + scale_round)
                  >> scale_bits;
        }
      }
     break;

    case D45_PRED:
      for (i = 0; i < bh; ++i) {
        for (j = 0; j < bw; ++j) {
          int k = i * interstride + j;
          int scale = (weights1d[i * size_scale] +
                       weights1d[j * size_scale]) >> 1;
#if CONFIG_MASKED_COMPOUND
          int m = mask[i * bw + j];
          if (use_masked_interintra && get_mask_bits_interintra(bsize))
              interpred[k] = (intrapred[i * intrastride + j] * m +
                              interpred[k] *
                              ((1 << MASK_WEIGHT_BITS_INTERINTRA) - m) +
                              (1 << (MASK_WEIGHT_BITS_INTERINTRA - 1))) >>
                              MASK_WEIGHT_BITS_INTERINTRA;
          else
#endif
          interpred[k] =
              ((scale_max - scale) * interpred[k] +
                  scale * intrapred[i * intrastride + j] + scale_round)
                  >> scale_bits;
        }
      }
     break;

    case TM_PRED:
    case DC_PRED:
    default:
      for (i = 0; i < bh; ++i) {
        for (j = 0; j < bw; ++j) {
          int k = i * interstride + j;
#if CONFIG_MASKED_COMPOUND
          int m = mask[i * bw + j];
          if (use_masked_interintra && get_mask_bits_interintra(bsize))
              interpred[k] = (intrapred[i * intrastride + j] * m +
                              interpred[k] *
                              ((1 << MASK_WEIGHT_BITS_INTERINTRA) - m) +
                              (1 << (MASK_WEIGHT_BITS_INTERINTRA - 1))) >>
                              MASK_WEIGHT_BITS_INTERINTRA;
          else
#endif
          interpred[k] = (interpred[k] + intrapred[i * intrastride + j]) >> 1;
        }
      }
      break;
  }
}

// Break down rectangular intra prediction for joint spatio-temporal prediction
// into two square intra predictions.
static void build_intra_predictors_for_interintra(uint8_t *src, int src_stride,
                                           uint8_t *pred_ptr, int stride,
                                           MB_PREDICTION_MODE mode,
                                           int bw, int bh,
                                           int up_available, int left_available,
                                           int right_available) {
  if (bw == bh) {
    build_intra_predictors(src, src_stride, pred_ptr, stride,
                           mode, intra_size_log2_for_interintra(bw),
                           up_available, left_available, right_available);
  } else if (bw < bh) {
    uint8_t *src_bottom = src + bw * src_stride;
    uint8_t *pred_ptr_bottom = pred_ptr + bw * stride;
    build_intra_predictors(src, src_stride, pred_ptr, stride,
                           mode, intra_size_log2_for_interintra(bw),
                           up_available, left_available, right_available);
    build_intra_predictors_for_2nd_block_interintra(src_bottom, src_stride,
                                                    pred_ptr_bottom, stride,
                                       mode, intra_size_log2_for_interintra(bw),
                                       1, left_available, right_available, 1);
  } else {
    uint8_t *src_right = src + bh;
    uint8_t *pred_ptr_right = pred_ptr + bh;
    build_intra_predictors(src, src_stride, pred_ptr, stride,
                           mode, intra_size_log2_for_interintra(bh),
                           up_available, left_available, right_available);
    build_intra_predictors_for_2nd_block_interintra(src_right, src_stride,
                                                    pred_ptr_right, stride,
                                       mode, intra_size_log2_for_interintra(bh),
                                       up_available, 1, right_available, 0);
  }
}

void vp9_build_interintra_predictors_sby(MACROBLOCKD *xd,
                                         uint8_t *ypred,
                                         int ystride,
                                         BLOCK_SIZE_TYPE bsize) {
  const struct macroblockd_plane* const pd = &xd->plane[0];
  const int bw = plane_block_width(bsize, pd);
  const int bh = plane_block_height(bsize, pd);
  uint8_t intrapredictor[4096];
  build_intra_predictors_for_interintra(
      xd->plane[0].dst.buf, xd->plane[0].dst.stride,
      intrapredictor, bw,
      xd->mode_info_context->mbmi.interintra_mode, bw, bh,
      xd->up_available, xd->left_available, xd->right_available);
  combine_interintra(xd->mode_info_context->mbmi.interintra_mode,
#if CONFIG_MASKED_COMPOUND
                     xd->mode_info_context->mbmi.use_masked_interintra,
                     xd->mode_info_context->mbmi.interintra_mask_index,
                     bsize,
#endif
                     ypred, ystride, intrapredictor, bw, bw, bh);
}

void vp9_build_interintra_predictors_sbuv(MACROBLOCKD *xd,
                                          uint8_t *upred,
                                          uint8_t *vpred,
                                          int uvstride,
                                          BLOCK_SIZE_TYPE bsize) {
  int bwl = b_width_log2(bsize), bw = 2 << bwl;
  int bhl = b_height_log2(bsize), bh = 2 << bhl;
  uint8_t uintrapredictor[1024];
  uint8_t vintrapredictor[1024];
  build_intra_predictors_for_interintra(
      xd->plane[1].dst.buf, xd->plane[1].dst.stride,
      uintrapredictor, bw,
      xd->mode_info_context->mbmi.interintra_uv_mode, bw, bh,
      xd->up_available, xd->left_available, xd->right_available);
  build_intra_predictors_for_interintra(
      xd->plane[2].dst.buf, xd->plane[1].dst.stride,
      vintrapredictor, bw,
      xd->mode_info_context->mbmi.interintra_uv_mode, bw, bh,
      xd->up_available, xd->left_available, xd->right_available);
  combine_interintra(xd->mode_info_context->mbmi.interintra_uv_mode,
#if CONFIG_MASKED_COMPOUND
                     xd->mode_info_context->mbmi.use_masked_interintra,
                     xd->mode_info_context->mbmi.interintra_uv_mask_index,
                     bsize,
#endif
                     upred, uvstride, uintrapredictor, bw, bw, bh);
  combine_interintra(xd->mode_info_context->mbmi.interintra_uv_mode,
#if CONFIG_MASKED_COMPOUND
                     xd->mode_info_context->mbmi.use_masked_interintra,
                     xd->mode_info_context->mbmi.interintra_uv_mask_index,
                     bsize,
#endif
                     vpred, uvstride, vintrapredictor, bw, bw, bh);
}

void vp9_build_interintra_predictors(MACROBLOCKD *xd,
                                     uint8_t *ypred,
                                     uint8_t *upred,
                                     uint8_t *vpred,
                                     int ystride, int uvstride,
                                     BLOCK_SIZE_TYPE bsize) {
  vp9_build_interintra_predictors_sby(xd, ypred, ystride, bsize);
  vp9_build_interintra_predictors_sbuv(xd, upred, vpred, uvstride, bsize);
}
#endif
