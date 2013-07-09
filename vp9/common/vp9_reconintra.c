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


void vp9_d27_predictor_c(uint8_t *ypred_ptr, int y_stride, int bs,
                         uint8_t *yabove_row, uint8_t *yleft_col) {
  int r, c;
  // first column
  for (r = 0; r < bs - 1; ++r) {
      ypred_ptr[r * y_stride] = ROUND_POWER_OF_TWO(yleft_col[r] +
                                                   yleft_col[r + 1], 1);
  }
  ypred_ptr[(bs - 1) * y_stride] = yleft_col[bs - 1];
  ypred_ptr++;
  // second column
  for (r = 0; r < bs - 2; ++r) {
      ypred_ptr[r * y_stride] = ROUND_POWER_OF_TWO(yleft_col[r] +
                                                   yleft_col[r + 1] * 2 +
                                                   yleft_col[r + 2], 2);
  }
  ypred_ptr[(bs - 2) * y_stride] = ROUND_POWER_OF_TWO(yleft_col[bs - 2] +
                                                      yleft_col[bs - 1] * 3,
                                                      2);
  ypred_ptr[(bs - 1) * y_stride] = yleft_col[bs - 1];
  ypred_ptr++;

  // rest of last row
  for (c = 0; c < bs - 2; ++c) {
    ypred_ptr[(bs - 1) * y_stride + c] = yleft_col[bs - 1];
  }

  for (r = bs - 2; r >= 0; --r) {
    for (c = 0; c < bs - 2; ++c) {
      ypred_ptr[r * y_stride + c] = ypred_ptr[(r + 1) * y_stride + c - 2];
    }
  }
}

void vp9_d63_predictor_c(uint8_t *ypred_ptr, int y_stride, int bs,
                         uint8_t *yabove_row, uint8_t *yleft_col) {
  int r, c;
  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      if (r & 1) {
        ypred_ptr[c] = ROUND_POWER_OF_TWO(yabove_row[r/2 + c] +
                                          yabove_row[r/2 + c + 1] * 2 +
                                          yabove_row[r/2 + c + 2], 2);
      } else {
        ypred_ptr[c] =ROUND_POWER_OF_TWO(yabove_row[r/2 + c] +
                                         yabove_row[r/2+ c + 1], 1);
      }
    }
    ypred_ptr += y_stride;
  }
}

void vp9_d45_predictor_c(uint8_t *ypred_ptr, int y_stride, int bs,
                         uint8_t *yabove_row, uint8_t *yleft_col) {
  int r, c;
  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      if (r + c + 2 < bs * 2)
        ypred_ptr[c] = ROUND_POWER_OF_TWO(yabove_row[r + c] +
                                          yabove_row[r + c + 1] * 2 +
                                          yabove_row[r + c + 2], 2);
      else
        ypred_ptr[c] = yabove_row[bs * 2 - 1];
    }
    ypred_ptr += y_stride;
  }
}

void vp9_d117_predictor_c(uint8_t *ypred_ptr, int y_stride, int bs,
                          uint8_t *yabove_row, uint8_t *yleft_col) {
  int r, c;
  // first row
  for (c = 0; c < bs; c++)
    ypred_ptr[c] = ROUND_POWER_OF_TWO(yabove_row[c - 1] + yabove_row[c], 1);
  ypred_ptr += y_stride;

  // second row
  ypred_ptr[0] = ROUND_POWER_OF_TWO(yleft_col[0] +
                                    yabove_row[-1] * 2 +
                                    yabove_row[0], 2);
  for (c = 1; c < bs; c++)
    ypred_ptr[c] = ROUND_POWER_OF_TWO(yabove_row[c - 2] +
                                      yabove_row[c - 1] * 2 +
                                      yabove_row[c], 2);
  ypred_ptr += y_stride;

  // the rest of first col
  ypred_ptr[0] = ROUND_POWER_OF_TWO(yabove_row[-1] +
                                    yleft_col[0] * 2 +
                                    yleft_col[1], 2);
  for (r = 3; r < bs; ++r)
    ypred_ptr[(r-2) * y_stride] = ROUND_POWER_OF_TWO(yleft_col[r - 3] +
                                                     yleft_col[r - 2] * 2 +
                                                     yleft_col[r - 1], 2);
  // the rest of the block
  for (r = 2; r < bs; ++r) {
    for (c = 1; c < bs; c++)
      ypred_ptr[c] = ypred_ptr[-2 * y_stride + c - 1];
    ypred_ptr += y_stride;
  }
}


void vp9_d135_predictor_c(uint8_t *ypred_ptr, int y_stride, int bs,
                          uint8_t *yabove_row, uint8_t *yleft_col) {
  int r, c;
  ypred_ptr[0] = ROUND_POWER_OF_TWO(yleft_col[0] +
                                    yabove_row[-1] * 2 +
                                    yabove_row[0], 2);
  for (c = 1; c < bs; c++)
    ypred_ptr[c] = ROUND_POWER_OF_TWO(yabove_row[c - 2] +
                                      yabove_row[c - 1] * 2 +
                                      yabove_row[c], 2);

  ypred_ptr[y_stride] = ROUND_POWER_OF_TWO(yabove_row[-1] +
                                           yleft_col[0] * 2 +
                                           yleft_col[1], 2);
  for (r = 2; r < bs; ++r)
    ypred_ptr[r * y_stride] = ROUND_POWER_OF_TWO(yleft_col[r - 2] +
                                                 yleft_col[r - 1] * 2 +
                                                 yleft_col[r], 2);

  ypred_ptr += y_stride;
  for (r = 1; r < bs; ++r) {
    for (c = 1; c < bs; c++)
      ypred_ptr[c] = ypred_ptr[-y_stride + c - 1];
    ypred_ptr += y_stride;
  }
}

void vp9_d153_predictor_c(uint8_t *ypred_ptr, int y_stride, int bs,
                          uint8_t *yabove_row, uint8_t *yleft_col) {
  int r, c;
  ypred_ptr[0] = ROUND_POWER_OF_TWO(yabove_row[-1] + yleft_col[0], 1);
  for (r = 1; r < bs; r++)
    ypred_ptr[r * y_stride] =
        ROUND_POWER_OF_TWO(yleft_col[r - 1] + yleft_col[r], 1);
  ypred_ptr++;

  ypred_ptr[0] = ROUND_POWER_OF_TWO(yleft_col[0] +
                                    yabove_row[-1] * 2 +
                                    yabove_row[0], 2);
  ypred_ptr[y_stride] = ROUND_POWER_OF_TWO(yabove_row[-1] +
                                           yleft_col[0] * 2 +
                                           yleft_col[1], 2);
  for (r = 2; r < bs; r++)
    ypred_ptr[r * y_stride] = ROUND_POWER_OF_TWO(yleft_col[r - 2] +
                                                 yleft_col[r - 1] * 2 +
                                                 yleft_col[r], 2);
  ypred_ptr++;

  for (c = 0; c < bs - 2; c++)
    ypred_ptr[c] = ROUND_POWER_OF_TWO(yabove_row[c - 1] +
                                      yabove_row[c] * 2 +
                                      yabove_row[c + 1], 2);
  ypred_ptr += y_stride;
  for (r = 1; r < bs; ++r) {
    for (c = 0; c < bs - 2; c++)
      ypred_ptr[c] = ypred_ptr[-y_stride + c - 2];
    ypred_ptr += y_stride;
  }
}

void vp9_v_predictor_c(uint8_t *ypred_ptr, int y_stride, int bs,
                       uint8_t *yabove_row, uint8_t *yleft_col) {
  int r;

  for (r = 0; r < bs; r++) {
    vpx_memcpy(ypred_ptr, yabove_row, bs);
    ypred_ptr += y_stride;
  }
}

void vp9_h_predictor_c(uint8_t *ypred_ptr, int y_stride, int bs,
                       uint8_t *yabove_row, uint8_t *yleft_col) {
  int r;

  for (r = 0; r < bs; r++) {
    vpx_memset(ypred_ptr, yleft_col[r], bs);
    ypred_ptr += y_stride;
  }
}

void vp9_tm_predictor_c(uint8_t *ypred_ptr, int y_stride, int bs,
                        uint8_t *yabove_row, uint8_t *yleft_col) {
  int r, c;
  int ytop_left = yabove_row[-1];

  for (r = 0; r < bs; r++) {
    for (c = 0; c < bs; c++)
      ypred_ptr[c] = clip_pixel(yleft_col[r] + yabove_row[c] - ytop_left);
    ypred_ptr += y_stride;
  }
}

void vp9_dc_128_predictor_c(uint8_t *ypred_ptr, int y_stride, int bs,
                            uint8_t *yabove_row, uint8_t *yleft_col) {
  int r;

  for (r = 0; r < bs; r++) {
    vpx_memset(ypred_ptr, 128, bs);
    ypred_ptr += y_stride;
  }
}

void vp9_dc_left_predictor_c(uint8_t *ypred_ptr, int y_stride, int bs,
                             uint8_t *yabove_row, uint8_t *yleft_col) {
  int i, r;
  int expected_dc = 128;
  int average = 0;
  const int count = bs;

  for (i = 0; i < bs; i++)
    average += yleft_col[i];
  expected_dc = (average + (count >> 1)) / count;

  for (r = 0; r < bs; r++) {
    vpx_memset(ypred_ptr, expected_dc, bs);
    ypred_ptr += y_stride;
  }
}

void vp9_dc_top_predictor_c(uint8_t *ypred_ptr, int y_stride, int bs,
                            uint8_t *yabove_row, uint8_t *yleft_col) {
  int i, r;
  int expected_dc = 128;
  int average = 0;
  const int count = bs;

  for (i = 0; i < bs; i++)
    average += yabove_row[i];
  expected_dc = (average + (count >> 1)) / count;

  for (r = 0; r < bs; r++) {
    vpx_memset(ypred_ptr, expected_dc, bs);
    ypred_ptr += y_stride;
  }
}

void vp9_dc_predictor_c(uint8_t *ypred_ptr, int y_stride, int bs,
                        uint8_t *yabove_row, uint8_t *yleft_col) {
  int i, r;
  int expected_dc = 128;
  int average = 0;
  const int count = 2 * bs;

  for (i = 0; i < bs; i++)
    average += yabove_row[i];
  for (i = 0; i < bs; i++)
    average += yleft_col[i];
  expected_dc = (average + (count >> 1)) / count;

  for (r = 0; r < bs; r++) {
    vpx_memset(ypred_ptr, expected_dc, bs);
    ypred_ptr += y_stride;
  }
}

typedef void (*intra_pred_fn)(uint8_t *ypred_ptr, int y_stride, int bs,
                              uint8_t *yabove_row, uint8_t *yleft_col);

static void build_intra_predictors(uint8_t *src, int src_stride,
                                   uint8_t *ypred_ptr, int y_stride,
                                   MB_PREDICTION_MODE mode, int bs,
                                   int up_available, int left_available,
                                   int right_available) {
  int i;
  DECLARE_ALIGNED_ARRAY(16, uint8_t, yleft_col, 64);
  DECLARE_ALIGNED_ARRAY(16, uint8_t, yabove_data, 128 + 16);
  uint8_t *yabove_row = yabove_data + 16;
  static const intra_pred_fn pred[VP9_INTRA_MODES] = {
    NULL, vp9_v_predictor, vp9_h_predictor, vp9_d45_predictor,
    vp9_d135_predictor, vp9_d117_predictor, vp9_d153_predictor,
    vp9_d27_predictor, vp9_d63_predictor, vp9_tm_predictor
  };

  // 127 127 127 .. 127 127 127 127 127 127
  // 129  A   B  ..  Y   Z
  // 129  C   D  ..  W   X
  // 129  E   F  ..  U   V
  // 129  G   H  ..  S   T   T   T   T   T
  // ..

  if (left_available) {
    for (i = 0; i < bs; i++)
      yleft_col[i] = src[i * src_stride - 1];
  } else {
    vpx_memset(yleft_col, 129, bs);
  }

  if (up_available) {
    uint8_t *yabove_ptr = src - src_stride;
    if (bs == 4 && right_available && left_available) {
      yabove_row = yabove_ptr;
    } else {
      vpx_memcpy(yabove_row, yabove_ptr, bs);
      if (bs == 4 && right_available)
        vpx_memcpy(yabove_row + bs, yabove_ptr + bs, bs);
      else
        vpx_memset(yabove_row + bs, yabove_row[bs - 1], bs);
      yabove_row[-1] = left_available ? yabove_ptr[-1] : 129;
    }
  } else {
    vpx_memset(yabove_row, 127, bs * 2);
    yabove_row[-1] = 127;
  }

  if (mode == DC_PRED) {
    if (left_available) {
      if (up_available) {
        vp9_dc_predictor(ypred_ptr, y_stride, bs, yabove_row, yleft_col);
      } else {
        vp9_dc_left_predictor(ypred_ptr, y_stride, bs, yabove_row, yleft_col);
      }
    } else if (up_available) {
      vp9_dc_top_predictor(ypred_ptr, y_stride, bs, yabove_row, yleft_col);
    } else {
      vp9_dc_128_predictor(ypred_ptr, y_stride, bs, yabove_row, yleft_col);
    }
  } else {
    pred[mode](ypred_ptr, y_stride, bs, yabove_row, yleft_col);
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
  const int txfm_block_size = 4 << tx_size;

  assert(bwl >= 0);
  build_intra_predictors(reference, ref_stride,
                         predictor, pre_stride,
                         mode,
                         txfm_block_size,
                         have_top, have_left,
                         have_right);
}
