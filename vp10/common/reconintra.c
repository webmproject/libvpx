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

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"

#if CONFIG_VP9_HIGHBITDEPTH
#include "vpx_dsp/vpx_dsp_common.h"
#endif  // CONFIG_VP9_HIGHBITDEPTH
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem.h"
#include "vpx_ports/vpx_once.h"

#include "vp10/common/reconintra.h"
#include "vp10/common/onyxc_int.h"

enum {
  NEED_LEFT = 1 << 1,
  NEED_ABOVE = 1 << 2,
  NEED_ABOVERIGHT = 1 << 3,
  NEED_ABOVELEFT = 1 << 4,
  NEED_BOTTOMLEFT = 1 << 5,
};

static const uint8_t extend_modes[INTRA_MODES] = {
  NEED_ABOVE | NEED_LEFT,                   // DC
  NEED_ABOVE,                               // V
  NEED_LEFT,                                // H
  NEED_ABOVE | NEED_ABOVERIGHT,             // D45
  NEED_LEFT | NEED_ABOVE | NEED_ABOVELEFT,  // D135
  NEED_LEFT | NEED_ABOVE | NEED_ABOVELEFT,  // D117
  NEED_LEFT | NEED_ABOVE | NEED_ABOVELEFT,  // D153
  NEED_LEFT | NEED_BOTTOMLEFT,              // D207
  NEED_ABOVE | NEED_ABOVERIGHT,             // D63
  NEED_LEFT | NEED_ABOVE | NEED_ABOVELEFT,  // TM
};

static const uint8_t orders_64x64[1] = { 0 };
static const uint8_t orders_64x32[2] = { 0, 1 };
static const uint8_t orders_32x64[2] = { 0, 1 };
static const uint8_t orders_32x32[4] = {
  0, 1,
  2, 3,
};
static const uint8_t orders_32x16[8] = {
  0, 2,
  1, 3,
  4, 6,
  5, 7,
};
static const uint8_t orders_16x32[8] = {
  0, 1, 2, 3,
  4, 5, 6, 7,
};
static const uint8_t orders_16x16[16] = {
  0,   1,  4,  5,
  2,   3,  6,  7,
  8,   9, 12, 13,
  10, 11, 14, 15,
};
static const uint8_t orders_16x8[32] = {
  0,   2,  8, 10,
  1,   3,  9, 11,
  4,   6, 12, 14,
  5,   7, 13, 15,
  16, 18, 24, 26,
  17, 19, 25, 27,
  20, 22, 28, 30,
  21, 23, 29, 31,
};
static const uint8_t orders_8x16[32] = {
  0,   1,  2,  3,  8,  9, 10, 11,
  4,   5,  6,  7, 12, 13, 14, 15,
  16, 17, 18, 19, 24, 25, 26, 27,
  20, 21, 22, 23, 28, 29, 30, 31,
};
static const uint8_t orders_8x8[64] = {
  0,   1,  4,  5, 16, 17, 20, 21,
  2,   3,  6,  7, 18, 19, 22, 23,
  8,   9, 12, 13, 24, 25, 28, 29,
  10, 11, 14, 15, 26, 27, 30, 31,
  32, 33, 36, 37, 48, 49, 52, 53,
  34, 35, 38, 39, 50, 51, 54, 55,
  40, 41, 44, 45, 56, 57, 60, 61,
  42, 43, 46, 47, 58, 59, 62, 63,
};
static const uint8_t *const orders[BLOCK_SIZES] = {
  orders_8x8, orders_8x8, orders_8x8, orders_8x8,
  orders_8x16, orders_16x8, orders_16x16,
  orders_16x32, orders_32x16, orders_32x32,
  orders_32x64, orders_64x32, orders_64x64,
};

static int vp10_has_right(BLOCK_SIZE bsize, int mi_row, int mi_col,
                          int right_available,
                          TX_SIZE txsz, int y, int x, int ss_x) {
  const int wl = mi_width_log2_lookup[bsize];
  const int w = VPXMAX(num_4x4_blocks_wide_lookup[bsize] >> ss_x, 1);
  const int step = 1 << txsz;

  // Handle block size 4x8 and 4x4
  if (ss_x == 0 && num_4x4_blocks_wide_lookup[bsize] < 2 && x == 0)
    return 1;

  if (y == 0) {
    const int hl = mi_height_log2_lookup[bsize];
    const uint8_t *order = orders[bsize];
    int my_order, tr_order;

    if (x + step < w)
      return 1;

    mi_row = (mi_row & 7) >> hl;
    mi_col = (mi_col & 7) >> wl;

    if (mi_row == 0)
      return right_available;

    if (((mi_col + 1) << wl) >= 8)
      return 0;

    my_order = order[((mi_row + 0) << (3 - wl)) + mi_col + 0];
    tr_order = order[((mi_row - 1) << (3 - wl)) + mi_col + 1];

    return my_order > tr_order && right_available;
  } else {
    return x + step < w;
  }
}

static int vp10_has_bottom(BLOCK_SIZE bsize, int mi_row, int mi_col,
                           int bottom_available, TX_SIZE txsz,
                           int y, int x, int ss_y) {
  if (x == 0) {
    const int wl = mi_width_log2_lookup[bsize];
    const int hl = mi_height_log2_lookup[bsize];
    const int h = 1 << (hl + 1 - ss_y);
    const int step = 1 << txsz;
    const uint8_t *order = orders[bsize];
    int my_order, bl_order;

    // Handle block size 8x4 and 4x4
    if (ss_y == 0 && num_4x4_blocks_high_lookup[bsize] < 2 && y == 0)
      return 1;

    if (y + step < h)
      return 1;

    mi_row = (mi_row & 7) >> hl;
    mi_col = (mi_col & 7) >> wl;

    if (mi_col == 0)
      return bottom_available &&
             (mi_row << (hl + !ss_y)) + y + step < (8 << !ss_y);

    if (((mi_row + 1) << hl) >= 8)
      return 0;

    my_order = order[((mi_row + 0) << (3 - wl)) + mi_col + 0];
    bl_order = order[((mi_row + 1) << (3 - wl)) + mi_col - 1];

    return bl_order < my_order && bottom_available;
  } else {
    return 0;
  }
}

typedef void (*intra_pred_fn)(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left);

static intra_pred_fn pred[INTRA_MODES][TX_SIZES];
static intra_pred_fn dc_pred[2][2][TX_SIZES];

#if CONFIG_VP9_HIGHBITDEPTH
typedef void (*intra_high_pred_fn)(uint16_t *dst, ptrdiff_t stride,
                                   const uint16_t *above, const uint16_t *left,
                                   int bd);
static intra_high_pred_fn pred_high[INTRA_MODES][4];
static intra_high_pred_fn dc_pred_high[2][2][4];
#endif  // CONFIG_VP9_HIGHBITDEPTH

static void vp10_init_intra_predictors_internal(void) {
#define INIT_NO_4X4(p, type) \
  p[TX_8X8] = vpx_##type##_predictor_8x8; \
  p[TX_16X16] = vpx_##type##_predictor_16x16; \
  p[TX_32X32] = vpx_##type##_predictor_32x32

#define INIT_ALL_SIZES(p, type) \
  p[TX_4X4] = vpx_##type##_predictor_4x4; \
  INIT_NO_4X4(p, type)

  INIT_ALL_SIZES(pred[V_PRED], v);
  INIT_ALL_SIZES(pred[H_PRED], h);
  INIT_ALL_SIZES(pred[D207_PRED], d207e);
  INIT_ALL_SIZES(pred[D45_PRED], d45e);
  INIT_ALL_SIZES(pred[D63_PRED], d63e);
  INIT_ALL_SIZES(pred[D117_PRED], d117);
  INIT_ALL_SIZES(pred[D135_PRED], d135);
  INIT_ALL_SIZES(pred[D153_PRED], d153);
  INIT_ALL_SIZES(pred[TM_PRED], tm);

  INIT_ALL_SIZES(dc_pred[0][0], dc_128);
  INIT_ALL_SIZES(dc_pred[0][1], dc_top);
  INIT_ALL_SIZES(dc_pred[1][0], dc_left);
  INIT_ALL_SIZES(dc_pred[1][1], dc);

#if CONFIG_VP9_HIGHBITDEPTH
  INIT_ALL_SIZES(pred_high[V_PRED], highbd_v);
  INIT_ALL_SIZES(pred_high[H_PRED], highbd_h);
  INIT_ALL_SIZES(pred_high[D207_PRED], highbd_d207e);
  INIT_ALL_SIZES(pred_high[D45_PRED], highbd_d45e);
  INIT_ALL_SIZES(pred_high[D63_PRED], highbd_d63e);
  INIT_ALL_SIZES(pred_high[D117_PRED], highbd_d117);
  INIT_ALL_SIZES(pred_high[D135_PRED], highbd_d135);
  INIT_ALL_SIZES(pred_high[D153_PRED], highbd_d153);
  INIT_ALL_SIZES(pred_high[TM_PRED], highbd_tm);

  INIT_ALL_SIZES(dc_pred_high[0][0], highbd_dc_128);
  INIT_ALL_SIZES(dc_pred_high[0][1], highbd_dc_top);
  INIT_ALL_SIZES(dc_pred_high[1][0], highbd_dc_left);
  INIT_ALL_SIZES(dc_pred_high[1][1], highbd_dc);
#endif  // CONFIG_VP9_HIGHBITDEPTH

#undef intra_pred_allsizes
}

#if CONFIG_EXT_INTRA
#define PI 3.14159265
#define FILTER_INTRA_PREC_BITS 10
#define FILTER_INTRA_ROUND_VAL 511

static const uint8_t ext_intra_extend_modes[FILTER_INTRA_MODES] = {
  NEED_LEFT | NEED_ABOVE,      // FILTER_DC
  NEED_LEFT | NEED_ABOVE,      // FILTER_V
  NEED_LEFT | NEED_ABOVE,      // FILTER_H
  NEED_LEFT | NEED_ABOVE,      // FILTER_D45
  NEED_LEFT | NEED_ABOVE,      // FILTER_D135
  NEED_LEFT | NEED_ABOVE,      // FILTER_D117
  NEED_LEFT | NEED_ABOVE,      // FILTER_D153
  NEED_LEFT | NEED_ABOVE,      // FILTER_D207
  NEED_LEFT | NEED_ABOVE,      // FILTER_D63
  NEED_LEFT | NEED_ABOVE,      // FILTER_TM
};

static int intra_subpel_interp(int base, int shift, const uint8_t *ref,
                               int ref_start_idx, int ref_end_idx,
                               INTRA_FILTER filter_type) {
  int val, k, idx, filter_idx = 0;
  const int16_t *filter = NULL;

  if (filter_type == INTRA_FILTER_LINEAR) {
    val = ref[base] * (256 - shift) + ref[base + 1] * shift;
    val = ROUND_POWER_OF_TWO(val, 8);
  } else {
    filter_idx = ROUND_POWER_OF_TWO(shift, 8 - SUBPEL_BITS);
    filter = vp10_intra_filter_kernels[filter_type][filter_idx];

    if (filter_idx < (1 << SUBPEL_BITS)) {
      val = 0;
      for (k = 0; k < SUBPEL_TAPS; ++k) {
        idx = base + 1 - (SUBPEL_TAPS / 2) + k;
        idx = VPXMAX(VPXMIN(idx, ref_end_idx), ref_start_idx);
        val += ref[idx] * filter[k];
      }
      val = ROUND_POWER_OF_TWO(val, FILTER_BITS);
    } else {
      val = ref[base + 1];
    }
  }

  return val;
}

// Directional prediction, zone 1: 0 < angle < 90
static void dr_prediction_z1(uint8_t *dst, ptrdiff_t stride, int bs,
                             const uint8_t *above, const uint8_t *left,
                             int dx, int dy, INTRA_FILTER filter_type) {
  int r, c, x, base, shift, val;

  (void)left;
  (void)dy;
  assert(dy == 1);
  assert(dx < 0);

  if (filter_type != INTRA_FILTER_LINEAR) {
    const int pad_size = SUBPEL_TAPS >> 1;
    int len;
    DECLARE_ALIGNED(16, uint8_t, buf[SUBPEL_SHIFTS][64]);
    DECLARE_ALIGNED(16, uint8_t, src[64 + SUBPEL_TAPS]);
    uint8_t flags[SUBPEL_SHIFTS];

    memset(flags, 0, SUBPEL_SHIFTS * sizeof(flags[0]));
    memset(src, above[0], pad_size * sizeof(above[0]));
    memcpy(src + pad_size, above, 2 * bs * sizeof(above[0]));
    memset(src + pad_size + 2 * bs, above[2 * bs - 1],
           pad_size * sizeof(above[0]));
    flags[0] = 1;
    x = -dx;
    for (r = 0; r < bs; ++r, dst += stride, x -= dx) {
      base = x >> 8;
      shift = x & 0xFF;
      shift = ROUND_POWER_OF_TWO(shift, 8 - SUBPEL_BITS);
      if (shift == SUBPEL_SHIFTS) {
        base += 1;
        shift = 0;
      }
      len = VPXMIN(bs, 2 * bs - 1 - base);
      if (len <= 0) {
        int i;
        for (i = r; i < bs; ++i) {
          memset(dst, above[2 * bs - 1], bs * sizeof(dst[0]));
          dst += stride;
        }
        return;
      }

      if (len <= (bs >> 1) && !flags[shift]) {
        base = x >> 8;
        shift = x & 0xFF;
        for (c = 0; c < len; ++c) {
          val = intra_subpel_interp(base, shift, above, 0, 2 * bs - 1,
                                    filter_type);
          dst[c] = clip_pixel(val);
          ++base;
        }
      } else {
        if (!flags[shift]) {
          vpx_convolve8_horiz(src + pad_size, 2 * bs, buf[shift], 2 * bs,
                              vp10_intra_filter_kernels[filter_type][shift], 16,
                              NULL, 16, 2 * bs, 2 * bs < 16 ? 2 : 1);
          flags[shift] = 1;
        }
        memcpy(dst, shift == 0 ? src + pad_size + base : &buf[shift][base],
            len * sizeof(dst[0]));
      }

      if (len < bs)
        memset(dst + len, above[2 * bs - 1], (bs - len) * sizeof(dst[0]));
    }
    return;
  }

  // For linear filter, C code is faster.
  x = -dx;
  for (r = 0; r < bs; ++r, dst += stride, x -= dx) {
    base = x >> 8;
    shift = x & 0xFF;

    if (base >= 2 * bs - 1) {
      int i;
      for (i = r; i < bs; ++i) {
        memset(dst, above[2 * bs - 1], bs * sizeof(dst[0]));
        dst += stride;
      }
      return;
    }

    for (c = 0; c < bs; ++c, ++base) {
      if (base < 2 * bs - 1) {
        val = above[base] * (256 - shift) + above[base + 1] * shift;
        val = ROUND_POWER_OF_TWO(val, 8);
        dst[c] = clip_pixel(val);
      } else {
        dst[c] = above[2 * bs - 1];
      }
    }
  }
}

// Directional prediction, zone 2: 90 < angle < 180
static void dr_prediction_z2(uint8_t *dst, ptrdiff_t stride, int bs,
                             const uint8_t *above, const uint8_t *left,
                             int dx, int dy, INTRA_FILTER filter_type) {
  int r, c, x, y, shift1, shift2, val, base1, base2;

  assert(dx > 0);
  assert(dy > 0);

  x = -dx;
  for (r = 0; r < bs; ++r, x -= dx, dst += stride) {
    base1 = x >> 8;
    y = (r << 8) - dy;
    for (c = 0; c < bs; ++c, ++base1, y -= dy) {
      if (base1 >= -1) {
        shift1 = x & 0xFF;
        val = intra_subpel_interp(base1, shift1, above, -1, bs - 1,
                                  filter_type);
      } else {
        base2 = y >> 8;
        if (base2 >= 0) {
          shift2 = y & 0xFF;
          val = intra_subpel_interp(base2, shift2, left, 0, bs - 1,
                                    filter_type);
        } else {
          val = left[0];
        }
      }
      dst[c] = clip_pixel(val);
    }
  }
}

// Directional prediction, zone 3: 180 < angle < 270
static void dr_prediction_z3(uint8_t *dst, ptrdiff_t stride, int bs,
                             const uint8_t *above, const uint8_t *left,
                             int dx, int dy, INTRA_FILTER filter_type) {
  int r, c, y, base, shift, val;

  (void)above;
  (void)dx;

  assert(dx == 1);
  assert(dy < 0);

  if (filter_type != INTRA_FILTER_LINEAR) {
    const int pad_size = SUBPEL_TAPS >> 1;
    int len, i;
    DECLARE_ALIGNED(16, uint8_t, buf[64][4 * SUBPEL_SHIFTS]);
    DECLARE_ALIGNED(16, uint8_t, src[(64 + SUBPEL_TAPS) * 4]);
    uint8_t flags[SUBPEL_SHIFTS];

    memset(flags, 0, SUBPEL_SHIFTS * sizeof(flags[0]));
    for (i = 0; i < pad_size; ++i)
      src[4 * i] = left[0];
    for (i = 0; i < 2 * bs; ++i)
      src[4 * (i + pad_size)] = left[i];
    for (i = 0; i < pad_size; ++i)
      src[4 * (i + 2 * bs + pad_size)] = left[2 * bs - 1];
    flags[0] = 1;
    y = -dy;
    for (c = 0; c < bs; ++c, y -= dy) {
      base = y >> 8;
      shift = y & 0xFF;
      shift = ROUND_POWER_OF_TWO(shift, 8 - SUBPEL_BITS);
      if (shift == SUBPEL_SHIFTS) {
        base += 1;
        shift = 0;
      }
      len = VPXMIN(bs, 2 * bs - 1 - base);

      if (len <= 0) {
        for (r = 0; r < bs; ++r) {
          dst[r * stride + c] = left[ 2 * bs - 1];
        }
        continue;
      }

      if (len <= (bs >> 1) && !flags[shift]) {
        base = y >> 8;
        shift = y & 0xFF;
        for (r = 0; r < len; ++r) {
          val = intra_subpel_interp(base, shift, left, 0, 2 * bs - 1,
                                    filter_type);
          dst[r * stride + c] = clip_pixel(val);
          ++base;
        }
      } else {
        if (!flags[shift]) {
          vpx_convolve8_vert(src + 4 * pad_size, 4,
                             buf[0] + 4 * shift, 4 * SUBPEL_SHIFTS, NULL, 16,
                             vp10_intra_filter_kernels[filter_type][shift], 16,
                             2 * bs < 16 ? 4 : 4, 2 * bs);
          flags[shift] = 1;
        }

        if (shift == 0) {
          for (r = 0; r < len; ++r) {
            dst[r * stride + c] = left[r + base];
          }
        } else {
          for (r = 0; r < len; ++r) {
            dst[r * stride + c] = buf[r + base][4 * shift];
          }
        }
      }

      if (len < bs) {
        for (r = len; r < bs; ++r) {
          dst[r * stride + c] = left[ 2 * bs - 1];
        }
      }
    }
    return;
  }

  // For linear filter, C code is faster.
  y = -dy;
  for (c = 0; c < bs; ++c, y -= dy) {
    base = y >> 8;
    shift = y & 0xFF;

    for (r = 0; r < bs; ++r, ++base) {
      if (base < 2 * bs - 1) {
        val = left[base] * (256 - shift) + left[base + 1] * shift;
        val = ROUND_POWER_OF_TWO(val, 8);
        dst[r * stride + c] = clip_pixel(val);
      } else {
        for (; r < bs; ++r)
          dst[r * stride + c] = left[2 * bs - 1];
        break;
      }
    }
  }
}

static void dr_predictor(uint8_t *dst, ptrdiff_t stride, TX_SIZE tx_size,
                         const uint8_t *above, const uint8_t *left, int angle,
                         INTRA_FILTER filter_type) {
  double t = 0;
  int dx, dy;
  int bs = 4 << tx_size;

  if (angle != 90 && angle != 180)
    t = tan(angle * PI / 180.0);
  if (angle > 0 && angle < 90) {
    dx = -((int)(256 / t));
    dy = 1;
    dr_prediction_z1(dst, stride, bs, above, left, dx, dy, filter_type);
  } else if (angle > 90 && angle < 180) {
    t = -t;
    dx = (int)(256 / t);
    dy = (int)(256 * t);
    dr_prediction_z2(dst, stride, bs, above, left, dx, dy, filter_type);
  } else if (angle > 180 && angle < 270) {
    dx = 1;
    dy = -((int)(256 * t));
    dr_prediction_z3(dst, stride, bs, above, left, dx, dy, filter_type);
  } else if (angle == 90) {
    pred[V_PRED][tx_size](dst, stride, above, left);
  } else if (angle == 180) {
    pred[H_PRED][tx_size](dst, stride, above, left);
  }
}

static int filter_intra_taps_4[TX_SIZES][INTRA_MODES][4] = {
    {
        {735, 881, -537, -54},
        {1005, 519, -488, -11},
        {383, 990, -343, -6},
        {442, 805, -542, 319},
        {658, 616, -133, -116},
        {875, 442, -141, -151},
        {386, 741, -23, -80},
        {390, 1027, -446, 51},
        {679, 606, -523, 262},
        {903, 922, -778, -23},
    },
    {
        {648, 803, -444, 16},
        {972, 620, -576, 7},
        {561, 967, -499, -5},
        {585, 762, -468, 144},
        {596, 619, -182, -9},
        {895, 459, -176, -153},
        {557, 722, -126, -129},
        {601, 839, -523, 105},
        {562, 709, -499, 251},
        {803, 872, -695, 43},
    },
    {
        {423, 728, -347, 111},
        {963, 685, -665, 23},
        {281, 1024, -480, 216},
        {640, 596, -437, 78},
        {429, 669, -259, 99},
        {740, 646, -415, 23},
        {568, 771, -346, 40},
        {404, 833, -486, 209},
        {398, 712, -423, 307},
        {939, 935, -887, 17},
    },
    {
        {477, 737, -393, 150},
        {881, 630, -546, 67},
        {506, 984, -443, -20},
        {114, 459, -270, 528},
        {433, 528, 14, 3},
        {837, 470, -301, -30},
        {181, 777, 89, -107},
        {-29, 716, -232, 259},
        {589, 646, -495, 255},
        {740, 884, -728, 77},
    },
};

static void filter_intra_predictors_4tap(uint8_t *dst, ptrdiff_t stride, int bs,
                                         const uint8_t *above,
                                         const uint8_t *left,
                                         int mode) {
  int k, r, c;
  int pred[33][65];
  int mean, ipred;
  const TX_SIZE tx_size = (bs == 32) ? TX_32X32 :
      ((bs == 16) ? TX_16X16 : ((bs == 8) ? TX_8X8 : (TX_4X4)));
  const int c0 = filter_intra_taps_4[tx_size][mode][0];
  const int c1 = filter_intra_taps_4[tx_size][mode][1];
  const int c2 = filter_intra_taps_4[tx_size][mode][2];
  const int c3 = filter_intra_taps_4[tx_size][mode][3];

  k = 0;
  mean = 0;
  while (k < bs) {
    mean = mean + (int)left[k];
    mean = mean + (int)above[k];
    k++;
  }
  mean = (mean + bs) / (2 * bs);

  for (r = 0; r < bs; ++r)
    pred[r + 1][0] = (int)left[r] - mean;

  for (c = 0; c < 2 * bs + 1; ++c)
    pred[0][c] = (int)above[c - 1] - mean;

  for (r = 1; r < bs + 1; ++r)
    for (c = 1; c < 2 * bs + 1 - r; ++c) {
      ipred = c0 * pred[r - 1][c] + c1 * pred[r][c - 1] +
          c2 * pred[r - 1][c - 1] + c3 * pred[r - 1][c + 1];
      pred[r][c] = ipred < 0 ?
          -((-ipred + FILTER_INTRA_ROUND_VAL) >> FILTER_INTRA_PREC_BITS) :
          ((ipred + FILTER_INTRA_ROUND_VAL) >> FILTER_INTRA_PREC_BITS);
    }

  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      ipred = pred[r + 1][c + 1] + mean;
      dst[c] = clip_pixel(ipred);
    }
    dst += stride;
  }
}

static void dc_filter_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                               const uint8_t *above, const uint8_t *left) {
  filter_intra_predictors_4tap(dst, stride, bs, above, left, DC_PRED);
}

static void v_filter_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                               const uint8_t *above, const uint8_t *left) {
  filter_intra_predictors_4tap(dst, stride, bs, above, left, V_PRED);
}

static void h_filter_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                               const uint8_t *above, const uint8_t *left) {
  filter_intra_predictors_4tap(dst, stride, bs, above, left, H_PRED);
}

static void d45_filter_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                 const uint8_t *above, const uint8_t *left) {
  filter_intra_predictors_4tap(dst, stride, bs, above, left, D45_PRED);
}

static void d135_filter_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                  const uint8_t *above, const uint8_t *left) {
  filter_intra_predictors_4tap(dst, stride, bs, above, left, D135_PRED);
}

static void d117_filter_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                  const uint8_t *above, const uint8_t *left) {
  filter_intra_predictors_4tap(dst, stride, bs, above, left, D117_PRED);
}

static void d153_filter_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                  const uint8_t *above, const uint8_t *left) {
  filter_intra_predictors_4tap(dst, stride, bs, above, left, D153_PRED);
}

static void d207_filter_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                  const uint8_t *above, const uint8_t *left) {
  filter_intra_predictors_4tap(dst, stride, bs, above, left, D207_PRED);
}

static void d63_filter_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                 const uint8_t *above, const uint8_t *left) {
  filter_intra_predictors_4tap(dst, stride, bs, above, left, D63_PRED);
}

static void tm_filter_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                const uint8_t *above, const uint8_t *left) {
  filter_intra_predictors_4tap(dst, stride, bs, above, left, TM_PRED);
}

static void (*filter_intra_predictors[EXT_INTRA_MODES])(uint8_t *dst,
    ptrdiff_t stride, int bs, const uint8_t *above, const uint8_t *left) = {
        dc_filter_predictor, v_filter_predictor, h_filter_predictor,
        d45_filter_predictor, d135_filter_predictor, d117_filter_predictor,
        d153_filter_predictor, d207_filter_predictor, d63_filter_predictor,
        tm_filter_predictor,
};

#if CONFIG_VP9_HIGHBITDEPTH
static int highbd_intra_subpel_interp(int base, int shift, const uint16_t *ref,
                                      int ref_start_idx, int ref_end_idx,
                                      INTRA_FILTER filter_type) {
  int val, k, idx, filter_idx = 0;
  const int16_t *filter = NULL;

  if (filter_type == INTRA_FILTER_LINEAR) {
    val = ref[base] * (256 - shift) + ref[base + 1] * shift;
    val = ROUND_POWER_OF_TWO(val, 8);
  } else {
    filter_idx = ROUND_POWER_OF_TWO(shift, 8 - SUBPEL_BITS);
    filter = vp10_intra_filter_kernels[filter_type][filter_idx];

    if (filter_idx < (1 << SUBPEL_BITS)) {
      val = 0;
      for (k = 0; k < SUBPEL_TAPS; ++k) {
        idx = base + 1 - (SUBPEL_TAPS / 2) + k;
        idx = VPXMAX(VPXMIN(idx, ref_end_idx), ref_start_idx);
        val += ref[idx] * filter[k];
      }
      val = ROUND_POWER_OF_TWO(val, FILTER_BITS);
    } else {
      val = ref[base + 1];
    }
  }

  return val;
}

// Directional prediction, zone 1: 0 < angle < 90
static void highbd_dr_prediction_z1(uint16_t *dst, ptrdiff_t stride, int bs,
                                    const uint16_t *above, const uint16_t *left,
                                    int dx, int dy, int bd,
                                    INTRA_FILTER filter_type) {
  int r, c, x, y, base, shift, val;

  (void)left;
  (void)dy;
  assert(dy == 1);
  assert(dx < 0);

  for (r = 0; r < bs; ++r) {
    y = r + 1;
    for (c = 0; c < bs; ++c) {
      x = (c << 8) - y * dx;
      base = x >> 8;
      shift = x - (base << 8);
      if (base < 2 * bs - 1) {
        val = highbd_intra_subpel_interp(base, shift, above, 0, 2 * bs - 1,
                                         filter_type);
        dst[c] = clip_pixel_highbd(val, bd);
      } else {
        dst[c] = above[2 * bs - 1];
      }
    }
    dst += stride;
  }
}

// Directional prediction, zone 2: 90 < angle < 180
static void highbd_dr_prediction_z2(uint16_t *dst, ptrdiff_t stride, int bs,
                                    const uint16_t *above, const uint16_t *left,
                                    int dx, int dy, int bd,
                                    INTRA_FILTER filter_type) {
  int r, c, x, y, shift, val, base;

  assert(dx > 0);
  assert(dy > 0);

  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      y = r + 1;
      x = (c << 8) - y * dx;
      base = x >> 8;
      if (base >= -1) {
        shift = x - (base << 8);
        val = highbd_intra_subpel_interp(base, shift, above, -1, bs - 1,
                                         filter_type);
      } else {
        x = c + 1;
        y = (r << 8) - x * dy;
        base = y >> 8;
        if (base >= 0) {
          shift = y - (base  << 8);
          val = highbd_intra_subpel_interp(base, shift, left, 0, bs - 1,
                                           filter_type);
        } else {
          val = left[0];
        }
      }
      dst[c] = clip_pixel_highbd(val, bd);
    }
    dst += stride;
  }
}

// Directional prediction, zone 3: 180 < angle < 270
static void highbd_dr_prediction_z3(uint16_t *dst, ptrdiff_t stride, int bs,
                                    const uint16_t *above, const uint16_t *left,
                                    int dx, int dy, int bd,
                                    INTRA_FILTER filter_type) {
  int r, c, x, y, base, shift, val;

  (void)above;
  (void)dx;
  assert(dx == 1);
  assert(dy < 0);

  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      x = c + 1;
      y = (r << 8) - x * dy;
      base = y >> 8;
      shift = y - (base << 8);
      if (base < 2 * bs - 1) {
        val = highbd_intra_subpel_interp(base, shift, left, 0, 2 * bs - 1,
                                         filter_type);
        dst[c] = clip_pixel_highbd(val, bd);
      } else {
        dst[c] = left[2 * bs - 1];
      }
    }
    dst += stride;
  }
}

static INLINE void highbd_v_predictor(uint16_t *dst, ptrdiff_t stride,
                                      int bs, const uint16_t *above,
                                      const uint16_t *left, int bd) {
  int r;
  (void) left;
  (void) bd;
  for (r = 0; r < bs; r++) {
    memcpy(dst, above, bs * sizeof(uint16_t));
    dst += stride;
  }
}

static INLINE void highbd_h_predictor(uint16_t *dst, ptrdiff_t stride,
                                      int bs, const uint16_t *above,
                                      const uint16_t *left, int bd) {
  int r;
  (void) above;
  (void) bd;
  for (r = 0; r < bs; r++) {
    vpx_memset16(dst, left[r], bs);
    dst += stride;
  }
}

static void highbd_dr_predictor(uint16_t *dst, ptrdiff_t stride, int bs,
                                const uint16_t *above, const uint16_t *left,
                                int angle, int bd, INTRA_FILTER filter) {
  double t = 0;
  int dx, dy;

  if (angle != 90 && angle != 180)
    t = tan(angle * PI / 180.0);
  if (angle > 0 && angle < 90) {
    dx = -((int)(256 / t));
    dy = 1;
    highbd_dr_prediction_z1(dst, stride, bs, above, left, dx, dy, bd, filter);
  } else if (angle > 90 && angle < 180) {
    t = -t;
    dx = (int)(256 / t);
    dy = (int)(256 * t);
    highbd_dr_prediction_z2(dst, stride, bs, above, left, dx, dy, bd, filter);
  } else if (angle > 180 && angle < 270) {
    dx = 1;
    dy = -((int)(256 * t));
    highbd_dr_prediction_z3(dst, stride, bs, above, left, dx, dy, bd, filter);
  } else if (angle == 90) {
    highbd_v_predictor(dst, stride, bs, above, left, bd);
  } else if (angle == 180) {
    highbd_h_predictor(dst, stride, bs, above, left, bd);
  }
}

static void highbd_filter_intra_predictors_4tap(uint16_t *dst, ptrdiff_t stride,
                                                int bs, const uint16_t *above,
                                                const uint16_t *left, int mode,
                                                int bd) {
  int k, r, c;
  int pred[33][65];
  int mean, ipred;
  const TX_SIZE tx_size = (bs == 32) ? TX_32X32 :
      ((bs == 16) ? TX_16X16 : ((bs == 8) ? TX_8X8 : (TX_4X4)));
  const int c0 = filter_intra_taps_4[tx_size][mode][0];
  const int c1 = filter_intra_taps_4[tx_size][mode][1];
  const int c2 = filter_intra_taps_4[tx_size][mode][2];
  const int c3 = filter_intra_taps_4[tx_size][mode][3];

  k = 0;
  mean = 0;
  while (k < bs) {
    mean = mean + (int)left[k];
    mean = mean + (int)above[k];
    k++;
  }
  mean = (mean + bs) / (2 * bs);

  for (r = 0; r < bs; ++r)
    pred[r + 1][0] = (int)left[r] - mean;

  for (c = 0; c < 2 * bs + 1; ++c)
    pred[0][c] = (int)above[c - 1] - mean;

  for (r = 1; r < bs + 1; ++r)
    for (c = 1; c < 2 * bs + 1 - r; ++c) {
      ipred = c0 * pred[r - 1][c] + c1 * pred[r][c - 1] +
          c2 * pred[r - 1][c - 1] + c3 * pred[r - 1][c + 1];
      pred[r][c] = ipred < 0 ?
          -((-ipred + FILTER_INTRA_ROUND_VAL) >> FILTER_INTRA_PREC_BITS) :
          ((ipred + FILTER_INTRA_ROUND_VAL) >> FILTER_INTRA_PREC_BITS);
    }

  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      ipred = pred[r + 1][c + 1] + mean;
      dst[c] = clip_pixel_highbd(ipred, bd);
    }
    dst += stride;
  }
}

static void highbd_dc_filter_predictor(uint16_t *dst, ptrdiff_t stride,
                                       int bs, const uint16_t *above,
                                       const uint16_t *left, int bd) {
  highbd_filter_intra_predictors_4tap(dst, stride, bs, above, left, DC_PRED,
                                      bd);
}

static void highbd_v_filter_predictor(uint16_t *dst, ptrdiff_t stride,
                                      int bs, const uint16_t *above,
                                      const uint16_t *left, int bd) {
  highbd_filter_intra_predictors_4tap(dst, stride, bs, above, left, V_PRED,
                                      bd);
}

static void highbd_h_filter_predictor(uint16_t *dst, ptrdiff_t stride,
                                      int bs, const uint16_t *above,
                                      const uint16_t *left, int bd) {
  highbd_filter_intra_predictors_4tap(dst, stride, bs, above, left, H_PRED,
                                      bd);
}

static void highbd_d45_filter_predictor(uint16_t *dst, ptrdiff_t stride,
                                        int bs, const uint16_t *above,
                                        const uint16_t *left, int bd) {
  highbd_filter_intra_predictors_4tap(dst, stride, bs, above, left, D45_PRED,
                                      bd);
}

static void highbd_d135_filter_predictor(uint16_t *dst, ptrdiff_t stride,
                                         int bs, const uint16_t *above,
                                         const uint16_t *left, int bd) {
  highbd_filter_intra_predictors_4tap(dst, stride, bs, above, left, D135_PRED,
                                      bd);
}

static void highbd_d117_filter_predictor(uint16_t *dst, ptrdiff_t stride,
                                         int bs, const uint16_t *above,
                                         const uint16_t *left, int bd) {
  highbd_filter_intra_predictors_4tap(dst, stride, bs, above, left, D117_PRED,
                                      bd);
}

static void highbd_d153_filter_predictor(uint16_t *dst, ptrdiff_t stride,
                                         int bs, const uint16_t *above,
                                         const uint16_t *left, int bd) {
  highbd_filter_intra_predictors_4tap(dst, stride, bs, above, left, D153_PRED,
                                      bd);
}

static void highbd_d207_filter_predictor(uint16_t *dst, ptrdiff_t stride,
                                         int bs, const uint16_t *above,
                                         const uint16_t *left, int bd) {
  highbd_filter_intra_predictors_4tap(dst, stride, bs, above, left, D207_PRED,
                                      bd);
}

static void highbd_d63_filter_predictor(uint16_t *dst, ptrdiff_t stride,
                                        int bs, const uint16_t *above,
                                        const uint16_t *left, int bd) {
  highbd_filter_intra_predictors_4tap(dst, stride, bs, above, left, D63_PRED,
                                      bd);
}

static void highbd_tm_filter_predictor(uint16_t *dst, ptrdiff_t stride,
                                       int bs, const uint16_t *above,
                                       const uint16_t *left, int bd) {
  highbd_filter_intra_predictors_4tap(dst, stride, bs, above, left, TM_PRED,
                                      bd);
}

static void (*highbd_filter_intra_predictors[EXT_INTRA_MODES])(uint16_t *dst,
    ptrdiff_t stride, int bs, const uint16_t *above, const uint16_t *left,
    int bd) = {
        highbd_dc_filter_predictor, highbd_v_filter_predictor,
        highbd_h_filter_predictor, highbd_d45_filter_predictor,
        highbd_d135_filter_predictor, highbd_d117_filter_predictor,
        highbd_d153_filter_predictor, highbd_d207_filter_predictor,
        highbd_d63_filter_predictor, highbd_tm_filter_predictor,
};
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_EXT_INTRA

#if CONFIG_VP9_HIGHBITDEPTH
static void build_intra_predictors_high(const MACROBLOCKD *xd,
                                        const uint8_t *ref8,
                                        int ref_stride,
                                        uint8_t *dst8,
                                        int dst_stride,
                                        PREDICTION_MODE mode,
                                        TX_SIZE tx_size,
                                        int n_top_px, int n_topright_px,
                                        int n_left_px, int n_bottomleft_px,
                                        int plane) {
  int i;
  uint16_t *dst = CONVERT_TO_SHORTPTR(dst8);
  uint16_t *ref = CONVERT_TO_SHORTPTR(ref8);
  DECLARE_ALIGNED(16, uint16_t, left_col[64]);
  DECLARE_ALIGNED(16, uint16_t, above_data[64 + 16]);
  uint16_t *above_row = above_data + 16;
  const uint16_t *const_above_row = above_row;
  const int bs = 4 << tx_size;
  int need_left = extend_modes[mode] & NEED_LEFT;
  int need_above = extend_modes[mode] & NEED_ABOVE;
  const uint16_t *above_ref = ref - ref_stride;
  int base = 128 << (xd->bd - 8);
  // 127 127 127 .. 127 127 127 127 127 127
  // 129  A   B  ..  Y   Z
  // 129  C   D  ..  W   X
  // 129  E   F  ..  U   V
  // 129  G   H  ..  S   T   T   T   T   T

#if CONFIG_EXT_INTRA
  const EXT_INTRA_MODE_INFO *ext_intra_mode_info =
      &xd->mi[0]->mbmi.ext_intra_mode_info;
  const EXT_INTRA_MODE ext_intra_mode =
      ext_intra_mode_info->ext_intra_mode[plane != 0];
  int p_angle = 0;

  if (mode != DC_PRED && mode != TM_PRED &&
      xd->mi[0]->mbmi.sb_type >= BLOCK_8X8) {
    p_angle = mode_to_angle_map[mode] +
        xd->mi[0]->mbmi.angle_delta[plane != 0] * ANGLE_STEP;
    if (p_angle <= 90)
      need_above = 1, need_left = 0;
    else if (p_angle < 180)
      need_above = 1, need_left = 1;
    else
      need_above = 0, need_left = 1;
  }

  if (ext_intra_mode_info->use_ext_intra_mode[plane != 0]) {
    EXT_INTRA_MODE ext_intra_mode =
        ext_intra_mode_info->ext_intra_mode[plane != 0];
    need_left = ext_intra_extend_modes[ext_intra_mode] & NEED_LEFT;
    need_above = ext_intra_extend_modes[ext_intra_mode] & NEED_ABOVE;
  }
#endif  // CONFIG_EXT_INTRA

  (void) plane;

  // NEED_LEFT
  if (need_left) {
#if CONFIG_EXT_INTRA
    int need_bottom;
    if (ext_intra_mode_info->use_ext_intra_mode[plane != 0]) {
        need_bottom = 0;
    } else if (mode != DC_PRED && mode != TM_PRED &&
        xd->mi[0]->mbmi.sb_type >= BLOCK_8X8) {
        need_bottom = p_angle > 180;
    } else {
      need_bottom = !!(extend_modes[mode] & NEED_BOTTOMLEFT);
    }
#else
    const int need_bottom = !!(extend_modes[mode] & NEED_BOTTOMLEFT);
#endif  // CONFIG_EXT_INTRA
    i = 0;
    if (n_left_px > 0) {
      for (; i < n_left_px; i++)
        left_col[i] = ref[i * ref_stride - 1];
      if (need_bottom && n_bottomleft_px > 0) {
        assert(i == bs);
        for (; i < bs + n_bottomleft_px; i++)
          left_col[i] = ref[i * ref_stride - 1];
      }
      if (i < (bs << need_bottom))
        vpx_memset16(&left_col[i], left_col[i - 1], (bs << need_bottom) - i);
    } else {
      vpx_memset16(left_col, base + 1, bs << need_bottom);
    }
  }

  // NEED_ABOVE
  if (need_above) {
#if CONFIG_EXT_INTRA
    int need_right;
    if (ext_intra_mode_info->use_ext_intra_mode[plane != 0]) {
      need_right = 1;
    } else if (mode != DC_PRED && mode != TM_PRED &&
        xd->mi[0]->mbmi.sb_type >= BLOCK_8X8) {
      need_right = p_angle < 90;
    } else {
      need_right = !!(extend_modes[mode] & NEED_ABOVERIGHT);
    }
#else
    const int need_right = !!(extend_modes[mode] & NEED_ABOVERIGHT);
#endif  // CONFIG_EXT_INTRA
    if (n_top_px > 0) {
      memcpy(above_row, above_ref, n_top_px * 2);
      i = n_top_px;
      if (need_right && n_topright_px > 0) {
        assert(n_top_px == bs);
        memcpy(above_row + bs, above_ref + bs, n_topright_px * 2);
        i += n_topright_px;
      }
      if (i < (bs << need_right))
        vpx_memset16(&above_row[i], above_row[i - 1], (bs << need_right) - i);
    } else {
      vpx_memset16(above_row, base - 1, bs << need_right);
    }
  }

#if CONFIG_EXT_INTRA
  if (ext_intra_mode_info->use_ext_intra_mode[plane != 0] ||
      (extend_modes[mode] & NEED_ABOVELEFT) ||
      (mode != DC_PRED && mode != TM_PRED &&
        xd->mi[0]->mbmi.sb_type >= BLOCK_8X8)) {
    above_row[-1] = n_top_px > 0 ?
        (n_left_px > 0 ? above_ref[-1] : base + 1) : base - 1;
  }
#else
  if ((extend_modes[mode] & NEED_ABOVELEFT)) {
    above_row[-1] = n_top_px > 0 ?
        (n_left_px > 0 ? above_ref[-1] : base + 1) : base - 1;
  }
#endif  // CONFIG_EXT_INTRA

#if CONFIG_EXT_INTRA
  if (ext_intra_mode_info->use_ext_intra_mode[plane != 0]) {
    highbd_filter_intra_predictors[ext_intra_mode](dst, dst_stride, bs,
        const_above_row, left_col, xd->bd);
    return;
  }

  if (mode != DC_PRED && mode != TM_PRED &&
      xd->mi[0]->mbmi.sb_type >= BLOCK_8X8) {
    INTRA_FILTER filter = INTRA_FILTER_LINEAR;
    if (plane == 0 && pick_intra_filter(p_angle))
      filter = xd->mi[0]->mbmi.intra_filter;
    highbd_dr_predictor(dst, dst_stride, bs, const_above_row, left_col,
                        p_angle, xd->bd, filter);
    return;
  }
#endif  // CONFIG_EXT_INTRA

  // predict
  if (mode == DC_PRED) {
    dc_pred_high[n_left_px > 0][n_top_px > 0][tx_size](dst, dst_stride,
                                                       const_above_row,
                                                       left_col, xd->bd);
  } else {
    pred_high[mode][tx_size](dst, dst_stride, const_above_row, left_col,
                             xd->bd);
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

static void build_intra_predictors(const MACROBLOCKD *xd, const uint8_t *ref,
                                   int ref_stride, uint8_t *dst, int dst_stride,
                                   PREDICTION_MODE mode, TX_SIZE tx_size,
                                   int n_top_px, int n_topright_px,
                                   int n_left_px, int n_bottomleft_px,
                                   int plane) {
  int i;
  DECLARE_ALIGNED(16, uint8_t, left_col[64]);
  const uint8_t *above_ref = ref - ref_stride;
  DECLARE_ALIGNED(16, uint8_t, above_data[64 + 16]);
  uint8_t *above_row = above_data + 16;
  const uint8_t *const_above_row = above_row;
  const int bs = 4 << tx_size;
  int need_left = extend_modes[mode] & NEED_LEFT;
  int need_above = extend_modes[mode] & NEED_ABOVE;
#if CONFIG_EXT_INTRA
  const EXT_INTRA_MODE_INFO *ext_intra_mode_info =
      &xd->mi[0]->mbmi.ext_intra_mode_info;
  const EXT_INTRA_MODE ext_intra_mode =
      ext_intra_mode_info->ext_intra_mode[plane != 0];
  int p_angle = 0;

  if (mode != DC_PRED && mode != TM_PRED &&
      xd->mi[0]->mbmi.sb_type >= BLOCK_8X8) {
    p_angle = mode_to_angle_map[mode] +
        xd->mi[0]->mbmi.angle_delta[plane != 0] * ANGLE_STEP;
    if (p_angle <= 90)
      need_above = 1, need_left = 0;
    else if (p_angle < 180)
      need_above = 1, need_left = 1;
    else
      need_above = 0, need_left = 1;
  }

  if (ext_intra_mode_info->use_ext_intra_mode[plane != 0]) {
    EXT_INTRA_MODE ext_intra_mode =
        ext_intra_mode_info->ext_intra_mode[plane != 0];
    need_left = ext_intra_extend_modes[ext_intra_mode] & NEED_LEFT;
    need_above = ext_intra_extend_modes[ext_intra_mode] & NEED_ABOVE;
  }
#endif  // CONFIG_EXT_INTRA

  // 127 127 127 .. 127 127 127 127 127 127
  // 129  A   B  ..  Y   Z
  // 129  C   D  ..  W   X
  // 129  E   F  ..  U   V
  // 129  G   H  ..  S   T   T   T   T   T
  // ..

  (void) xd;
  (void) plane;
  assert(n_top_px >= 0);
  assert(n_topright_px >= 0);
  assert(n_left_px >= 0);
  assert(n_bottomleft_px >= 0);

  // NEED_LEFT
  if (need_left) {
#if CONFIG_EXT_INTRA
    int need_bottom;
    if (ext_intra_mode_info->use_ext_intra_mode[plane != 0]) {
      need_bottom = 0;
    } else if (mode != DC_PRED && mode != TM_PRED &&
        xd->mi[0]->mbmi.sb_type >= BLOCK_8X8) {
      need_bottom = p_angle > 180;
    } else {
      need_bottom = !!(extend_modes[mode] & NEED_BOTTOMLEFT);
    }
#else
    const int need_bottom = !!(extend_modes[mode] & NEED_BOTTOMLEFT);
#endif  // CONFIG_EXT_INTRA
    i = 0;
    if (n_left_px > 0) {
      for (; i < n_left_px; i++)
        left_col[i] = ref[i * ref_stride - 1];
      if (need_bottom && n_bottomleft_px > 0) {
        assert(i == bs);
        for (; i < bs + n_bottomleft_px; i++)
          left_col[i] = ref[i * ref_stride - 1];
      }
      if (i < (bs << need_bottom))
        memset(&left_col[i], left_col[i - 1], (bs << need_bottom) - i);
    } else {
      memset(left_col, 129, bs << need_bottom);
    }
  }

  // NEED_ABOVE
  if (need_above) {
#if CONFIG_EXT_INTRA
    int need_right;
    if (ext_intra_mode_info->use_ext_intra_mode[plane != 0]) {
      need_right = 1;
    } else if (mode != DC_PRED && mode != TM_PRED &&
        xd->mi[0]->mbmi.sb_type >= BLOCK_8X8) {
      need_right = p_angle < 90;
    } else {
      need_right = !!(extend_modes[mode] & NEED_ABOVERIGHT);
    }
#else
    const int need_right = !!(extend_modes[mode] & NEED_ABOVERIGHT);
#endif  // CONFIG_EXT_INTRA
    if (n_top_px > 0) {
      memcpy(above_row, above_ref, n_top_px);
      i = n_top_px;
      if (need_right && n_topright_px > 0) {
        assert(n_top_px == bs);
        memcpy(above_row + bs, above_ref + bs, n_topright_px);
        i += n_topright_px;
      }
      if (i < (bs << need_right))
        memset(&above_row[i], above_row[i - 1], (bs << need_right) - i);
    } else {
      memset(above_row, 127, bs << need_right);
    }
  }

#if CONFIG_EXT_INTRA
  if (ext_intra_mode_info->use_ext_intra_mode[plane != 0] ||
      (extend_modes[mode] & NEED_ABOVELEFT) ||
      (mode != DC_PRED && mode != TM_PRED &&
          xd->mi[0]->mbmi.sb_type >= BLOCK_8X8)) {
    above_row[-1] = n_top_px > 0 ? (n_left_px > 0 ? above_ref[-1] : 129) : 127;
  }
#else
  if ((extend_modes[mode] & NEED_ABOVELEFT)) {
    above_row[-1] = n_top_px > 0 ? (n_left_px > 0 ? above_ref[-1] : 129) : 127;
  }
#endif  // CONFIG_EXT_INTRA

#if CONFIG_EXT_INTRA
  if (ext_intra_mode_info->use_ext_intra_mode[plane != 0]) {
    filter_intra_predictors[ext_intra_mode](dst, dst_stride, bs,
        const_above_row, left_col);
    return;
  }

  if (mode != DC_PRED && mode != TM_PRED &&
      xd->mi[0]->mbmi.sb_type >= BLOCK_8X8) {
    INTRA_FILTER filter = INTRA_FILTER_LINEAR;
    if (plane == 0 && pick_intra_filter(p_angle))
      filter = xd->mi[0]->mbmi.intra_filter;
    dr_predictor(dst, dst_stride, tx_size, const_above_row, left_col, p_angle,
                 filter);
    return;
  }
#endif  // CONFIG_EXT_INTRA

  // predict
  if (mode == DC_PRED) {
    dc_pred[n_left_px > 0][n_top_px > 0][tx_size](dst, dst_stride,
                                                  const_above_row, left_col);
  } else {
    pred[mode][tx_size](dst, dst_stride, const_above_row, left_col);
  }
}

void vp10_predict_intra_block(const MACROBLOCKD *xd, int bwl_in, int bhl_in,
                              TX_SIZE tx_size, PREDICTION_MODE mode,
                              const uint8_t *ref, int ref_stride,
                              uint8_t *dst, int dst_stride,
                              int col_off, int row_off, int plane) {
  const int txw = (1 << tx_size);
  const int have_top = row_off || xd->up_available;
  const int have_left = col_off || xd->left_available;
  const int x = col_off * 4;
  const int y = row_off * 4;
  const int bw = VPXMAX(2, 1 << bwl_in);
  const int bh = VPXMAX(2, 1 << bhl_in);
  const int mi_row = -xd->mb_to_top_edge >> (3 + MI_SIZE_LOG2);
  const int mi_col = -xd->mb_to_left_edge >> (3 + MI_SIZE_LOG2);
  const BLOCK_SIZE bsize = xd->mi[0]->mbmi.sb_type;
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const int right_available =
      mi_col + (bw >> !pd->subsampling_x) < xd->tile.mi_col_end;
  const int have_right = vp10_has_right(bsize, mi_row, mi_col,
                                        right_available,
                                        tx_size, row_off, col_off,
                                        pd->subsampling_x);
  const int have_bottom = vp10_has_bottom(bsize, mi_row, mi_col,
                                          xd->mb_to_bottom_edge > 0,
                                          tx_size, row_off, col_off,
                                          pd->subsampling_y);
  const int wpx = 4 * bw;
  const int hpx = 4 * bh;
  const int txpx = 4 * txw;
  // Distance between the right edge of this prediction block to
  // the frame right edge
  const int xr = (xd->mb_to_right_edge >> (3 + pd->subsampling_x)) +
      (wpx - x - txpx);
  // Distance between the bottom edge of this prediction block to
  // the frame bottom edge
  const int yd = (xd->mb_to_bottom_edge >> (3 + pd->subsampling_y)) +
      (hpx - y - txpx);

  if (xd->mi[0]->mbmi.palette_mode_info.palette_size[plane != 0] > 0) {
    const int bs = 4 * (1 << tx_size);
    const int stride = 4 * (1 << bwl_in);
    int r, c;
    uint8_t *map = NULL;
#if CONFIG_VP9_HIGHBITDEPTH
    uint16_t *palette = xd->mi[0]->mbmi.palette_mode_info.palette_colors +
        plane * PALETTE_MAX_SIZE;
#else
    uint8_t *palette = xd->mi[0]->mbmi.palette_mode_info.palette_colors +
        plane * PALETTE_MAX_SIZE;
#endif  // CONFIG_VP9_HIGHBITDEPTH

    map = xd->plane[plane != 0].color_index_map;

#if CONFIG_VP9_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
      uint16_t *dst16 = CONVERT_TO_SHORTPTR(dst);
      for (r = 0; r < bs; ++r)
        for (c = 0; c < bs; ++c)
          dst16[r * dst_stride + c] =
              palette[map[(r + y) * stride + c + x]];
    } else {
      for (r = 0; r < bs; ++r)
        for (c = 0; c < bs; ++c)
          dst[r * dst_stride + c] =
              (uint8_t)(palette[map[(r + y) * stride + c + x]]);
    }
#else
    for (r = 0; r < bs; ++r)
      for (c = 0; c < bs; ++c)
        dst[r * dst_stride + c] = palette[map[(r + y) * stride + c + x]];
#endif  // CONFIG_VP9_HIGHBITDEPTH
    return;
  }

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    build_intra_predictors_high(xd, ref, ref_stride, dst, dst_stride, mode,
                                tx_size,
                                have_top ? VPXMIN(txpx, xr + txpx) : 0,
                                have_top && have_right ? VPXMIN(txpx, xr) : 0,
                                have_left ? VPXMIN(txpx, yd + txpx) : 0,
                                have_bottom && have_left ? VPXMIN(txpx, yd) : 0,
                                plane);
    return;
  }
#endif
  build_intra_predictors(xd, ref, ref_stride, dst, dst_stride, mode,
                         tx_size,
                         have_top ? VPXMIN(txpx, xr + txpx) : 0,
                         have_top && have_right ? VPXMIN(txpx, xr) : 0,
                         have_left ? VPXMIN(txpx, yd + txpx) : 0,
                         have_bottom && have_left ? VPXMIN(txpx, yd) : 0,
                         plane);
}

void vp10_init_intra_predictors(void) {
  once(vp10_init_intra_predictors_internal);
}
