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
#include "vpx_ports/system_state.h"

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

static const uint8_t orders_128x128[1] = { 0 };
static const uint8_t orders_128x64[2] = { 0, 1 };
static const uint8_t orders_64x128[2] = { 0, 1 };
static const uint8_t orders_64x64[4] = {
  0, 1,
  2, 3,
};
static const uint8_t orders_64x32[8] = {
  0, 2,
  1, 3,
  4, 6,
  5, 7,
};
static const uint8_t orders_32x64[8] = {
  0, 1, 2, 3,
  4, 5, 6, 7,
};
static const uint8_t orders_32x32[16] = {
  0,   1,  4,  5,
  2,   3,  6,  7,
  8,   9, 12, 13,
  10, 11, 14, 15,
};
static const uint8_t orders_32x16[32] = {
  0,   2,  8, 10,
  1,   3,  9, 11,
  4,   6, 12, 14,
  5,   7, 13, 15,
  16, 18, 24, 26,
  17, 19, 25, 27,
  20, 22, 28, 30,
  21, 23, 29, 31,
};
static const uint8_t orders_16x32[32] = {
  0,   1,  2,  3,  8,  9, 10, 11,
  4,   5,  6,  7, 12, 13, 14, 15,
  16, 17, 18, 19, 24, 25, 26, 27,
  20, 21, 22, 23, 28, 29, 30, 31,
};
static const uint8_t orders_16x16[64] = {
  0,   1,  4,  5, 16, 17, 20, 21,
  2,   3,  6,  7, 18, 19, 22, 23,
  8,   9, 12, 13, 24, 25, 28, 29,
  10, 11, 14, 15, 26, 27, 30, 31,
  32, 33, 36, 37, 48, 49, 52, 53,
  34, 35, 38, 39, 50, 51, 54, 55,
  40, 41, 44, 45, 56, 57, 60, 61,
  42, 43, 46, 47, 58, 59, 62, 63,
};

#if CONFIG_EXT_PARTITION
static const uint8_t orders_16x8[128] = {
  0,   2,  8, 10,  32,  34,  40,  42,
  1,   3,  9, 11,  33,  35,  41,  43,
  4,   6, 12, 14,  36,  38,  44,  46,
  5,   7, 13, 15,  37,  39,  45,  47,
  16, 18, 24, 26,  48,  50,  56,  58,
  17, 19, 25, 27,  49,  51,  57,  59,
  20, 22, 28, 30,  52,  54,  60,  62,
  21, 23, 29, 31,  53,  55,  61,  63,
  64, 66, 72, 74,  96,  98, 104, 106,
  65, 67, 73, 75,  97,  99, 105, 107,
  68, 70, 76, 78, 100, 102, 108, 110,
  69, 71, 77, 79, 101, 103, 109, 111,
  80, 82, 88, 90, 112, 114, 120, 122,
  81, 83, 89, 91, 113, 115, 121, 123,
  84, 86, 92, 94, 116, 118, 124, 126,
  85, 87, 93, 95, 117, 119, 125, 127,
};
static const uint8_t orders_8x16[128] = {
  0,   1,  2,  3,  8,  9, 10, 11,  32,  33,  34,  35,  40,  41,  42,  43,
  4,   5,  6,  7, 12, 13, 14, 15,  36,  37,  38,  39,  44,  45,  46,  47,
  16, 17, 18, 19, 24, 25, 26, 27,  48,  49,  50,  51,  56,  57,  58,  59,
  20, 21, 22, 23, 28, 29, 30, 31,  52,  53,  54,  55,  60,  61,  62,  63,
  64, 65, 66, 67, 72, 73, 74, 75,  96,  97,  98,  99, 104, 105, 106, 107,
  68, 69, 70, 71, 76, 77, 78, 79, 100, 101, 102, 103, 108, 109, 110, 111,
  80, 81, 82, 83, 88, 89, 90, 91, 112, 113, 114, 115, 120, 121, 122, 123,
  84, 85, 86, 87, 92, 93, 94, 95, 116, 117, 118, 119, 124, 125, 126, 127,
};
static const uint8_t orders_8x8[256] = {
0,     1,   4,   5,  16,  17,  20,  21,  64,  65,  68,  69,  80,  81,  84,  85,
2,     3,   6,   7,  18,  19,  22,  23,  66,  67,  70,  71,  82,  83,  86,  87,
8,     9,  12,  13,  24,  25,  28,  29,  72,  73,  76,  77,  88,  89,  92,  93,
10,   11,  14,  15,  26,  27,  30,  31,  74,  75,  78,  79,  90,  91,  94,  95,
32,   33,  36,  37,  48,  49,  52,  53,  96,  97, 100, 101, 112, 113, 116, 117,
34,   35,  38,  39,  50,  51,  54,  55,  98,  99, 102, 103, 114, 115, 118, 119,
40,   41,  44,  45,  56,  57,  60,  61, 104, 105, 108, 109, 120, 121, 124, 125,
42,   43,  46,  47,  58,  59,  62,  63, 106, 107, 110, 111, 122, 123, 126, 127,
128, 129, 132, 133, 144, 145, 148, 149, 192, 193, 196, 197, 208, 209, 212, 213,
130, 131, 134, 135, 146, 147, 150, 151, 194, 195, 198, 199, 210, 211, 214, 215,
136, 137, 140, 141, 152, 153, 156, 157, 200, 201, 204, 205, 216, 217, 220, 221,
138, 139, 142, 143, 154, 155, 158, 159, 202, 203, 206, 207, 218, 219, 222, 223,
160, 161, 164, 165, 176, 177, 180, 181, 224, 225, 228, 229, 240, 241, 244, 245,
162, 163, 166, 167, 178, 179, 182, 183, 226, 227, 230, 231, 242, 243, 246, 247,
168, 169, 172, 173, 184, 185, 188, 189, 232, 233, 236, 237, 248, 249, 252, 253,
170, 171, 174, 175, 186, 187, 190, 191, 234, 235, 238, 239, 250, 251, 254, 255,
};

static const uint8_t *const orders[BLOCK_SIZES] = {
  //                              4X4
                                  orders_8x8,
  // 4X8,         8X4,            8X8
  orders_8x8,     orders_8x8,     orders_8x8,
  // 8X16,        16X8,           16X16
  orders_8x16,    orders_16x8,    orders_16x16,
  // 16X32,       32X16,          32X32
  orders_16x32,   orders_32x16,   orders_32x32,
  // 32X64,       64X32,          64X64
  orders_32x64,   orders_64x32,   orders_64x64,
  // 64x128,      128x64,         128x128
  orders_64x128,  orders_128x64,  orders_128x128
};
#else
static const uint8_t *const orders[BLOCK_SIZES] = {
  //                              4X4
                                  orders_16x16,
  // 4X8,         8X4,            8X8
  orders_16x16,   orders_16x16,   orders_16x16,
  // 8X16,        16X8,           16X16
  orders_16x32,   orders_32x16,   orders_32x32,
  // 16X32,       32X16,          32X32
  orders_32x64,   orders_64x32,   orders_64x64,
  // 32X64,       64X32,          64X64
  orders_64x128,  orders_128x64,  orders_128x128
};
#endif  // CONFIG_EXT_PARTITION

#if CONFIG_EXT_PARTITION_TYPES
static const uint8_t orders_verta_64x64[4] = {
  0, 2,
  1, 2,
};
static const uint8_t orders_verta_32x32[16] = {
  0,   2,  4,  6,
  1,   2,  5,  6,
  8,  10, 12, 14,
  9,  10, 13, 14,
};
static const uint8_t orders_verta_16x16[64] = {
  0,   2,  4,  6, 16, 18, 20, 22,
  1,   2,  5,  6, 17, 18, 21, 22,
  8,  10, 12, 14, 24, 26, 28, 30,
  9,  10, 13, 14, 25, 26, 29, 30,
  32, 34, 36, 38, 48, 50, 52, 54,
  33, 34, 37, 38, 49, 50, 53, 54,
  40, 42, 44, 46, 56, 58, 60, 62,
  41, 42, 45, 46, 57, 58, 61, 62,
};
#if CONFIG_EXT_PARTITION
static const uint8_t orders_verta_8x8[256] = {
0,     2,   4,   6,  16,  18,  20,  22,  64,  66,  68,  70,  80,  82,  84,  86,
1,     2,   5,   6,  17,  18,  21,  22,  65,  66,  69,  70,  81,  82,  85,  86,
8,    10,  12,  14,  24,  26,  28,  30,  72,  74,  76,  78,  88,  90,  92,  94,
9,    10,  13,  14,  25,  26,  29,  30,  73,  74,  77,  78,  89,  90,  93,  94,
32,   34,  36,  38,  48,  50,  52,  54,  96,  98, 100, 102, 112, 114, 116, 118,
33,   34,  37,  38,  49,  50,  53,  54,  97,  98, 101, 102, 113, 114, 117, 118,
40,   42,  44,  46,  56,  58,  60,  62, 104, 106, 108, 110, 120, 122, 124, 126,
41,   42,  45,  46,  57,  58,  61,  62, 105, 106, 109, 110, 121, 122, 125, 126,
128, 130, 132, 134, 144, 146, 148, 150, 192, 194, 196, 198, 208, 210, 212, 214,
129, 130, 133, 134, 145, 146, 149, 150, 193, 194, 197, 198, 209, 210, 213, 214,
136, 138, 140, 142, 152, 154, 156, 158, 200, 202, 204, 206, 216, 218, 220, 222,
137, 138, 141, 142, 153, 154, 157, 158, 201, 202, 205, 206, 217, 218, 221, 222,
160, 162, 164, 166, 176, 178, 180, 182, 224, 226, 228, 230, 240, 242, 244, 246,
161, 162, 165, 166, 177, 178, 181, 182, 225, 226, 229, 230, 241, 242, 245, 246,
168, 170, 172, 174, 184, 186, 188, 190, 232, 234, 236, 238, 248, 250, 252, 254,
169, 170, 173, 174, 185, 186, 189, 190, 233, 234, 237, 238, 249, 250, 253, 254,
};
static const uint8_t *const orders_verta[BLOCK_SIZES] = {
  //                                  4X4
                                      orders_verta_8x8,
  // 4X8,           8X4,              8X8
  orders_verta_8x8, orders_verta_8x8, orders_verta_8x8,
  // 8X16,          16X8,             16X16
  orders_8x16,      orders_16x8,      orders_verta_16x16,
  // 16X32,         32X16,            32X32
  orders_16x32,     orders_32x16,     orders_verta_32x32,
  // 32X64,         64X32,            64X64
  orders_32x64,     orders_64x32,     orders_verta_64x64,
  // 64x128,        128x64,           128x128
  orders_64x128,    orders_128x64,    orders_128x128
};
#else
static const uint8_t *const orders_verta[BLOCK_SIZES] = {
  //                                      4X4
                                          orders_verta_16x16,
  // 4X8,             8X4,                8X8
  orders_verta_16x16, orders_verta_16x16, orders_verta_16x16,
  // 8X16,            16X8,               16X16
  orders_16x32,       orders_32x16,       orders_verta_32x32,
  // 16X32,           32X16,              32X32
  orders_32x64,       orders_64x32,       orders_verta_64x64,
  // 32X64,           64X32,              64X64
  orders_64x128,      orders_128x64,      orders_128x128
};
#endif  // CONFIG_EXT_PARTITION
#endif  // CONFIG_EXT_PARTITION_TYPES

static int vp10_has_right(BLOCK_SIZE bsize, int mi_row, int mi_col,
                          int right_available,
#if CONFIG_EXT_PARTITION_TYPES
                          PARTITION_TYPE partition,
#endif
                          TX_SIZE txsz, int y, int x, int ss_x) {
  const int wl = mi_width_log2_lookup[bsize];
  const int w = VPXMAX(num_4x4_blocks_wide_lookup[bsize] >> ss_x, 1);
  const int step = 1 << txsz;

  if (!right_available) {
    return 0;
  } else {
    // Handle block size 4x8 and 4x4
    if (ss_x == 0 && num_4x4_blocks_wide_lookup[bsize] < 2 && x == 0)
      return 1;

    if (y == 0) {
      const int hl = mi_height_log2_lookup[bsize];
      const uint8_t *order;
      int my_order, tr_order;
#if CONFIG_EXT_PARTITION_TYPES
      if (partition == PARTITION_VERT_A)
        order = orders_verta[bsize];
      else
#endif  // CONFIG_EXT_PARTITION_TYPES
      order = orders[bsize];

      if (x + step < w)
        return 1;

      mi_row = (mi_row & MAX_MIB_MASK) >> hl;
      mi_col = (mi_col & MAX_MIB_MASK) >> wl;

      // If top row of coding unit
      if (mi_row == 0)
        return 1;

      // If rightmost column of coding unit
      if (((mi_col + 1) << wl) >= MAX_MIB_SIZE)
        return 0;

      my_order = order[((mi_row + 0) << (MAX_MIB_SIZE_LOG2 - wl)) + mi_col + 0];
      tr_order = order[((mi_row - 1) << (MAX_MIB_SIZE_LOG2 - wl)) + mi_col + 1];

      return my_order > tr_order;
    } else {
      return x + step < w;
    }
  }
}

static int vp10_has_bottom(BLOCK_SIZE bsize, int mi_row, int mi_col,
                           int bottom_available, TX_SIZE txsz,
                           int y, int x, int ss_y) {
  if (!bottom_available || x != 0) {
    return 0;
  } else {
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

    mi_row = (mi_row & MAX_MIB_MASK) >> hl;
    mi_col = (mi_col & MAX_MIB_MASK) >> wl;

    if (mi_col == 0)
      return (mi_row << (hl + !ss_y)) + y + step < (MAX_MIB_SIZE << !ss_y);

    if (((mi_row + 1) << hl) >= MAX_MIB_SIZE)
      return 0;

    my_order = order[((mi_row + 0) << (MAX_MIB_SIZE_LOG2 - wl)) + mi_col + 0];
    bl_order = order[((mi_row + 1) << (MAX_MIB_SIZE_LOG2 - wl)) + mi_col - 1];

    return bl_order < my_order;
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
#define FILTER_INTRA_PREC_BITS 10

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
    DECLARE_ALIGNED(16, uint8_t, buf[SUBPEL_SHIFTS][MAX_SB_SIZE]);
    DECLARE_ALIGNED(16, uint8_t, src[MAX_SB_SIZE + SUBPEL_TAPS]);
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
          const int16_t *filter = vp10_intra_filter_kernels[filter_type][shift];
          vpx_convolve8_horiz(src + pad_size, 2 * bs, buf[shift], 2 * bs,
                              filter, 16,
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
    DECLARE_ALIGNED(16, uint8_t, buf[MAX_SB_SIZE][4 * SUBPEL_SHIFTS]);
    DECLARE_ALIGNED(16, uint8_t, src[(MAX_SB_SIZE + SUBPEL_TAPS) * 4]);
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
          const int16_t *filter = vp10_intra_filter_kernels[filter_type][shift];
          vpx_convolve8_vert(src + 4 * pad_size, 4,
                             buf[0] + 4 * shift, 4 * SUBPEL_SHIFTS, NULL, 16,
                             filter, 16,
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
  const int dx = (int)dr_intra_derivative[angle][0];
  const int dy = (int)dr_intra_derivative[angle][1];
  const int bs = 4 * num_4x4_blocks_wide_txsize_lookup[tx_size];
  assert(angle > 0 && angle < 270);

  if (angle > 0 && angle < 90) {
    dr_prediction_z1(dst, stride, bs, above, left, dx, dy, filter_type);
  } else if (angle > 90 && angle < 180) {
    dr_prediction_z2(dst, stride, bs, above, left, dx, dy, filter_type);
  } else if (angle > 180 && angle < 270) {
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
      pred[r][c] = ROUND_POWER_OF_TWO_SIGNED(ipred, FILTER_INTRA_PREC_BITS);
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
  const int dx = (int)dr_intra_derivative[angle][0];
  const int dy = (int)dr_intra_derivative[angle][1];
  assert(angle > 0 && angle < 270);

  if (angle > 0 && angle < 90) {
    highbd_dr_prediction_z1(dst, stride, bs, above, left, dx, dy, bd, filter);
  } else if (angle > 90 && angle < 180) {
    highbd_dr_prediction_z2(dst, stride, bs, above, left, dx, dy, bd, filter);
  } else if (angle > 180 && angle < 270) {
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
      pred[r][c] = ROUND_POWER_OF_TWO_SIGNED(ipred, FILTER_INTRA_PREC_BITS);
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
  DECLARE_ALIGNED(16, uint16_t, left_col[MAX_SB_SIZE]);
  DECLARE_ALIGNED(16, uint16_t, above_data[MAX_SB_SIZE + 16]);
  uint16_t *above_row = above_data + 16;
  const uint16_t *const_above_row = above_row;
  const int bs = 4 * num_4x4_blocks_wide_txsize_lookup[tx_size];
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
  assert(n_top_px >= 0);
  assert(n_topright_px >= 0);
  assert(n_left_px >= 0);
  assert(n_bottomleft_px >= 0);

  if ((!need_above && n_left_px == 0) || (!need_left && n_top_px == 0)) {
    int i;
    const int val = (n_left_px == 0) ? base + 1 : base - 1;
    for (i = 0; i < bs; ++i) {
      vpx_memset16(dst, val, bs);
      dst += dst_stride;
    }
    return;
  }

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
    if (plane == 0 && vp10_is_intra_filter_switchable(p_angle))
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
  DECLARE_ALIGNED(16, uint8_t, left_col[MAX_SB_SIZE]);
  const uint8_t *above_ref = ref - ref_stride;
  DECLARE_ALIGNED(16, uint8_t, above_data[MAX_SB_SIZE + 16]);
  uint8_t *above_row = above_data + 16;
  const uint8_t *const_above_row = above_row;
  const int bs = 4 * num_4x4_blocks_wide_txsize_lookup[tx_size];
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

  if ((!need_above && n_left_px == 0) || (!need_left && n_top_px == 0)) {
    int i;
    const int val = (n_left_px == 0) ? 129 : 127;
    for (i = 0; i < bs; ++i) {
      memset(dst, val, bs);
      dst += dst_stride;
    }
    return;
  }

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
    if (plane == 0 && vp10_is_intra_filter_switchable(p_angle))
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
  const int txw = num_4x4_blocks_wide_txsize_lookup[tx_size];
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
      mi_col + (1 << mi_width_log2_lookup[bsize]) < xd->tile.mi_col_end;
#if CONFIG_EXT_PARTITION_TYPES
  const PARTITION_TYPE partition = xd->mi[0]->mbmi.partition;
#endif
  const int have_right = vp10_has_right(bsize, mi_row, mi_col,
                                          right_available,
#if CONFIG_EXT_PARTITION_TYPES
                                          partition,
#endif
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
    const int bs = 4 * num_4x4_blocks_wide_txsize_lookup[tx_size];
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
