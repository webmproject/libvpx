/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 *
 */

#include <math.h>

#include "./aom_config.h"
#include "./aom_dsp_rtcd.h"
#include "av1/common/onyxc_int.h"
#include "av1/common/restoration.h"
#include "aom_dsp/aom_dsp_common.h"
#include "aom_mem/aom_mem.h"
#include "aom_ports/mem.h"

#define BILATERAL_PARAM_PRECISION 16
#define BILATERAL_AMP_RANGE 256
#define BILATERAL_AMP_RANGE_SYM (2 * BILATERAL_AMP_RANGE + 1)

static uint8_t bilateral_filter_coeffs_r_kf[BILATERAL_LEVELS_KF]
                                           [BILATERAL_AMP_RANGE_SYM];
static uint8_t bilateral_filter_coeffs_r[BILATERAL_LEVELS]
                                        [BILATERAL_AMP_RANGE_SYM];
static uint8_t bilateral_filter_coeffs_s_kf[BILATERAL_LEVELS_KF]
                                           [RESTORATION_WIN][RESTORATION_WIN];
static uint8_t bilateral_filter_coeffs_s[BILATERAL_LEVELS][RESTORATION_WIN]
                                        [RESTORATION_WIN];

typedef struct bilateral_params {
  int sigma_x;  // spatial variance x
  int sigma_y;  // spatial variance y
  int sigma_r;  // range variance
} BilateralParamsType;

static BilateralParamsType bilateral_level_to_params_arr[BILATERAL_LEVELS] = {
  //  Values are rounded to 1/16 th precision
  { 8, 9, 30 },   { 9, 8, 30 },   { 9, 11, 32 },  { 11, 9, 32 },
  { 14, 14, 36 }, { 18, 18, 36 }, { 24, 24, 40 }, { 32, 32, 40 },
};

static BilateralParamsType
    bilateral_level_to_params_arr_kf[BILATERAL_LEVELS_KF] = {
      //  Values are rounded to 1/16 th precision
      { 8, 8, 30 },   { 9, 9, 32 },   { 10, 10, 32 }, { 12, 12, 32 },
      { 14, 14, 32 }, { 18, 18, 36 }, { 24, 24, 40 }, { 30, 30, 44 },
      { 36, 36, 48 }, { 42, 42, 48 }, { 48, 48, 48 }, { 48, 48, 56 },
      { 56, 56, 48 }, { 56, 56, 56 }, { 56, 56, 64 }, { 64, 64, 48 },
    };

typedef void (*restore_func_type)(uint8_t *data8, int width, int height,
                                  int stride, RestorationInternal *rst,
                                  uint8_t *tmpdata8, int tmpstride);
#if CONFIG_AOM_HIGHBITDEPTH
typedef void (*restore_func_highbd_type)(uint8_t *data8, int width, int height,
                                         int stride, RestorationInternal *rst,
                                         uint8_t *tmpdata8, int tmpstride,
                                         int bit_depth);
#endif  // CONFIG_AOM_HIGHBITDEPTH

static INLINE BilateralParamsType av1_bilateral_level_to_params(int index,
                                                                int kf) {
  return kf ? bilateral_level_to_params_arr_kf[index]
            : bilateral_level_to_params_arr[index];
}

void av1_loop_restoration_precal() {
  int i;
  for (i = 0; i < BILATERAL_LEVELS_KF; i++) {
    const BilateralParamsType param = av1_bilateral_level_to_params(i, 1);
    const int sigma_x = param.sigma_x;
    const int sigma_y = param.sigma_y;
    const int sigma_r = param.sigma_r;
    const double sigma_r_d = (double)sigma_r / BILATERAL_PARAM_PRECISION;
    const double sigma_x_d = (double)sigma_x / BILATERAL_PARAM_PRECISION;
    const double sigma_y_d = (double)sigma_y / BILATERAL_PARAM_PRECISION;

    uint8_t *fr = bilateral_filter_coeffs_r_kf[i] + BILATERAL_AMP_RANGE;
    int j, x, y;
    for (j = 0; j <= BILATERAL_AMP_RANGE; j++) {
      fr[j] = (uint8_t)(0.5 +
                        RESTORATION_FILT_STEP *
                            exp(-(j * j) / (2 * sigma_r_d * sigma_r_d)));
      fr[-j] = fr[j];
    }
    for (y = -RESTORATION_HALFWIN; y <= RESTORATION_HALFWIN; y++) {
      for (x = -RESTORATION_HALFWIN; x <= RESTORATION_HALFWIN; x++) {
        bilateral_filter_coeffs_s_kf[i][y + RESTORATION_HALFWIN]
                                    [x + RESTORATION_HALFWIN] = (uint8_t)(
                                        0.5 +
                                        RESTORATION_FILT_STEP *
                                            exp(-(x * x) / (2 * sigma_x_d *
                                                            sigma_x_d) -
                                                (y * y) / (2 * sigma_y_d *
                                                           sigma_y_d)));
      }
    }
  }
  for (i = 0; i < BILATERAL_LEVELS; i++) {
    const BilateralParamsType param = av1_bilateral_level_to_params(i, 0);
    const int sigma_x = param.sigma_x;
    const int sigma_y = param.sigma_y;
    const int sigma_r = param.sigma_r;
    const double sigma_r_d = (double)sigma_r / BILATERAL_PARAM_PRECISION;
    const double sigma_x_d = (double)sigma_x / BILATERAL_PARAM_PRECISION;
    const double sigma_y_d = (double)sigma_y / BILATERAL_PARAM_PRECISION;

    uint8_t *fr = bilateral_filter_coeffs_r[i] + BILATERAL_AMP_RANGE;
    int j, x, y;
    for (j = 0; j <= BILATERAL_AMP_RANGE; j++) {
      fr[j] = (uint8_t)(0.5 +
                        RESTORATION_FILT_STEP *
                            exp(-(j * j) / (2 * sigma_r_d * sigma_r_d)));
      fr[-j] = fr[j];
    }
    for (y = -RESTORATION_HALFWIN; y <= RESTORATION_HALFWIN; y++) {
      for (x = -RESTORATION_HALFWIN; x <= RESTORATION_HALFWIN; x++) {
        bilateral_filter_coeffs_s[i][y + RESTORATION_HALFWIN]
                                 [x + RESTORATION_HALFWIN] = (uint8_t)(
                                     0.5 +
                                     RESTORATION_FILT_STEP *
                                         exp(-(x * x) /
                                                 (2 * sigma_x_d * sigma_x_d) -
                                             (y * y) /
                                                 (2 * sigma_y_d * sigma_y_d)));
      }
    }
  }
}

int av1_bilateral_level_bits(const AV1_COMMON *const cm) {
  return cm->frame_type == KEY_FRAME ? BILATERAL_LEVEL_BITS_KF
                                     : BILATERAL_LEVEL_BITS;
}

void av1_loop_restoration_init(RestorationInternal *rst, RestorationInfo *rsi,
                               int kf, int width, int height) {
  int i, tile_idx;
  rst->rsi = rsi;
  rst->keyframe = kf;
  rst->subsampling_x = 0;
  rst->subsampling_y = 0;
  rst->ntiles =
      av1_get_rest_ntiles(width, height, &rst->tile_width, &rst->tile_height,
                          &rst->nhtiles, &rst->nvtiles);
  if (rsi->frame_restoration_type == RESTORE_WIENER) {
    for (tile_idx = 0; tile_idx < rst->ntiles; ++tile_idx) {
      rsi->wiener_info[tile_idx].vfilter[RESTORATION_HALFWIN] =
          rsi->wiener_info[tile_idx].hfilter[RESTORATION_HALFWIN] =
              RESTORATION_FILT_STEP;
      for (i = 0; i < RESTORATION_HALFWIN; ++i) {
        rsi->wiener_info[tile_idx].vfilter[RESTORATION_WIN - 1 - i] =
            rsi->wiener_info[tile_idx].vfilter[i];
        rsi->wiener_info[tile_idx].hfilter[RESTORATION_WIN - 1 - i] =
            rsi->wiener_info[tile_idx].hfilter[i];
        rsi->wiener_info[tile_idx].vfilter[RESTORATION_HALFWIN] -=
            2 * rsi->wiener_info[tile_idx].vfilter[i];
        rsi->wiener_info[tile_idx].hfilter[RESTORATION_HALFWIN] -=
            2 * rsi->wiener_info[tile_idx].hfilter[i];
      }
    }
  } else if (rsi->frame_restoration_type == RESTORE_SWITCHABLE) {
    for (tile_idx = 0; tile_idx < rst->ntiles; ++tile_idx) {
      if (rsi->restoration_type[tile_idx] == RESTORE_WIENER) {
        rsi->wiener_info[tile_idx].vfilter[RESTORATION_HALFWIN] =
            rsi->wiener_info[tile_idx].hfilter[RESTORATION_HALFWIN] =
                RESTORATION_FILT_STEP;
        for (i = 0; i < RESTORATION_HALFWIN; ++i) {
          rsi->wiener_info[tile_idx].vfilter[RESTORATION_WIN - 1 - i] =
              rsi->wiener_info[tile_idx].vfilter[i];
          rsi->wiener_info[tile_idx].hfilter[RESTORATION_WIN - 1 - i] =
              rsi->wiener_info[tile_idx].hfilter[i];
          rsi->wiener_info[tile_idx].vfilter[RESTORATION_HALFWIN] -=
              2 * rsi->wiener_info[tile_idx].vfilter[i];
          rsi->wiener_info[tile_idx].hfilter[RESTORATION_HALFWIN] -=
              2 * rsi->wiener_info[tile_idx].hfilter[i];
        }
      }
    }
  }
}

static void loop_bilateral_filter_tile(uint8_t *data, int tile_idx, int width,
                                       int height, int stride,
                                       RestorationInternal *rst,
                                       uint8_t *tmpdata, int tmpstride) {
  int i, j, subtile_idx;
  int h_start, h_end, v_start, v_end;
  const int tile_width = rst->tile_width >> rst->subsampling_x;
  const int tile_height = rst->tile_height >> rst->subsampling_y;

  for (subtile_idx = 0; subtile_idx < BILATERAL_SUBTILES; ++subtile_idx) {
    uint8_t *data_p, *tmpdata_p;
    const int level = rst->rsi->bilateral_info[tile_idx].level[subtile_idx];
    uint8_t(*wx_lut)[RESTORATION_WIN];
    uint8_t *wr_lut_;

    if (level < 0) continue;
    wr_lut_ = (rst->keyframe ? bilateral_filter_coeffs_r_kf[level]
                             : bilateral_filter_coeffs_r[level]) +
              BILATERAL_AMP_RANGE;
    wx_lut = rst->keyframe ? bilateral_filter_coeffs_s_kf[level]
                           : bilateral_filter_coeffs_s[level];

    av1_get_rest_tile_limits(tile_idx, subtile_idx, BILATERAL_SUBTILE_BITS,
                             rst->nhtiles, rst->nvtiles, tile_width,
                             tile_height, width, height, 1, 1, &h_start, &h_end,
                             &v_start, &v_end);

    data_p = data + h_start + v_start * stride;
    tmpdata_p = tmpdata + h_start + v_start * tmpstride;

    for (i = 0; i < (v_end - v_start); ++i) {
      for (j = 0; j < (h_end - h_start); ++j) {
        int x, y, wt;
        int64_t flsum = 0, wtsum = 0;
        uint8_t *data_p2 = data_p + j - RESTORATION_HALFWIN * stride;
        for (y = -RESTORATION_HALFWIN; y <= RESTORATION_HALFWIN; ++y) {
          for (x = -RESTORATION_HALFWIN; x <= RESTORATION_HALFWIN; ++x) {
            wt = (int)wx_lut[y + RESTORATION_HALFWIN][x + RESTORATION_HALFWIN] *
                 (int)wr_lut_[data_p2[x] - data_p[j]];
            wtsum += (int64_t)wt;
            flsum += (int64_t)wt * data_p2[x];
          }
          data_p2 += stride;
        }
        if (wtsum > 0)
          tmpdata_p[j] = clip_pixel((int)((flsum + wtsum / 2) / wtsum));
        else
          tmpdata_p[j] = data_p[j];
      }
      tmpdata_p += tmpstride;
      data_p += stride;
    }
    for (i = v_start; i < v_end; ++i) {
      memcpy(data + i * stride + h_start, tmpdata + i * tmpstride + h_start,
             (h_end - h_start) * sizeof(*data));
    }
  }
}

static void loop_bilateral_filter(uint8_t *data, int width, int height,
                                  int stride, RestorationInternal *rst,
                                  uint8_t *tmpdata, int tmpstride) {
  int tile_idx;
  for (tile_idx = 0; tile_idx < rst->ntiles; ++tile_idx) {
    loop_bilateral_filter_tile(data, tile_idx, width, height, stride, rst,
                               tmpdata, tmpstride);
  }
}

uint8_t hor_sym_filter(uint8_t *d, int *hfilter) {
  int32_t s =
      (1 << (RESTORATION_FILT_BITS - 1)) + d[0] * hfilter[RESTORATION_HALFWIN];
  int i;
  for (i = 1; i <= RESTORATION_HALFWIN; ++i)
    s += (d[i] + d[-i]) * hfilter[RESTORATION_HALFWIN + i];
  return clip_pixel(s >> RESTORATION_FILT_BITS);
}

uint8_t ver_sym_filter(uint8_t *d, int stride, int *vfilter) {
  int32_t s =
      (1 << (RESTORATION_FILT_BITS - 1)) + d[0] * vfilter[RESTORATION_HALFWIN];
  int i;
  for (i = 1; i <= RESTORATION_HALFWIN; ++i)
    s += (d[i * stride] + d[-i * stride]) * vfilter[RESTORATION_HALFWIN + i];
  return clip_pixel(s >> RESTORATION_FILT_BITS);
}

static void loop_wiener_filter_tile(uint8_t *data, int tile_idx, int width,
                                    int height, int stride,
                                    RestorationInternal *rst, uint8_t *tmpdata,
                                    int tmpstride) {
  const int tile_width = rst->tile_width >> rst->subsampling_x;
  const int tile_height = rst->tile_height >> rst->subsampling_y;
  int i, j;
  int h_start, h_end, v_start, v_end;
  uint8_t *data_p, *tmpdata_p;

  if (rst->rsi->wiener_info[tile_idx].level == 0) return;
  // Filter row-wise
  av1_get_rest_tile_limits(tile_idx, 0, 0, rst->nhtiles, rst->nvtiles,
                           tile_width, tile_height, width, height, 1, 0,
                           &h_start, &h_end, &v_start, &v_end);
  data_p = data + h_start + v_start * stride;
  tmpdata_p = tmpdata + h_start + v_start * tmpstride;
  for (i = 0; i < (v_end - v_start); ++i) {
    for (j = 0; j < (h_end - h_start); ++j) {
      *tmpdata_p++ =
          hor_sym_filter(data_p++, rst->rsi->wiener_info[tile_idx].hfilter);
    }
    data_p += stride - (h_end - h_start);
    tmpdata_p += tmpstride - (h_end - h_start);
  }
  // Filter col-wise
  av1_get_rest_tile_limits(tile_idx, 0, 0, rst->nhtiles, rst->nvtiles,
                           tile_width, tile_height, width, height, 0, 1,
                           &h_start, &h_end, &v_start, &v_end);
  data_p = data + h_start + v_start * stride;
  tmpdata_p = tmpdata + h_start + v_start * tmpstride;
  for (i = 0; i < (v_end - v_start); ++i) {
    for (j = 0; j < (h_end - h_start); ++j) {
      *data_p++ = ver_sym_filter(tmpdata_p++, tmpstride,
                                 rst->rsi->wiener_info[tile_idx].vfilter);
    }
    data_p += stride - (h_end - h_start);
    tmpdata_p += tmpstride - (h_end - h_start);
  }
}

static void loop_wiener_filter(uint8_t *data, int width, int height, int stride,
                               RestorationInternal *rst, uint8_t *tmpdata,
                               int tmpstride) {
  int i, tile_idx;
  uint8_t *data_p, *tmpdata_p;

  // Initialize tmp buffer
  data_p = data;
  tmpdata_p = tmpdata;
  for (i = 0; i < height; ++i) {
    memcpy(tmpdata_p, data_p, sizeof(*data_p) * width);
    data_p += stride;
    tmpdata_p += tmpstride;
  }
  for (tile_idx = 0; tile_idx < rst->ntiles; ++tile_idx) {
    loop_wiener_filter_tile(data, tile_idx, width, height, stride, rst, tmpdata,
                            tmpstride);
  }
}

static void loop_switchable_filter(uint8_t *data, int width, int height,
                                   int stride, RestorationInternal *rst,
                                   uint8_t *tmpdata, int tmpstride) {
  int i, tile_idx;
  uint8_t *data_p, *tmpdata_p;

  // Initialize tmp buffer
  data_p = data;
  tmpdata_p = tmpdata;
  for (i = 0; i < height; ++i) {
    memcpy(tmpdata_p, data_p, sizeof(*data_p) * width);
    data_p += stride;
    tmpdata_p += tmpstride;
  }
  for (tile_idx = 0; tile_idx < rst->ntiles; ++tile_idx) {
    if (rst->rsi->restoration_type[tile_idx] == RESTORE_BILATERAL) {
      loop_bilateral_filter_tile(data, tile_idx, width, height, stride, rst,
                                 tmpdata, tmpstride);
    } else if (rst->rsi->restoration_type[tile_idx] == RESTORE_WIENER) {
      loop_wiener_filter_tile(data, tile_idx, width, height, stride, rst,
                              tmpdata, tmpstride);
    }
  }
}

#if CONFIG_AOM_HIGHBITDEPTH
static void loop_bilateral_filter_tile_highbd(uint16_t *data, int tile_idx,
                                              int width, int height, int stride,
                                              RestorationInternal *rst,
                                              uint16_t *tmpdata, int tmpstride,
                                              int bit_depth) {
  const int tile_width = rst->tile_width >> rst->subsampling_x;
  const int tile_height = rst->tile_height >> rst->subsampling_y;
  int i, j, subtile_idx;
  int h_start, h_end, v_start, v_end;
  const int shift = bit_depth - 8;

  for (subtile_idx = 0; subtile_idx < BILATERAL_SUBTILES; ++subtile_idx) {
    uint16_t *data_p, *tmpdata_p;
    const int level = rst->rsi->bilateral_info[tile_idx].level[subtile_idx];
    uint8_t(*wx_lut)[RESTORATION_WIN];
    uint8_t *wr_lut_;

    if (level < 0) continue;
    wr_lut_ = (rst->keyframe ? bilateral_filter_coeffs_r_kf[level]
                             : bilateral_filter_coeffs_r[level]) +
              BILATERAL_AMP_RANGE;
    wx_lut = rst->keyframe ? bilateral_filter_coeffs_s_kf[level]
                           : bilateral_filter_coeffs_s[level];
    av1_get_rest_tile_limits(tile_idx, subtile_idx, BILATERAL_SUBTILE_BITS,
                             rst->nhtiles, rst->nvtiles, tile_width,
                             tile_height, width, height, 1, 1, &h_start, &h_end,
                             &v_start, &v_end);

    data_p = data + h_start + v_start * stride;
    tmpdata_p = tmpdata + h_start + v_start * tmpstride;

    for (i = 0; i < (v_end - v_start); ++i) {
      for (j = 0; j < (h_end - h_start); ++j) {
        int x, y, wt;
        int64_t flsum = 0, wtsum = 0;
        uint16_t *data_p2 = data_p + j - RESTORATION_HALFWIN * stride;
        for (y = -RESTORATION_HALFWIN; y <= RESTORATION_HALFWIN; ++y) {
          for (x = -RESTORATION_HALFWIN; x <= RESTORATION_HALFWIN; ++x) {
            wt = (int)wx_lut[y + RESTORATION_HALFWIN][x + RESTORATION_HALFWIN] *
                 (int)wr_lut_[(data_p2[x] >> shift) - (data_p[j] >> shift)];
            wtsum += (int64_t)wt;
            flsum += (int64_t)wt * data_p2[x];
          }
          data_p2 += stride;
        }
        if (wtsum > 0)
          tmpdata_p[j] =
              clip_pixel_highbd((int)((flsum + wtsum / 2) / wtsum), bit_depth);
        else
          tmpdata_p[j] = data_p[j];
      }
      tmpdata_p += tmpstride;
      data_p += stride;
    }
    for (i = v_start; i < v_end; ++i) {
      memcpy(data + i * stride + h_start, tmpdata + i * tmpstride + h_start,
             (h_end - h_start) * sizeof(*data));
    }
  }
}

static void loop_bilateral_filter_highbd(uint8_t *data8, int width, int height,
                                         int stride, RestorationInternal *rst,
                                         uint8_t *tmpdata8, int tmpstride,
                                         int bit_depth) {
  int tile_idx;
  uint16_t *data = CONVERT_TO_SHORTPTR(data8);
  uint16_t *tmpdata = CONVERT_TO_SHORTPTR(tmpdata8);

  for (tile_idx = 0; tile_idx < rst->ntiles; ++tile_idx) {
    loop_bilateral_filter_tile_highbd(data, tile_idx, width, height, stride,
                                      rst, tmpdata, tmpstride, bit_depth);
  }
}

uint16_t hor_sym_filter_highbd(uint16_t *d, int *hfilter, int bd) {
  int32_t s =
      (1 << (RESTORATION_FILT_BITS - 1)) + d[0] * hfilter[RESTORATION_HALFWIN];
  int i;
  for (i = 1; i <= RESTORATION_HALFWIN; ++i)
    s += (d[i] + d[-i]) * hfilter[RESTORATION_HALFWIN + i];
  return clip_pixel_highbd(s >> RESTORATION_FILT_BITS, bd);
}

uint16_t ver_sym_filter_highbd(uint16_t *d, int stride, int *vfilter, int bd) {
  int32_t s =
      (1 << (RESTORATION_FILT_BITS - 1)) + d[0] * vfilter[RESTORATION_HALFWIN];
  int i;
  for (i = 1; i <= RESTORATION_HALFWIN; ++i)
    s += (d[i * stride] + d[-i * stride]) * vfilter[RESTORATION_HALFWIN + i];
  return clip_pixel_highbd(s >> RESTORATION_FILT_BITS, bd);
}

static void loop_wiener_filter_tile_highbd(uint16_t *data, int tile_idx,
                                           int width, int height, int stride,
                                           RestorationInternal *rst,
                                           uint16_t *tmpdata, int tmpstride,
                                           int bit_depth) {
  const int tile_width = rst->tile_width >> rst->subsampling_x;
  const int tile_height = rst->tile_height >> rst->subsampling_y;
  int h_start, h_end, v_start, v_end;
  int i, j;
  uint16_t *data_p, *tmpdata_p;

  if (rst->rsi->wiener_info[tile_idx].level == 0) return;
  // Filter row-wise
  av1_get_rest_tile_limits(tile_idx, 0, 0, rst->nhtiles, rst->nvtiles,
                           tile_width, tile_height, width, height, 1, 0,
                           &h_start, &h_end, &v_start, &v_end);
  data_p = data + h_start + v_start * stride;
  tmpdata_p = tmpdata + h_start + v_start * tmpstride;
  for (i = 0; i < (v_end - v_start); ++i) {
    for (j = 0; j < (h_end - h_start); ++j) {
      *tmpdata_p++ = hor_sym_filter_highbd(
          data_p++, rst->rsi->wiener_info[tile_idx].hfilter, bit_depth);
    }
    data_p += stride - (h_end - h_start);
    tmpdata_p += tmpstride - (h_end - h_start);
  }
  // Filter col-wise
  av1_get_rest_tile_limits(tile_idx, 0, 0, rst->nhtiles, rst->nvtiles,
                           tile_width, tile_height, width, height, 0, 1,
                           &h_start, &h_end, &v_start, &v_end);
  data_p = data + h_start + v_start * stride;
  tmpdata_p = tmpdata + h_start + v_start * tmpstride;
  for (i = 0; i < (v_end - v_start); ++i) {
    for (j = 0; j < (h_end - h_start); ++j) {
      *data_p++ = ver_sym_filter_highbd(tmpdata_p++, tmpstride,
                                        rst->rsi->wiener_info[tile_idx].vfilter,
                                        bit_depth);
    }
    data_p += stride - (h_end - h_start);
    tmpdata_p += tmpstride - (h_end - h_start);
  }
}

static void loop_wiener_filter_highbd(uint8_t *data8, int width, int height,
                                      int stride, RestorationInternal *rst,
                                      uint8_t *tmpdata8, int tmpstride,
                                      int bit_depth) {
  uint16_t *data = CONVERT_TO_SHORTPTR(data8);
  uint16_t *tmpdata = CONVERT_TO_SHORTPTR(tmpdata8);
  int i, tile_idx;
  uint16_t *data_p, *tmpdata_p;

  // Initialize tmp buffer
  data_p = data;
  tmpdata_p = tmpdata;
  for (i = 0; i < height; ++i) {
    memcpy(tmpdata_p, data_p, sizeof(*data_p) * width);
    data_p += stride;
    tmpdata_p += tmpstride;
  }
  for (tile_idx = 0; tile_idx < rst->ntiles; ++tile_idx) {
    loop_wiener_filter_tile_highbd(data, tile_idx, width, height, stride, rst,
                                   tmpdata, tmpstride, bit_depth);
  }
}

static void loop_switchable_filter_highbd(uint8_t *data8, int width, int height,
                                          int stride, RestorationInternal *rst,
                                          uint8_t *tmpdata8, int tmpstride,
                                          int bit_depth) {
  uint16_t *data = CONVERT_TO_SHORTPTR(data8);
  uint16_t *tmpdata = CONVERT_TO_SHORTPTR(tmpdata8);
  int i, tile_idx;
  uint16_t *data_p, *tmpdata_p;

  // Initialize tmp buffer
  data_p = data;
  tmpdata_p = tmpdata;
  for (i = 0; i < height; ++i) {
    memcpy(tmpdata_p, data_p, sizeof(*data_p) * width);
    data_p += stride;
    tmpdata_p += tmpstride;
  }
  for (tile_idx = 0; tile_idx < rst->ntiles; ++tile_idx) {
    if (rst->rsi->restoration_type[tile_idx] == RESTORE_BILATERAL) {
      loop_bilateral_filter_tile_highbd(data, tile_idx, width, height, stride,
                                        rst, tmpdata, tmpstride, bit_depth);
    } else if (rst->rsi->restoration_type[tile_idx] == RESTORE_WIENER) {
      loop_wiener_filter_tile_highbd(data, tile_idx, width, height, stride, rst,
                                     tmpdata, tmpstride, bit_depth);
    }
  }
}
#endif  // CONFIG_AOM_HIGHBITDEPTH

void av1_loop_restoration_rows(YV12_BUFFER_CONFIG *frame, AV1_COMMON *cm,
                               int start_mi_row, int end_mi_row, int y_only) {
  const int ywidth = frame->y_crop_width;
  const int ystride = frame->y_stride;
  const int uvwidth = frame->uv_crop_width;
  const int uvstride = frame->uv_stride;
  const int ystart = start_mi_row << MI_SIZE_LOG2;
  const int uvstart = ystart >> cm->subsampling_y;
  int yend = end_mi_row << MI_SIZE_LOG2;
  int uvend = yend >> cm->subsampling_y;
  restore_func_type restore_func =
      cm->rst_internal.rsi->frame_restoration_type == RESTORE_BILATERAL
          ? loop_bilateral_filter
          : (cm->rst_internal.rsi->frame_restoration_type == RESTORE_WIENER
                 ? loop_wiener_filter
                 : loop_switchable_filter);
#if CONFIG_AOM_HIGHBITDEPTH
  restore_func_highbd_type restore_func_highbd =
      cm->rst_internal.rsi->frame_restoration_type == RESTORE_BILATERAL
          ? loop_bilateral_filter_highbd
          : (cm->rst_internal.rsi->frame_restoration_type == RESTORE_WIENER
                 ? loop_wiener_filter_highbd
                 : loop_switchable_filter_highbd);
#endif  // CONFIG_AOM_HIGHBITDEPTH
  YV12_BUFFER_CONFIG tmp_buf;

  if (cm->rst_internal.rsi->frame_restoration_type == RESTORE_NONE) return;

  memset(&tmp_buf, 0, sizeof(YV12_BUFFER_CONFIG));

  yend = AOMMIN(yend, cm->height);
  uvend = AOMMIN(uvend, cm->subsampling_y ? (cm->height + 1) >> 1 : cm->height);

  if (aom_realloc_frame_buffer(
          &tmp_buf, cm->width, cm->height, cm->subsampling_x, cm->subsampling_y,
#if CONFIG_AOM_HIGHBITDEPTH
          cm->use_highbitdepth,
#endif
          AOM_BORDER_IN_PIXELS, cm->byte_alignment, NULL, NULL, NULL) < 0)
    aom_internal_error(&cm->error, AOM_CODEC_MEM_ERROR,
                       "Failed to allocate tmp restoration buffer");

#if CONFIG_AOM_HIGHBITDEPTH
  if (cm->use_highbitdepth)
    restore_func_highbd(frame->y_buffer + ystart * ystride, ywidth,
                        yend - ystart, ystride, &cm->rst_internal,
                        tmp_buf.y_buffer + ystart * tmp_buf.y_stride,
                        tmp_buf.y_stride, cm->bit_depth);
  else
#endif  // CONFIG_AOM_HIGHBITDEPTH
    restore_func(frame->y_buffer + ystart * ystride, ywidth, yend - ystart,
                 ystride, &cm->rst_internal,
                 tmp_buf.y_buffer + ystart * tmp_buf.y_stride,
                 tmp_buf.y_stride);
  if (!y_only) {
    cm->rst_internal.subsampling_x = cm->subsampling_x;
    cm->rst_internal.subsampling_y = cm->subsampling_y;
#if CONFIG_AOM_HIGHBITDEPTH
    if (cm->use_highbitdepth) {
      restore_func_highbd(frame->u_buffer + uvstart * uvstride, uvwidth,
                          uvend - uvstart, uvstride, &cm->rst_internal,
                          tmp_buf.u_buffer + uvstart * tmp_buf.uv_stride,
                          tmp_buf.uv_stride, cm->bit_depth);
      restore_func_highbd(frame->v_buffer + uvstart * uvstride, uvwidth,
                          uvend - uvstart, uvstride, &cm->rst_internal,
                          tmp_buf.v_buffer + uvstart * tmp_buf.uv_stride,
                          tmp_buf.uv_stride, cm->bit_depth);
    } else {
#endif  // CONFIG_AOM_HIGHBITDEPTH
      restore_func(frame->u_buffer + uvstart * uvstride, uvwidth,
                   uvend - uvstart, uvstride, &cm->rst_internal,
                   tmp_buf.u_buffer + uvstart * tmp_buf.uv_stride,
                   tmp_buf.uv_stride);
      restore_func(frame->v_buffer + uvstart * uvstride, uvwidth,
                   uvend - uvstart, uvstride, &cm->rst_internal,
                   tmp_buf.v_buffer + uvstart * tmp_buf.uv_stride,
                   tmp_buf.uv_stride);
#if CONFIG_AOM_HIGHBITDEPTH
    }
#endif  // CONFIG_AOM_HIGHBITDEPTH
  }
  aom_free_frame_buffer(&tmp_buf);
}

void av1_loop_restoration_frame(YV12_BUFFER_CONFIG *frame, AV1_COMMON *cm,
                                RestorationInfo *rsi, int y_only,
                                int partial_frame) {
  int start_mi_row, end_mi_row, mi_rows_to_filter;
  if (rsi->frame_restoration_type != RESTORE_NONE) {
    start_mi_row = 0;
    mi_rows_to_filter = cm->mi_rows;
    if (partial_frame && cm->mi_rows > 8) {
      start_mi_row = cm->mi_rows >> 1;
      start_mi_row &= 0xfffffff8;
      mi_rows_to_filter = AOMMAX(cm->mi_rows / 8, 8);
    }
    end_mi_row = start_mi_row + mi_rows_to_filter;
    av1_loop_restoration_init(&cm->rst_internal, rsi,
                              cm->frame_type == KEY_FRAME, cm->width,
                              cm->height);
    av1_loop_restoration_rows(frame, cm, start_mi_row, end_mi_row, y_only);
  }
}
