/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
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
#include "vp10/common/onyxc_int.h"
#include "vp10/common/restoration.h"
#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem.h"

#define RESTORATION_PARAM_PRECISION     16
#define RESTORATION_RANGE               256
#define RESTORATION_RANGE_SYM           (2 * RESTORATION_RANGE + 1)

static uint8_t restoration_filters_r_kf[RESTORATION_LEVELS_KF]
                                       [RESTORATION_RANGE_SYM];
static uint8_t restoration_filters_r[RESTORATION_LEVELS]
                                    [RESTORATION_RANGE_SYM];
static uint8_t restoration_filters_s_kf[RESTORATION_LEVELS_KF]
                                       [RESTORATION_WIN][RESTORATION_WIN];
static uint8_t restoration_filters_s[RESTORATION_LEVELS]
                                    [RESTORATION_WIN][RESTORATION_WIN];

typedef struct restoration_params {
  int sigma_x;  // spatial variance x
  int sigma_y;  // spatial variance y
  int sigma_r;  // range variance
} RestorationParamsType;

static RestorationParamsType
    restoration_level_to_params_arr[RESTORATION_LEVELS] = {
  // Values are rounded to 1/16 th precision
  {8, 9, 30},
  {9, 8, 30},
  {9, 11, 32},
  {11, 9, 32},
  {14, 14, 32},
  {18, 18, 36},
  {24, 24, 40},
  {32, 32, 40},
};

static RestorationParamsType
    restoration_level_to_params_arr_kf[RESTORATION_LEVELS_KF] = {
  // Values are rounded to 1/16 th precision
  {8, 8, 30},
  {9, 9, 32},
  {10, 10, 32},
  {12, 12, 32},
  {14, 14, 32},
  {18, 18, 36},
  {24, 24, 40},
  {30, 30, 44},
  {36, 36, 48},
  {42, 42, 48},
  {48, 48, 48},
  {48, 48, 56},
  {56, 56, 48},
  {56, 56, 56},
  {56, 56, 64},
  {64, 64, 48},
};

typedef void (*restore_func_type)(
    uint8_t *data8, int width, int height,
    int stride, RestorationInternal *rst,
    uint8_t *tmpdata8, int tmpstride);

#if CONFIG_VPX_HIGHBITDEPTH
typedef void (*restore_func_highbd_type)(
    uint8_t *data8, int width, int height,
    int stride, RestorationInternal *rst,
    uint8_t *tmpdata8, int tmpstride,
    int bit_depth);
#endif  // CONFIG_VPX_HIGHBITDEPTH

static INLINE RestorationParamsType vp10_restoration_level_to_params(
    int index, int kf) {
  return kf ? restoration_level_to_params_arr_kf[index] :
              restoration_level_to_params_arr[index];
}

void vp10_loop_restoration_precal() {
  int i;
  for (i = 0; i < RESTORATION_LEVELS_KF; i ++) {
    const RestorationParamsType param = vp10_restoration_level_to_params(i, 1);
    const int sigma_x = param.sigma_x;
    const int sigma_y = param.sigma_y;
    const int sigma_r = param.sigma_r;
    const double sigma_r_d = (double)sigma_r / RESTORATION_PARAM_PRECISION;
    const double sigma_x_d = (double)sigma_x / RESTORATION_PARAM_PRECISION;
    const double sigma_y_d = (double)sigma_y / RESTORATION_PARAM_PRECISION;

    uint8_t *fr = restoration_filters_r_kf[i] + RESTORATION_RANGE;
    int j, x, y;
    for (j = 0; j <= RESTORATION_RANGE; j++) {
      fr[j] = (uint8_t)(0.5 + RESTORATION_FILT_STEP *
                        exp(-(j * j) / (2 * sigma_r_d * sigma_r_d)));
      fr[-j] = fr[j];
    }
    for (y = -RESTORATION_HALFWIN; y <= RESTORATION_HALFWIN; y++) {
      for (x = -RESTORATION_HALFWIN; x <= RESTORATION_HALFWIN; x++) {
        restoration_filters_s_kf[i][y + RESTORATION_HALFWIN]
                                   [x + RESTORATION_HALFWIN] =
          (uint8_t)(0.5 + RESTORATION_FILT_STEP *
                    exp(-(x * x) / (2 * sigma_x_d * sigma_x_d)
                        -(y * y) / (2 * sigma_y_d * sigma_y_d)));
      }
    }
  }
  for (i = 0; i < RESTORATION_LEVELS; i ++) {
    const RestorationParamsType param = vp10_restoration_level_to_params(i, 0);
    const int sigma_x = param.sigma_x;
    const int sigma_y = param.sigma_y;
    const int sigma_r = param.sigma_r;
    const double sigma_r_d = (double)sigma_r / RESTORATION_PARAM_PRECISION;
    const double sigma_x_d = (double)sigma_x / RESTORATION_PARAM_PRECISION;
    const double sigma_y_d = (double)sigma_y / RESTORATION_PARAM_PRECISION;

    uint8_t *fr = restoration_filters_r[i] + RESTORATION_RANGE;
    int j, x, y;
    for (j = 0; j <= RESTORATION_RANGE; j++) {
      fr[j] = (uint8_t)(0.5 + RESTORATION_FILT_STEP *
                        exp(-(j * j) / (2 * sigma_r_d * sigma_r_d)));
      fr[-j] = fr[j];
    }
    for (y = -RESTORATION_HALFWIN; y <= RESTORATION_HALFWIN; y++) {
      for (x = -RESTORATION_HALFWIN; x <= RESTORATION_HALFWIN; x++) {
        restoration_filters_s[i][y + RESTORATION_HALFWIN]
                                [x + RESTORATION_HALFWIN] =
            (uint8_t)(0.5 + RESTORATION_FILT_STEP *
                      exp(-(x * x) / (2 * sigma_x_d * sigma_x_d)
                          -(y * y) / (2 * sigma_y_d * sigma_y_d)));
      }
    }
  }
}

int vp10_restoration_level_bits(const VP10_COMMON *const cm) {
  return cm->frame_type == KEY_FRAME ?
      RESTORATION_LEVEL_BITS_KF : RESTORATION_LEVEL_BITS;
}

void vp10_loop_restoration_init(RestorationInternal *rst,
                                RestorationInfo *rsi, int kf) {
  int i;
  rst->restoration_type = rsi->restoration_type;
  if (rsi->restoration_type == RESTORE_BILATERAL) {
    const int level = rsi->restoration_level;
    assert(level >= 0);
    rst->wr_lut = kf ? restoration_filters_r_kf[level] :
                       restoration_filters_r[level];
    for (i = 0; i < RESTORATION_WIN; i++)
      rst->wx_lut[i] = kf ? restoration_filters_s_kf[level][i] :
                            restoration_filters_s[level][i];
  } else if (rsi->restoration_type == RESTORE_WIENER) {
    rst->vfilter[RESTORATION_HALFWIN] = rst->hfilter[RESTORATION_HALFWIN] =
        RESTORATION_FILT_STEP;
    for (i = 0; i < RESTORATION_HALFWIN; ++i) {
      rst->vfilter[i] = rst->vfilter[RESTORATION_WIN - 1 - i] = rsi->vfilter[i];
      rst->hfilter[i] = rst->hfilter[RESTORATION_WIN - 1 - i] = rsi->hfilter[i];
      rst->vfilter[RESTORATION_HALFWIN] -= 2 * rsi->vfilter[i];
      rst->hfilter[RESTORATION_HALFWIN] -= 2 * rsi->hfilter[i];
    }
  }
}

static void loop_bilateral_filter(uint8_t *data, int width, int height,
                                  int stride, RestorationInternal *rst,
                                  uint8_t *tmpdata, int tmpstride) {
  int i, j;
  const uint8_t *wr_lut_ = rst->wr_lut + RESTORATION_RANGE;

  uint8_t *data_p = data + RESTORATION_HALFWIN * stride;
  uint8_t *tmpdata_p = tmpdata + RESTORATION_HALFWIN * tmpstride;
  for (i = RESTORATION_HALFWIN; i < height - RESTORATION_HALFWIN; ++i) {
    for (j = RESTORATION_HALFWIN; j < width - RESTORATION_HALFWIN; ++j) {
      int x, y;
      int flsum = 0, wtsum = 0, wt;
      uint8_t *data_p2 = data_p + j - RESTORATION_HALFWIN * stride;
      for (y = -RESTORATION_HALFWIN; y <= RESTORATION_HALFWIN; ++y) {
        for (x = -RESTORATION_HALFWIN; x <= RESTORATION_HALFWIN; ++x) {
          wt = (int)rst->wx_lut[y + RESTORATION_HALFWIN]
                               [x + RESTORATION_HALFWIN] *
               (int)wr_lut_[data_p2[x] - data_p[j]];
          wtsum += wt;
          flsum += wt * data_p2[x];
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

  for (i = RESTORATION_HALFWIN; i < height - RESTORATION_HALFWIN; ++i) {
    memcpy(data + i * stride + RESTORATION_HALFWIN,
           tmpdata + i * tmpstride + RESTORATION_HALFWIN,
           (width - RESTORATION_HALFWIN * 2) * sizeof(*data));
  }
}

uint8_t hor_sym_filter(uint8_t *d, int *hfilter) {
  int32_t s = (1 << (RESTORATION_FILT_BITS - 1)) +
      d[0] * hfilter[RESTORATION_HALFWIN];
  int i;
  for (i = 1; i <= RESTORATION_HALFWIN; ++i)
    s += (d[i] + d[-i]) * hfilter[RESTORATION_HALFWIN + i];
  return clip_pixel(s >> RESTORATION_FILT_BITS);
}

uint8_t ver_sym_filter(uint8_t *d, int stride, int *vfilter) {
  int32_t s = (1 << (RESTORATION_FILT_BITS - 1)) +
      d[0] * vfilter[RESTORATION_HALFWIN];
  int i;
  for (i = 1; i <= RESTORATION_HALFWIN; ++i)
    s += (d[i * stride] + d[-i * stride]) * vfilter[RESTORATION_HALFWIN + i];
  return clip_pixel(s >> RESTORATION_FILT_BITS);
}

static void loop_wiener_filter(uint8_t *data, int width, int height,
                               int stride, RestorationInternal *rst,
                               uint8_t *tmpdata, int tmpstride) {
  uint8_t *data_p = data;
  uint8_t *tmpdata_p = tmpdata;
  int i, j;

  for (i = 0; i < height; ++i) {
    memcpy(tmpdata_p, data_p, sizeof(*data_p) * RESTORATION_HALFWIN);
    data_p += RESTORATION_HALFWIN;
    tmpdata_p += RESTORATION_HALFWIN;
    for (j = RESTORATION_HALFWIN; j < width - RESTORATION_HALFWIN; ++j) {
      *tmpdata_p++ = hor_sym_filter(data_p++, rst->hfilter);
    }
    memcpy(tmpdata_p, data_p, sizeof(*data_p) * RESTORATION_HALFWIN);
    data_p += RESTORATION_HALFWIN - width + stride;
    tmpdata_p += RESTORATION_HALFWIN - width + tmpstride;
  }
  data_p = data;
  tmpdata_p = tmpdata;
  for (i = 0; i < RESTORATION_HALFWIN; ++i) {
    memcpy(data_p, tmpdata_p, sizeof(*data_p) * width);
    data_p += stride;
    tmpdata_p += tmpstride;
  }
  for (; i < height - RESTORATION_HALFWIN; ++i) {
    for (j = 0; j < width; ++j)
      *data_p++ = ver_sym_filter(tmpdata_p++, tmpstride, rst->vfilter);
    data_p += stride - width;
    tmpdata_p += tmpstride - width;
  }
  for (; i < height; ++i) {
    memcpy(data_p, tmpdata_p, sizeof(*data_p) * width);
    data_p += stride;
    tmpdata_p += tmpstride;
  }
}

#if CONFIG_VPX_HIGHBITDEPTH
static void loop_bilateral_filter_highbd(
    uint8_t *data8, int width, int height,
    int stride, RestorationInternal *rst,
    uint8_t *tmpdata8, int tmpstride, int bit_depth) {
  int i, j;
  const uint8_t *wr_lut_ = rst->wr_lut + RESTORATION_RANGE;

  uint16_t *data = CONVERT_TO_SHORTPTR(data8);
  uint16_t *tmpdata = CONVERT_TO_SHORTPTR(tmpdata8);
  uint16_t *data_p = data + RESTORATION_HALFWIN * stride;
  uint16_t *tmpdata_p = tmpdata + RESTORATION_HALFWIN * tmpstride;
  for (i = RESTORATION_HALFWIN; i < height - RESTORATION_HALFWIN; ++i) {
    for (j = RESTORATION_HALFWIN; j < width - RESTORATION_HALFWIN; ++j) {
      int x, y, diff_r;
      int flsum = 0, wtsum = 0, wt;
      uint16_t *data_p2 = data_p + j - RESTORATION_HALFWIN * stride;
      for (y = -RESTORATION_HALFWIN; y <= RESTORATION_HALFWIN; ++y) {
        for (x = -RESTORATION_HALFWIN; x <= RESTORATION_HALFWIN; ++x) {
          diff_r = (data_p2[x] - data_p[j]) >> (bit_depth - 8);
          assert(diff_r >= -RESTORATION_RANGE && diff_r <= RESTORATION_RANGE);
          wt = (int)rst->wx_lut[y + RESTORATION_HALFWIN]
                               [x + RESTORATION_HALFWIN] *
               (int)wr_lut_[diff_r];
          wtsum += wt;
          flsum += wt * data_p2[x];
        }
        data_p2 += stride;
      }
      if (wtsum > 0)
        tmpdata_p[j] = clip_pixel_highbd((int)((flsum + wtsum / 2) / wtsum),
                                         bit_depth);
      else
        tmpdata_p[j] = data_p[j];
    }
    tmpdata_p += tmpstride;
    data_p += stride;
  }

  for (i = RESTORATION_HALFWIN; i < height - RESTORATION_HALFWIN; ++i) {
    memcpy(data + i * stride + RESTORATION_HALFWIN,
           tmpdata + i * tmpstride + RESTORATION_HALFWIN,
           (width - RESTORATION_HALFWIN * 2) * sizeof(*data));
  }
}

uint16_t hor_sym_filter_highbd(uint16_t *d, int *hfilter, int bd) {
  int32_t s = (1 << (RESTORATION_FILT_BITS - 1)) +
      d[0] * hfilter[RESTORATION_HALFWIN];
  int i;
  for (i = 1; i <= RESTORATION_HALFWIN; ++i)
    s += (d[i] + d[-i]) * hfilter[RESTORATION_HALFWIN + i];
  return clip_pixel_highbd(s >> RESTORATION_FILT_BITS, bd);
}

uint16_t ver_sym_filter_highbd(uint16_t *d, int stride, int *vfilter, int bd) {
  int32_t s = (1 << (RESTORATION_FILT_BITS - 1)) +
      d[0] * vfilter[RESTORATION_HALFWIN];
  int i;
  for (i = 1; i <= RESTORATION_HALFWIN; ++i)
    s += (d[i * stride] + d[-i * stride]) * vfilter[RESTORATION_HALFWIN + i];
  return clip_pixel_highbd(s >> RESTORATION_FILT_BITS, bd);
}

static void loop_wiener_filter_highbd(uint8_t *data8, int width, int height,
                                      int stride, RestorationInternal *rst,
                                      uint8_t *tmpdata8, int tmpstride,
                                      int bit_depth) {
  uint16_t *data = CONVERT_TO_SHORTPTR(data8);
  uint16_t *tmpdata = CONVERT_TO_SHORTPTR(tmpdata8);
  uint16_t *data_p = data;
  uint16_t *tmpdata_p = tmpdata;
  int i, j;
  for (i = 0; i < height; ++i) {
    memcpy(tmpdata_p, data_p, sizeof(*data_p) * RESTORATION_HALFWIN);
    data_p += RESTORATION_HALFWIN;
    tmpdata_p += RESTORATION_HALFWIN;
    for (j = RESTORATION_HALFWIN; j < width - RESTORATION_HALFWIN; ++j) {
      *tmpdata_p++ = hor_sym_filter_highbd(data_p++, rst->hfilter, bit_depth);
    }
    memcpy(tmpdata_p, data_p, sizeof(*data_p) * RESTORATION_HALFWIN);
    data_p += RESTORATION_HALFWIN - width + stride;
    tmpdata_p += RESTORATION_HALFWIN - width + tmpstride;
  }
  data_p = data;
  tmpdata_p = tmpdata;
  for (i = 0; i < RESTORATION_HALFWIN; ++i) {
    memcpy(data_p, tmpdata_p, sizeof(*data_p) * width);
    data_p += stride;
    tmpdata_p += tmpstride;
  }
  for (; i < height - RESTORATION_HALFWIN; ++i) {
    for (j = 0; j < width; ++j)
      *data_p++ = ver_sym_filter_highbd(
          tmpdata_p++, tmpstride, rst->vfilter, bit_depth);
    data_p += stride - width;
    tmpdata_p += tmpstride - width;
  }
  for (; i < height; ++i) {
    memcpy(data_p, tmpdata_p, sizeof(*data_p) * width);
    data_p += stride;
    tmpdata_p += tmpstride;
  }
}
#endif  // CONFIG_VPX_HIGHBITDEPTH

void vp10_loop_restoration_rows(YV12_BUFFER_CONFIG *frame,
                                VP10_COMMON *cm,
                                int start_mi_row, int end_mi_row,
                                int y_only) {
  const int ywidth = frame->y_crop_width;
  const int ystride = frame->y_stride;
  const int uvwidth = frame->uv_crop_width;
  const int uvstride = frame->uv_stride;
  const int ystart = start_mi_row << MI_SIZE_LOG2;
  const int uvstart = ystart >> cm->subsampling_y;
  int yend = end_mi_row << MI_SIZE_LOG2;
  int uvend = yend >> cm->subsampling_y;
  restore_func_type restore_func =
      cm->rst_internal.restoration_type == RESTORE_BILATERAL ?
      loop_bilateral_filter : loop_wiener_filter;
#if CONFIG_VPX_HIGHBITDEPTH
  restore_func_highbd_type restore_func_highbd =
      cm->rst_internal.restoration_type == RESTORE_BILATERAL ?
      loop_bilateral_filter_highbd : loop_wiener_filter_highbd;
#endif  // CONFIG_VPX_HIGHBITDEPTH
  YV12_BUFFER_CONFIG *tmp_buf;

  yend = VPXMIN(yend, cm->height);
  uvend = VPXMIN(uvend, cm->subsampling_y ? (cm->height + 1) >> 1 : cm->height);

  if (vpx_realloc_frame_buffer(&cm->tmp_loop_buf, cm->width, cm->height,
                               cm->subsampling_x, cm->subsampling_y,
#if CONFIG_VPX_HIGHBITDEPTH
                               cm->use_highbitdepth,
#endif
                               VPX_DEC_BORDER_IN_PIXELS, cm->byte_alignment,
                               NULL, NULL, NULL) < 0)
    vpx_internal_error(&cm->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate tmp restoration buffer");

  tmp_buf = &cm->tmp_loop_buf;

#if CONFIG_VPX_HIGHBITDEPTH
  if (cm->use_highbitdepth)
    restore_func_highbd(
        frame->y_buffer + ystart * ystride,
        ywidth, yend - ystart, ystride, &cm->rst_internal,
        tmp_buf->y_buffer + ystart * tmp_buf->y_stride,
        tmp_buf->y_stride, cm->bit_depth);
  else
#endif  // CONFIG_VPX_HIGHBITDEPTH
    restore_func(
        frame->y_buffer + ystart * ystride,
        ywidth, yend - ystart, ystride, &cm->rst_internal,
        tmp_buf->y_buffer + ystart * tmp_buf->y_stride,
        tmp_buf->y_stride);
  if (!y_only) {
#if CONFIG_VPX_HIGHBITDEPTH
    if (cm->use_highbitdepth) {
      restore_func_highbd(
          frame->u_buffer + uvstart * uvstride,
          uvwidth, uvend - uvstart, uvstride, &cm->rst_internal,
          tmp_buf->u_buffer + uvstart * tmp_buf->uv_stride,
          tmp_buf->uv_stride, cm->bit_depth);
      restore_func_highbd(
          frame->v_buffer + uvstart * uvstride,
          uvwidth, uvend - uvstart, uvstride, &cm->rst_internal,
          tmp_buf->v_buffer + uvstart * tmp_buf->uv_stride,
          tmp_buf->uv_stride, cm->bit_depth);
    } else {
#endif  // CONFIG_VPX_HIGHBITDEPTH
      restore_func(
          frame->u_buffer + uvstart * uvstride,
          uvwidth, uvend - uvstart, uvstride, &cm->rst_internal,
          tmp_buf->u_buffer + uvstart * tmp_buf->uv_stride,
          tmp_buf->uv_stride);
      restore_func(
          frame->v_buffer + uvstart * uvstride,
          uvwidth, uvend - uvstart, uvstride, &cm->rst_internal,
          tmp_buf->v_buffer + uvstart * tmp_buf->uv_stride,
          tmp_buf->uv_stride);
#if CONFIG_VPX_HIGHBITDEPTH
    }
#endif  // CONFIG_VPX_HIGHBITDEPTH
  }
}

void vp10_loop_restoration_frame(YV12_BUFFER_CONFIG *frame,
                                 VP10_COMMON *cm,
                                 RestorationInfo *rsi,
                                 int y_only, int partial_frame) {
  int start_mi_row, end_mi_row, mi_rows_to_filter;
  if (rsi->restoration_type != RESTORE_NONE) {
    start_mi_row = 0;
    mi_rows_to_filter = cm->mi_rows;
    if (partial_frame && cm->mi_rows > 8) {
      start_mi_row = cm->mi_rows >> 1;
      start_mi_row &= 0xfffffff8;
      mi_rows_to_filter = VPXMAX(cm->mi_rows / 8, 8);
    }
    end_mi_row = start_mi_row + mi_rows_to_filter;
    vp10_loop_restoration_init(&cm->rst_internal, rsi,
                               cm->frame_type == KEY_FRAME);
    vp10_loop_restoration_rows(frame, cm, start_mi_row, end_mi_row, y_only);
  }
}
