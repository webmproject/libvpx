/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "vp9/common/vp9_sr_txfm.h"
#include "vp9/common/vp9_idct.h"
#include "vp9/common/vp9_idwt.h"

#if CONFIG_SR_MODE
int is_enable_srmode(BLOCK_SIZE bsize) {
  TX_SIZE max_tx_size = max_txsize_lookup[bsize];
  return (max_tx_size >= MIN_SR_TX_SIZE &&
          max_tx_size <= MAX_SR_TX_SIZE);
}

// Extending blocks border by copying the boundary pixels
// For convolution use
static void sr_extend(int16_t *src, int src_stride, int w, int h,
                      int16_t *src_ext, int src_ext_stride, int border) {
  int16_t *src_ext_ori = src_ext;
  int i, j;

  src_ext = src_ext_ori - border;
  for (i = 0; i < h; i ++) {
    for (j = 0; j < border; j ++)
      src_ext[j] = src[0];
    vpx_memcpy(src_ext + border, src, sizeof(int16_t) * w);
    for (j = 0; j < border; j ++)
      src_ext[border + w + j] = src[w - 1];
    src_ext += src_ext_stride;
    src += src_stride;
  }

  src_ext = src_ext_ori - border * src_ext_stride - border;
  for (i = 0; i < border; i ++)
    vpx_memcpy(src_ext + i * src_ext_stride, src_ext_ori - border,
               sizeof(int16_t) * (w + 2 * border));

  src_ext = src_ext_ori + h * src_ext_stride - border;
  for (i = 0; i < border; i ++)
    vpx_memcpy(src_ext + i * src_ext_stride,
               src_ext_ori + (h - 1) * src_ext_stride - border,
               sizeof(int16_t) * (w + 2 * border));
}

static void convolve_horiz(int16_t *src, int src_stride, int src_offset,
                           int16_t *dst, int dst_stride, int dst_offset,
                           int *x_filter, int filter_taps, int fil_offset,
                           int w, int h) {
  // src_offset, dst_offset: offsets to move to the next value (usually 1)
  // fil_offset: offsets from filter center to the first pixel of filter
  // (e.g: 3-tap filter (1, 2, 1),    fil_offset=1;
  //       4-tap filter (1, 2, 2, 1), fil_offset=1)
  int x, y;
  int round_offset = 1 << (UPSCALE_FILTER_SHIFT - 1);

  // Shift the buffer to the first pixel of filter
  src -= (fil_offset * src_offset);

  for (y = 0; y < h; ++y) {
    int16_t *src_x = src;
    for (x = 0; x < w; ++x) {
      int k, sum = 0;

      // If the filter is symmetric, can fold it first, then multiply
      for (k = 0; k < filter_taps; ++k)
        sum += src_x[k * src_offset] * x_filter[k];
      src_x += src_offset;
      dst[x * dst_offset] = (sum + round_offset) >> UPSCALE_FILTER_SHIFT;
    }
    src += src_stride;
    dst += dst_stride;
  }
}

static void convolve_vert(int16_t *src, int src_stride, int src_offset,
                          int16_t *dst, int dst_stride, int dst_offset,
                          int *y_filter, int filter_taps, int fil_offset,
                          int w, int h) {
  int x, y;
  int round_offset = 1 << (UPSCALE_FILTER_SHIFT - 1);

  // Shift the buffer to the first pixel of filter
  src -= src_stride * fil_offset;

  for (x = 0; x < w; ++x) {
    int16_t *src_y = src;
    for (y = 0; y < h; ++y) {
      int k, sum = 0;
      for (k = 0; k < filter_taps; ++k)
        sum += src_y[k * src_stride] * y_filter[k];
      src_y += src_stride;
      dst[y * dst_stride] = (sum + round_offset) >> UPSCALE_FILTER_SHIFT;
    }
    src += src_offset;
    dst += dst_offset;
  }
}

/* If changing filter taps, need to change some parameters below too */
int lp_filter[UPSCALE_FILTER_TAPS - 1] = {2, -14, -2, 43, 70, 43, -2, -14, 2};
int interpl_filter[UPSCALE_FILTER_TAPS] =
    // {1, -4, 10, -23, 80, 80, -23, 10, -4, 1};  // lanczos
    // {0, -1,  6, -19, 78, 78, -19,  6, -1, 0};  // laplacian
    {0, -4, 11, -23, 80, 80, -23, 11, -4, 0};  // DCT
    // {0, -1, -4,  14, 55, 55,  14, -4, -1, 0};  // freqmultiplier = 0.5

#if SR_USE_MULTI_F
// For multiple interplation filter options
int post_filter[SR_USFILTER_NUM_D][UPSCALE_FILTER_TAPS - 1] = {
      {2, 0, -7, 0, 137, 0, -7, 0, 2},
      {1, 0, -7, 0, 139, 0, -7, 0, 1},
      {4, 0, -5, 0, 134, 0, -5, 0, 4},
      {-2, 0, -8, 0, 149, 0, -8, 0, -2}
    };
#else
int post_filter[UPSCALE_FILTER_TAPS - 1] =
    {2, 0, -7, 0, 137, 0, -7, 0, 2};
#endif

#define buf_size (64 + UPSCALE_FILTER_TAPS*2)
static void sr_convolution(int16_t *src, int src_stride, int src_offset,
                           int16_t *dst, int dst_stride, int dst_offset,
                           int * fil_hor, int * fil_ver,
                           int fil_offset_h, int fil_offset_v,
                           int filter_taps,
                           int w, int h) {
  DECLARE_ALIGNED_ARRAY(16, int16_t, tmp_buf, buf_size * buf_size);
  int tmp_buf_stride = buf_size;
  int16_t *tmp_buf_ori = tmp_buf + fil_offset_v * tmp_buf_stride + fil_offset_h;

  convolve_horiz(
      src - fil_offset_v * src_stride, src_stride, src_offset,
      tmp_buf_ori - fil_offset_v * tmp_buf_stride, tmp_buf_stride, 1,
      fil_hor, filter_taps, fil_offset_h, w, h + filter_taps);

  convolve_vert(
      tmp_buf_ori, tmp_buf_stride, 1,
      dst, dst_stride, dst_offset,
      fil_ver, filter_taps, fil_offset_v, w, h);
}


void sr_lowpass(int16_t *src, int src_stride, int16_t *dst, int dst_stride,
                int w, int h) {
  int filter_taps = UPSCALE_FILTER_TAPS - 1;  // odd number of taps
  int border = (filter_taps - 1) >> 1;
  int fil_offset = border;  // see "fil_offset" in "convolve_horiz"

  DECLARE_ALIGNED_ARRAY(16, int16_t, src_ext, buf_size * buf_size);
  int src_ext_stride = buf_size;
  int16_t *src_ext_ori = src_ext + border * src_ext_stride + border;

  // extension
  sr_extend(src, src_stride, w, h, src_ext_ori, src_ext_stride, border);
  sr_convolution(src_ext_ori, src_ext_stride, 1, dst, dst_stride, 1,
                 lp_filter, lp_filter, fil_offset, fil_offset,
                 filter_taps, w, h);
}

static void sr_upsample(int16_t *src, int src_stride,
                        int16_t *dst, int dst_stride, int w, int h) {
  // Apply interpolation filter
  int filter_taps = UPSCALE_FILTER_TAPS;  // even number of taps
  int border = (filter_taps >> 1);  // maximum pixels needs to extended
  int fil_offset = border - 1;  // see "fil_offset" in "convolve_horiz"

  int i, j;
  DECLARE_ALIGNED_ARRAY(16, int16_t, src_ext, buf_size * buf_size);
  DECLARE_ALIGNED_ARRAY(16, int16_t, tmp_buf, buf_size * buf_size);
  int src_ext_stride = buf_size, tmp_buf_stride = buf_size;
  int16_t *src_ext_ori = src_ext + border * src_ext_stride + border;
  int16_t *tmp_buf_ori = tmp_buf + border * tmp_buf_stride + border;
  int16_t *dst_ori = dst;

  // Extend the buffer
  sr_extend(src, src_stride, w, h, src_ext_ori, src_ext_stride, border);

  // Keep the original pixels the same
  for (i = 0; i < h; i++) {
    for (j = 0; j < w; j++) {
      dst[j << 1] = src[j];
    }
    dst += (dst_stride << 1);
    src += src_stride;
  }

  convolve_horiz(src_ext_ori - border * src_ext_stride, src_ext_stride, 1,
                 tmp_buf_ori - border * tmp_buf_stride, tmp_buf_stride, 1,
                 interpl_filter, filter_taps, fil_offset, w, h + (2 * border));

  // Set the horizontally interpolated pixels
  dst = dst_ori;
  tmp_buf = tmp_buf_ori;
  for (i = 0; i < h; i++) {
    for (j = 0; j < w; j++) {
      dst[(j << 1) + 1] = tmp_buf[j];
    }
    dst += (dst_stride << 1);
    tmp_buf += tmp_buf_stride;
  }

  // Set the vertically interpolated pixels
  convolve_vert(src_ext_ori, src_ext_stride, 1,
                dst_ori + dst_stride, 2 * dst_stride, 2,
                interpl_filter, filter_taps, fil_offset, w, h);

  // Set the horizontally and vertically interpolated pixels
  convolve_vert(tmp_buf_ori, tmp_buf_stride, 1,
                dst_ori + dst_stride + 1, 2 * dst_stride, 2,
                interpl_filter, filter_taps, fil_offset, w, h);
}

#if SR_USE_MULTI_F
static void sr_post_filter(int16_t *src, int src_stride,
                           int16_t *dst, int dst_stride,
                           int w, int h, int f_hor, int f_ver) {
#else
static void sr_post_filter(int16_t *src, int src_stride,
                           int16_t *dst, int dst_stride, int w, int h) {
#endif
  int filter_taps = UPSCALE_FILTER_TAPS - 1;  // odd number of taps
  int border = (filter_taps - 1) >> 1;
  int fil_offset = border;

  DECLARE_ALIGNED_ARRAY(16, int16_t, src_ext, buf_size * buf_size);
  int src_ext_stride = buf_size;
  int16_t *src_ext_ori = src_ext + border * src_ext_stride + border;

  // extension
  sr_extend(src, src_stride, w, h, src_ext_ori, src_ext_stride, border);
  sr_convolution(src_ext_ori, src_ext_stride, 1, dst, dst_stride, 1,
#if SR_USE_MULTI_F
                 post_filter[f_hor], post_filter[f_ver],
                 fil_offset, fil_offset,
#else
                 post_filter, post_filter, fil_offset, fil_offset,
#endif
                 filter_taps, w, h);
}

#if SR_USE_MULTI_F
void sr_recon(int16_t *src, int src_stride, uint8_t *dst, int dst_stride,
              int w, int h, int f_hor, int f_ver) {
#else
void sr_recon(int16_t *src, int src_stride, uint8_t *dst, int dst_stride,
              int w, int h) {
#endif
  DECLARE_ALIGNED_ARRAY(16, int16_t, recon, 64 * 64);
  DECLARE_ALIGNED_ARRAY(16, int16_t, us_resi, 64 * 64);
  int i, j, us_resi_stride = 64;
  int recon_stride = 64;
#if USE_POST_F
  DECLARE_ALIGNED_ARRAY(16, int16_t, enh_recon, 64 * 64);
  int enh_recon_stride = 64;
#endif

  // Upsample residual
  sr_upsample(src, src_stride, us_resi, us_resi_stride, w/2, h/2);

  // Add upsampled residual to prediction
  for (i = 0; i < h; i++) {
    for (j = 0; j < w; j++) {
      recon[i * recon_stride + j] = dst[i * dst_stride + j] +
                                    us_resi[i * us_resi_stride + j];
    }
  }

#if USE_POST_F
  // Do super-resolution post processing to the reconstruction
#if SR_USE_MULTI_F
  sr_post_filter(recon, recon_stride, enh_recon, enh_recon_stride, w, h,
                 f_hor, f_ver);
#else  // SR_USE_MULTI_F
  sr_post_filter(recon, recon_stride, enh_recon, enh_recon_stride, w, h);
#endif  // SR_USE_MULTI_F

  for (i = 0; i < h; i++)
    for (j = 0; j < w; j++)
      dst[i * dst_stride + j] = clip_pixel(enh_recon[i * enh_recon_stride + j]);
#else  // USE_POST_F
  for (i = 0; i < h; i++)
    for (j = 0; j < w; j++)
      dst[i * dst_stride + j] = clip_pixel(recon[i * recon_stride + j]);
#endif  // USE_POST_F
}
#endif  // CONFIG_SR_MODE
