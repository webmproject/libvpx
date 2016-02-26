#include <assert.h>
#include <string.h>

#include "vp10/common/filter.h"
#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_ports/mem.h"

#define MAX_BLOCK_WIDTH (64)
#define MAX_BLOCK_HEIGHT (64)
#define MAX_STEP (32)
#define MAX_FILTER_TAP (12)

static void convolve_horiz(const uint8_t *src, int src_stride, uint8_t *dst,
                           int dst_stride, int w, int h,
                           const InterpFilterParams filter_params,
                           const int subpel_x_q4, int x_step_q4, int avg) {
  int x, y;
  int filter_size = filter_params.taps;
  src -= filter_size / 2 - 1;
  for (y = 0; y < h; ++y) {
    int x_q4 = subpel_x_q4;
    for (x = 0; x < w; ++x) {
      const uint8_t *const src_x = &src[x_q4 >> SUBPEL_BITS];
      const int16_t *x_filter =
          vp10_get_interp_filter_subpel_kernel(
              filter_params, x_q4 & SUBPEL_MASK);
      int k, sum = 0;
      for (k = 0; k < filter_size; ++k) sum += src_x[k] * x_filter[k];
      if (avg) {
        dst[x] = ROUND_POWER_OF_TWO(
            dst[x] + clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS)), 1);
      } else {
        dst[x] = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
      }
      x_q4 += x_step_q4;
    }
    src += src_stride;
    dst += dst_stride;
  }
}

static void convolve_vert(const uint8_t *src, int src_stride, uint8_t *dst,
                          int dst_stride, int w, int h,
                          const InterpFilterParams filter_params,
                          const int subpel_y_q4, int y_step_q4, int avg) {
  int x, y;
  int filter_size = filter_params.taps;
  src -= src_stride * (filter_size / 2 - 1);

  for (x = 0; x < w; ++x) {
    int y_q4 = subpel_y_q4;
    for (y = 0; y < h; ++y) {
      const uint8_t *const src_y = &src[(y_q4 >> SUBPEL_BITS) * src_stride];
      const int16_t *y_filter =
          vp10_get_interp_filter_subpel_kernel(
              filter_params, y_q4 & SUBPEL_MASK);
      int k, sum = 0;
      for (k = 0; k < filter_size; ++k)
        sum += src_y[k * src_stride] * y_filter[k];
      if (avg) {
        dst[y * dst_stride] = ROUND_POWER_OF_TWO(
            dst[y * dst_stride] +
                clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS)),
            1);
      } else {
        dst[y * dst_stride] = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
      }
      y_q4 += y_step_q4;
    }
    ++src;
    ++dst;
  }
}

static void convolve_copy(const uint8_t *src, int src_stride, uint8_t *dst,
                          int dst_stride, int w, int h, int avg) {
  if (avg == 0) {
    int r;
    for (r = 0; r < h; ++r) {
      memcpy(dst, src, w);
      src += src_stride;
      dst += dst_stride;
    }
  } else {
    int r, c;
    for (r = 0; r < h; ++r) {
      for (c = 0; c < w; ++c) {
        dst[c] = clip_pixel(ROUND_POWER_OF_TWO(dst[c] + src[c], 1));
      }
      src += src_stride;
      dst += dst_stride;
    }
  }
}

void vp10_convolve(const uint8_t *src, int src_stride, uint8_t *dst,
                   int dst_stride, int w, int h,
                   const InterpFilterParams filter_params,
                   const int subpel_x_q4, int x_step_q4, const int subpel_y_q4,
                   int y_step_q4, int avg) {
  int filter_size = filter_params.taps;
  int ignore_horiz = x_step_q4 == 16 && subpel_x_q4 == 0;
  int ignore_vert = y_step_q4 == 16 && subpel_y_q4 == 0;

  assert(w <= MAX_BLOCK_WIDTH);
  assert(h <= MAX_BLOCK_HEIGHT);
  assert(y_step_q4 <= MAX_STEP);
  assert(x_step_q4 <= MAX_STEP);
  assert(filter_params.taps <= MAX_FILTER_TAP);

  if (ignore_horiz && ignore_vert) {
    convolve_copy(src, src_stride, dst, dst_stride, w, h, avg);
  } else if (ignore_vert) {
    convolve_horiz(src, src_stride, dst, dst_stride, w, h, filter_params,
                   subpel_x_q4, x_step_q4, avg);
  } else if (ignore_horiz) {
    convolve_vert(src, src_stride, dst, dst_stride, w, h, filter_params,
                  subpel_y_q4, y_step_q4, avg);
  } else {
    // temp's size is set to (maximum possible intermediate_height) *
    // MAX_BLOCK_WIDTH
    uint8_t temp[((((MAX_BLOCK_HEIGHT - 1) * MAX_STEP + 15) >> SUBPEL_BITS) +
                  MAX_FILTER_TAP) *
                 MAX_BLOCK_WIDTH];
    int temp_stride = MAX_BLOCK_WIDTH;

    int intermediate_height =
        (((h - 1) * y_step_q4 + subpel_y_q4) >> SUBPEL_BITS) + filter_size;

    convolve_horiz(src - src_stride * (filter_size / 2 - 1), src_stride, temp,
                   temp_stride, w, intermediate_height, filter_params,
                   subpel_x_q4, x_step_q4, 0);
    convolve_vert(temp + temp_stride * (filter_size / 2 - 1), temp_stride, dst,
                  dst_stride, w, h, filter_params, subpel_y_q4, y_step_q4, avg);
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
static void highbd_convolve_horiz(const uint16_t *src, int src_stride,
                                  uint16_t *dst, int dst_stride, int w, int h,
                                  const InterpFilterParams filter_params,
                                  const int subpel_x_q4, int x_step_q4, int avg,
                                  int bd) {
  int x, y;
  int filter_size = filter_params.taps;
  src -= filter_size / 2 - 1;
  for (y = 0; y < h; ++y) {
    int x_q4 = subpel_x_q4;
    for (x = 0; x < w; ++x) {
      const uint16_t *const src_x = &src[x_q4 >> SUBPEL_BITS];
      const int16_t *x_filter =
          vp10_get_interp_filter_subpel_kernel(
              filter_params, x_q4 & SUBPEL_MASK);
      int k, sum = 0;
      for (k = 0; k < filter_size; ++k) sum += src_x[k] * x_filter[k];
      if (avg)
        dst[x] = ROUND_POWER_OF_TWO(
            dst[x] +
                clip_pixel_highbd(ROUND_POWER_OF_TWO(sum, FILTER_BITS), bd),
            1);
      else
        dst[x] = clip_pixel_highbd(ROUND_POWER_OF_TWO(sum, FILTER_BITS), bd);
      x_q4 += x_step_q4;
    }
    src += src_stride;
    dst += dst_stride;
  }
}

static void highbd_convolve_vert(const uint16_t *src, int src_stride,
                                 uint16_t *dst, int dst_stride, int w, int h,
                                 const InterpFilterParams filter_params,
                                 const int subpel_y_q4, int y_step_q4, int avg,
                                 int bd) {
  int x, y;
  int filter_size = filter_params.taps;
  src -= src_stride * (filter_size / 2 - 1);

  for (x = 0; x < w; ++x) {
    int y_q4 = subpel_y_q4;
    for (y = 0; y < h; ++y) {
      const uint16_t *const src_y = &src[(y_q4 >> SUBPEL_BITS) * src_stride];
      const int16_t *y_filter =
          vp10_get_interp_filter_subpel_kernel(
              filter_params, y_q4 & SUBPEL_MASK);
      int k, sum = 0;
      for (k = 0; k < filter_size; ++k)
        sum += src_y[k * src_stride] * y_filter[k];
      if (avg) {
        dst[y * dst_stride] = ROUND_POWER_OF_TWO(
            dst[y * dst_stride] +
                clip_pixel_highbd(ROUND_POWER_OF_TWO(sum, FILTER_BITS), bd),
            1);
      } else {
        dst[y * dst_stride] =
            clip_pixel_highbd(ROUND_POWER_OF_TWO(sum, FILTER_BITS), bd);
      }
      y_q4 += y_step_q4;
    }
    ++src;
    ++dst;
  }
}

static void highbd_convolve_copy(const uint16_t *src, int src_stride,
                                 uint16_t *dst, int dst_stride, int w, int h,
                                 int avg, int bd) {
  if (avg == 0) {
    int r;
    for (r = 0; r < h; ++r) {
      memcpy(dst, src, w * sizeof(*src));
      src += src_stride;
      dst += dst_stride;
    }
  } else {
    int r, c;
    for (r = 0; r < h; ++r) {
      for (c = 0; c < w; ++c) {
        dst[c] = clip_pixel_highbd(ROUND_POWER_OF_TWO(dst[c] + src[c], 1), bd);
      }
      src += src_stride;
      dst += dst_stride;
    }
  }
}

void vp10_highbd_convolve(const uint8_t *src8, int src_stride, uint8_t *dst8,
                          int dst_stride, int w, int h,
                          const InterpFilterParams filter_params,
                          const int subpel_x_q4, int x_step_q4,
                          const int subpel_y_q4, int y_step_q4, int avg,
                          int bd) {
  uint16_t *src = CONVERT_TO_SHORTPTR(src8);
  uint16_t *dst = CONVERT_TO_SHORTPTR(dst8);
  int filter_size = filter_params.taps;
  int ignore_horiz = x_step_q4 == 16 && subpel_x_q4 == 0;
  int ignore_vert = y_step_q4 == 16 && subpel_y_q4 == 0;

  assert(w <= MAX_BLOCK_WIDTH);
  assert(h <= MAX_BLOCK_HEIGHT);
  assert(y_step_q4 <= MAX_STEP);
  assert(x_step_q4 <= MAX_STEP);
  assert(filter_params.taps <= MAX_FILTER_TAP);

  if (ignore_horiz && ignore_vert) {
    highbd_convolve_copy(src, src_stride, dst, dst_stride, w, h, avg, bd);
  } else if (ignore_vert) {
    highbd_convolve_horiz(src, src_stride, dst, dst_stride, w, h, filter_params,
                          subpel_x_q4, x_step_q4, avg, bd);
  } else if (ignore_horiz) {
    highbd_convolve_vert(src, src_stride, dst, dst_stride, w, h, filter_params,
                         subpel_y_q4, y_step_q4, avg, bd);
  } else {
    // temp's size is set to (maximum possible intermediate_height) *
    // MAX_BLOCK_WIDTH
    uint16_t temp[((((MAX_BLOCK_HEIGHT - 1) * MAX_STEP + 15) >> SUBPEL_BITS) +
                   MAX_FILTER_TAP) *
                  MAX_BLOCK_WIDTH];
    int temp_stride = MAX_BLOCK_WIDTH;

    int intermediate_height =
        (((h - 1) * y_step_q4 + subpel_y_q4) >> SUBPEL_BITS) + filter_size;

    highbd_convolve_horiz(src - src_stride * (filter_size / 2 - 1), src_stride,
                          temp, temp_stride, w, intermediate_height,
                          filter_params, subpel_x_q4, x_step_q4, 0, bd);
    highbd_convolve_vert(temp + temp_stride * (filter_size / 2 - 1),
                         temp_stride, dst, dst_stride, w, h, filter_params,
                         subpel_y_q4, y_step_q4, avg, bd);
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
