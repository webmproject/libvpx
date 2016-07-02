#include <assert.h>
#include <string.h>

#include "./vp10_rtcd.h"
#include "vp10/common/vp10_convolve.h"
#include "vp10/common/filter.h"
#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_ports/mem.h"

#define MAX_BLOCK_WIDTH (MAX_SB_SIZE)
#define MAX_BLOCK_HEIGHT (MAX_SB_SIZE)
#define MAX_STEP (32)
#define MAX_FILTER_TAP (12)

void vp10_convolve_horiz_c(const uint8_t *src, int src_stride, uint8_t *dst,
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

void vp10_convolve_vert_c(const uint8_t *src, int src_stride, uint8_t *dst,
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
#if CONFIG_DUAL_FILTER
                   const INTERP_FILTER *interp_filter,
#else
                   const INTERP_FILTER interp_filter,
#endif
                   const int subpel_x_q4, int x_step_q4, const int subpel_y_q4,
                   int y_step_q4, int ref_idx) {
  int ignore_horiz = x_step_q4 == 16 && subpel_x_q4 == 0;
  int ignore_vert = y_step_q4 == 16 && subpel_y_q4 == 0;

  assert(w <= MAX_BLOCK_WIDTH);
  assert(h <= MAX_BLOCK_HEIGHT);
  assert(y_step_q4 <= MAX_STEP);
  assert(x_step_q4 <= MAX_STEP);

  if (ignore_horiz && ignore_vert) {
    convolve_copy(src, src_stride, dst, dst_stride, w, h, ref_idx);
  } else if (ignore_vert) {
#if CONFIG_DUAL_FILTER
    InterpFilterParams filter_params =
        vp10_get_interp_filter_params(interp_filter[1 + 2 * ref_idx]);
#else
    InterpFilterParams filter_params =
        vp10_get_interp_filter_params(interp_filter);
#endif
    assert(filter_params.taps <= MAX_FILTER_TAP);
    vp10_convolve_horiz(src, src_stride, dst, dst_stride, w, h, filter_params,
                        subpel_x_q4, x_step_q4, ref_idx);
  } else if (ignore_horiz) {
#if CONFIG_DUAL_FILTER
    InterpFilterParams filter_params =
        vp10_get_interp_filter_params(interp_filter[2 * ref_idx]);
#else
    InterpFilterParams filter_params =
        vp10_get_interp_filter_params(interp_filter);
#endif
    assert(filter_params.taps <= MAX_FILTER_TAP);
    vp10_convolve_vert(src, src_stride, dst, dst_stride, w, h, filter_params,
                       subpel_y_q4, y_step_q4, ref_idx);
  } else {
    // temp's size is set to (maximum possible intermediate_height) *
    // MAX_BLOCK_WIDTH
    uint8_t temp[((((MAX_BLOCK_HEIGHT - 1) * MAX_STEP + 15) >> SUBPEL_BITS) +
                  MAX_FILTER_TAP) *
                 MAX_BLOCK_WIDTH];
    int temp_stride = MAX_BLOCK_WIDTH;
#if CONFIG_DUAL_FILTER
    InterpFilterParams filter_params_x =
        vp10_get_interp_filter_params(interp_filter[1 + 2 * ref_idx]);
    InterpFilterParams filter_params_y =
        vp10_get_interp_filter_params(interp_filter[0 + 2 * ref_idx]);
    InterpFilterParams filter_params = filter_params_x;

    // The filter size implies the required number of reference pixels for
    // the second stage filtering. It is possible that the two directions
    // require different filter sizes.
    int filter_size = filter_params_y.taps;
#else
    InterpFilterParams filter_params =
        vp10_get_interp_filter_params(interp_filter);
    int filter_size = filter_params.taps;
#endif
    int intermediate_height =
        (((h - 1) * y_step_q4 + subpel_y_q4) >> SUBPEL_BITS) + filter_size;

    assert(filter_params.taps <= MAX_FILTER_TAP);

    vp10_convolve_horiz(src - src_stride * (filter_size / 2 - 1), src_stride,
                        temp, temp_stride, w, intermediate_height,
                        filter_params, subpel_x_q4, x_step_q4, 0);

#if CONFIG_DUAL_FILTER
    filter_params = filter_params_y;
#else
    filter_params = vp10_get_interp_filter_params(interp_filter);
#endif
    filter_size = filter_params.taps;
    assert(filter_params.taps <= MAX_FILTER_TAP);

    vp10_convolve_vert(temp + temp_stride * (filter_size / 2 - 1), temp_stride,
                       dst, dst_stride, w, h, filter_params,
                       subpel_y_q4, y_step_q4, ref_idx);
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
#if CONFIG_DUAL_FILTER
                          const INTERP_FILTER *interp_filter,
#else
                          const INTERP_FILTER interp_filter,
#endif
                          const int subpel_x_q4, int x_step_q4,
                          const int subpel_y_q4, int y_step_q4, int ref_idx,
                          int bd) {
  uint16_t *src = CONVERT_TO_SHORTPTR(src8);
  uint16_t *dst = CONVERT_TO_SHORTPTR(dst8);
  int ignore_horiz = x_step_q4 == 16 && subpel_x_q4 == 0;
  int ignore_vert = y_step_q4 == 16 && subpel_y_q4 == 0;

  assert(w <= MAX_BLOCK_WIDTH);
  assert(h <= MAX_BLOCK_HEIGHT);
  assert(y_step_q4 <= MAX_STEP);
  assert(x_step_q4 <= MAX_STEP);

  if (ignore_horiz && ignore_vert) {
    highbd_convolve_copy(src, src_stride, dst, dst_stride, w, h, ref_idx, bd);
  } else if (ignore_vert) {
#if CONFIG_DUAL_FILTER
    InterpFilterParams filter_params =
        vp10_get_interp_filter_params(interp_filter[1 + 2 * ref_idx]);
#else
    InterpFilterParams filter_params =
        vp10_get_interp_filter_params(interp_filter);
#endif
    highbd_convolve_horiz(src, src_stride, dst, dst_stride, w, h, filter_params,
                          subpel_x_q4, x_step_q4, ref_idx, bd);
  } else if (ignore_horiz) {
#if CONFIG_DUAL_FILTER
    InterpFilterParams filter_params =
        vp10_get_interp_filter_params(interp_filter[0 + 2 * ref_idx]);
#else
    InterpFilterParams filter_params =
        vp10_get_interp_filter_params(interp_filter);
#endif
    highbd_convolve_vert(src, src_stride, dst, dst_stride, w, h, filter_params,
                         subpel_y_q4, y_step_q4, ref_idx, bd);
  } else {
    // temp's size is set to (maximum possible intermediate_height) *
    // MAX_BLOCK_WIDTH
    uint16_t temp[((((MAX_BLOCK_HEIGHT - 1) * MAX_STEP + 15) >> SUBPEL_BITS) +
                   MAX_FILTER_TAP) *
                  MAX_BLOCK_WIDTH];
    int temp_stride = MAX_BLOCK_WIDTH;

#if CONFIG_DUAL_FILTER
    InterpFilterParams filter_params_x =
        vp10_get_interp_filter_params(interp_filter[1 + 2 * ref_idx]);
    InterpFilterParams filter_params_y =
        vp10_get_interp_filter_params(interp_filter[0 + 2 * ref_idx]);
    InterpFilterParams filter_params = filter_params_x;
    int filter_size = filter_params_y.taps;
#else
    InterpFilterParams filter_params =
        vp10_get_interp_filter_params(interp_filter);
    int filter_size = filter_params.taps;
#endif

    int intermediate_height =
        (((h - 1) * y_step_q4 + subpel_y_q4) >> SUBPEL_BITS) + filter_size;

    highbd_convolve_horiz(src - src_stride * (filter_size / 2 - 1), src_stride,
                          temp, temp_stride, w, intermediate_height,
                          filter_params, subpel_x_q4, x_step_q4, 0, bd);

#if CONFIG_DUAL_FILTER
    filter_params = filter_params_y;
#endif
    filter_size = filter_params.taps;
    assert(filter_params.taps <= MAX_FILTER_TAP);

    highbd_convolve_vert(temp + temp_stride * (filter_size / 2 - 1),
                         temp_stride, dst, dst_stride, w, h, filter_params,
                         subpel_y_q4, y_step_q4, ref_idx, bd);
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
