#ifndef VP10_COMMON_VP10_CONVOLVE_H_
#define VP10_COMMON_VP10_CONVOLVE_H_
#include "vp10/common/filter.h"

#ifdef __cplusplus
extern "C" {
#endif

void vp10_convolve(const uint8_t *src, int src_stride,
                   uint8_t *dst, int dst_stride,
                   int w, int h,
#if CONFIG_DUAL_FILTER
                   const INTERP_FILTER *interp_filter,
#else
                   const INTERP_FILTER interp_filter,
#endif
                   const int subpel_x,
                   const int subpel_y,
                   int xstep, int ystep, int avg);

#if CONFIG_VP9_HIGHBITDEPTH
void vp10_highbd_convolve(const uint8_t *src, int src_stride,
                   uint8_t *dst, int dst_stride,
                   int w, int h,
                   const InterpFilterParams filter_params,
                   const int subpel_x,
                   const int subpel_y,
                   int xstep, int ystep, int avg, int bd);
#endif  // CONFIG_VP9_HIGHBITDEPTH

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_COMMON_VP10_CONVOLVE_H_
