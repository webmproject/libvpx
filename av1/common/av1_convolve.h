#ifndef AV1_COMMON_AV1_CONVOLVE_H_
#define AV1_COMMON_AV1_CONVOLVE_H_
#include "av1/common/filter.h"

#ifdef __cplusplus
extern "C" {
#endif

void av1_convolve(const uint8_t *src, int src_stride, uint8_t *dst,
                  int dst_stride, int w, int h,
#if CONFIG_DUAL_FILTER
                  const INTERP_FILTER *interp_filter,
#else
                  const INTERP_FILTER interp_filter,
#endif
                  const int subpel_x, int xstep, const int subpel_y, int ystep,
                  int avg);

#if CONFIG_AOM_HIGHBITDEPTH
void av1_highbd_convolve(const uint8_t *src, int src_stride, uint8_t *dst,
                         int dst_stride, int w, int h,
#if CONFIG_DUAL_FILTER
                         const INTERP_FILTER *interp_filter,
#else
                         const INTERP_FILTER interp_filter,
#endif
                         const int subpel_x, int xstep, const int subpel_y,
                         int ystep, int avg, int bd);
#endif  // CONFIG_AOM_HIGHBITDEPTH

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_AV1_CONVOLVE_H_
