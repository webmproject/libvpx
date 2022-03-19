/*
 *  Copyright (c) 2022 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_dsp_rtcd.h"
#include "vpx_util/loongson_intrinsics.h"

#define HADD_SW_S32(in)                        \
  ({                                           \
    __m128i res0_m;                            \
    int32_t sum_m;                             \
                                               \
    res0_m = __lsx_vhaddw_d_w(in, in);         \
    res0_m = __lsx_vhaddw_q_d(res0_m, res0_m); \
    sum_m = __lsx_vpickve2gr_w(res0_m, 0);     \
    sum_m;                                     \
  })

#define CALC_MSE_AVG_B(src, ref, var, sub)                                \
  {                                                                       \
    __m128i src_l0_m, src_l1_m;                                           \
    __m128i res_l0_m, res_l1_m;                                           \
                                                                          \
    src_l0_m = __lsx_vilvl_b(src, ref);                                   \
    src_l1_m = __lsx_vilvh_b(src, ref);                                   \
    DUP2_ARG2(__lsx_vhsubw_hu_bu, src_l0_m, src_l0_m, src_l1_m, src_l1_m, \
              res_l0_m, res_l1_m);                                        \
    var = __lsx_vdp2add_w_h(var, res_l0_m, res_l0_m);                     \
    var = __lsx_vdp2add_w_h(var, res_l1_m, res_l1_m);                     \
    sub = __lsx_vadd_h(sub, res_l0_m);                                    \
    sub = __lsx_vadd_h(sub, res_l1_m);                                    \
  }

#define VARIANCE_LARGE_WxH(sse, diff, shift) \
  (sse) - (((int64_t)(diff) * (diff)) >> (shift))

static uint32_t sse_diff_32width_lsx(const uint8_t *src_ptr, int32_t src_stride,
                                     const uint8_t *ref_ptr, int32_t ref_stride,
                                     int32_t height, int32_t *diff) {
  int32_t ht_cnt = (height >> 2);
  __m128i avg = __lsx_vldi(0);
  __m128i src0, src1, ref0, ref1;
  __m128i vec;
  __m128i var = avg;

  for (; ht_cnt--;) {
    DUP2_ARG2(__lsx_vld, src_ptr, 0, src_ptr, 16, src0, src1);
    src_ptr += src_stride;
    DUP2_ARG2(__lsx_vld, ref_ptr, 0, ref_ptr, 16, ref0, ref1);
    ref_ptr += ref_stride;
    CALC_MSE_AVG_B(src0, ref0, var, avg);
    CALC_MSE_AVG_B(src1, ref1, var, avg);

    DUP2_ARG2(__lsx_vld, src_ptr, 0, src_ptr, 16, src0, src1);
    src_ptr += src_stride;
    DUP2_ARG2(__lsx_vld, ref_ptr, 0, ref_ptr, 16, ref0, ref1);
    ref_ptr += ref_stride;
    CALC_MSE_AVG_B(src0, ref0, var, avg);
    CALC_MSE_AVG_B(src1, ref1, var, avg);

    DUP2_ARG2(__lsx_vld, src_ptr, 0, src_ptr, 16, src0, src1);
    src_ptr += src_stride;
    DUP2_ARG2(__lsx_vld, ref_ptr, 0, ref_ptr, 16, ref0, ref1);
    ref_ptr += ref_stride;
    CALC_MSE_AVG_B(src0, ref0, var, avg);
    CALC_MSE_AVG_B(src1, ref1, var, avg);

    DUP2_ARG2(__lsx_vld, src_ptr, 0, src_ptr, 16, src0, src1);
    src_ptr += src_stride;
    DUP2_ARG2(__lsx_vld, ref_ptr, 0, ref_ptr, 16, ref0, ref1);
    ref_ptr += ref_stride;
    CALC_MSE_AVG_B(src0, ref0, var, avg);
    CALC_MSE_AVG_B(src1, ref1, var, avg);
  }

  vec = __lsx_vhaddw_w_h(avg, avg);
  *diff = HADD_SW_S32(vec);

  return HADD_SW_S32(var);
}

static uint32_t sse_diff_64x64_lsx(const uint8_t *src_ptr, int32_t src_stride,
                                   const uint8_t *ref_ptr, int32_t ref_stride,
                                   int32_t *diff) {
  int32_t ht_cnt = 32;
  __m128i avg0 = __lsx_vldi(0);
  __m128i src0, src1, src2, src3;
  __m128i ref0, ref1, ref2, ref3;
  __m128i vec0, vec1;
  __m128i avg1 = avg0;
  __m128i avg2 = avg0;
  __m128i avg3 = avg0;
  __m128i var = avg0;

  for (; ht_cnt--;) {
    DUP4_ARG2(__lsx_vld, src_ptr, 0, src_ptr, 16, src_ptr, 32, src_ptr, 48,
              src0, src1, src2, src3);
    src_ptr += src_stride;
    DUP4_ARG2(__lsx_vld, ref_ptr, 0, ref_ptr, 16, ref_ptr, 32, ref_ptr, 48,
              ref0, ref1, ref2, ref3);
    ref_ptr += ref_stride;

    CALC_MSE_AVG_B(src0, ref0, var, avg0);
    CALC_MSE_AVG_B(src1, ref1, var, avg1);
    CALC_MSE_AVG_B(src2, ref2, var, avg2);
    CALC_MSE_AVG_B(src3, ref3, var, avg3);
    DUP4_ARG2(__lsx_vld, src_ptr, 0, src_ptr, 16, src_ptr, 32, src_ptr, 48,
              src0, src1, src2, src3);
    src_ptr += src_stride;
    DUP4_ARG2(__lsx_vld, ref_ptr, 0, ref_ptr, 16, ref_ptr, 32, ref_ptr, 48,
              ref0, ref1, ref2, ref3);
    ref_ptr += ref_stride;
    CALC_MSE_AVG_B(src0, ref0, var, avg0);
    CALC_MSE_AVG_B(src1, ref1, var, avg1);
    CALC_MSE_AVG_B(src2, ref2, var, avg2);
    CALC_MSE_AVG_B(src3, ref3, var, avg3);
  }
  vec0 = __lsx_vhaddw_w_h(avg0, avg0);
  vec1 = __lsx_vhaddw_w_h(avg1, avg1);
  vec0 = __lsx_vadd_w(vec0, vec1);
  vec1 = __lsx_vhaddw_w_h(avg2, avg2);
  vec0 = __lsx_vadd_w(vec0, vec1);
  vec1 = __lsx_vhaddw_w_h(avg3, avg3);
  vec0 = __lsx_vadd_w(vec0, vec1);
  *diff = HADD_SW_S32(vec0);

  return HADD_SW_S32(var);
}

#define VARIANCE_32Wx32H(sse, diff) VARIANCE_LARGE_WxH(sse, diff, 10);
#define VARIANCE_64Wx64H(sse, diff) VARIANCE_LARGE_WxH(sse, diff, 12);

#define VPX_VARIANCE_WDXHT_LSX(wd, ht)                                         \
  uint32_t vpx_variance##wd##x##ht##_lsx(                                      \
      const uint8_t *src, int32_t src_stride, const uint8_t *ref,              \
      int32_t ref_stride, uint32_t *sse) {                                     \
    int32_t diff;                                                              \
                                                                               \
    *sse =                                                                     \
        sse_diff_##wd##width_lsx(src, src_stride, ref, ref_stride, ht, &diff); \
                                                                               \
    return VARIANCE_##wd##Wx##ht##H(*sse, diff);                               \
  }

VPX_VARIANCE_WDXHT_LSX(32, 32)

uint32_t vpx_variance64x64_lsx(const uint8_t *src, int32_t src_stride,
                               const uint8_t *ref, int32_t ref_stride,
                               uint32_t *sse) {
  int32_t diff;

  *sse = sse_diff_64x64_lsx(src, src_stride, ref, ref_stride, &diff);

  return VARIANCE_64Wx64H(*sse, diff);
}
