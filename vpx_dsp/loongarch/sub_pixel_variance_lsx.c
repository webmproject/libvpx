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
#include "vpx_ports/mem.h"
#include "vpx_util/loongson_intrinsics.h"
#include "vpx_dsp/variance.h"

#define HADD_SW_S32(in0, in1)                  \
  do {                                         \
    __m128i res0_m;                            \
                                               \
    res0_m = __lsx_vhaddw_d_w(in0, in0);       \
    res0_m = __lsx_vhaddw_q_d(res0_m, res0_m); \
    in1 = __lsx_vpickve2gr_w(res0_m, 0);       \
  } while (0)

#define HORIZ_2TAP_FILT_UH(in0, in1, mask, coeff, shift, in2) \
  do {                                                        \
    __m128i tmp0_m, tmp1_m;                                   \
                                                              \
    tmp0_m = __lsx_vshuf_b(in1, in0, mask);                   \
    tmp1_m = __lsx_vdp2_h_bu(tmp0_m, coeff);                  \
    in2 = __lsx_vsrari_h(tmp1_m, shift);                      \
  } while (0)

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

static const uint8_t bilinear_filters_lsx[8][2] = {
  { 128, 0 }, { 112, 16 }, { 96, 32 }, { 80, 48 },
  { 64, 64 }, { 48, 80 },  { 32, 96 }, { 16, 112 },
};

#define VARIANCE_LARGE_WxH(sse, diff, shift) \
  (sse) - (((int64_t)(diff) * (diff)) >> (shift))

static uint32_t sub_pixel_sse_diff_16width_h_lsx(
    const uint8_t *src, int32_t src_stride, const uint8_t *dst,
    int32_t dst_stride, const uint8_t *filter, int32_t height, int32_t *diff) {
  uint32_t loop_cnt = (height >> 2);
  int32_t res;
  __m128i src0, src1, src2, src3, src4, src5, src6, src7;
  __m128i dst0, dst1, dst2, dst3, filt0;
  __m128i vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7;
  __m128i out0, out1, out2, out3, out4, out5, out6, out7;
  __m128i vec, var = __lsx_vldi(0);
  __m128i avg = var;
  __m128i mask = { 0x0403030202010100, 0x0807070606050504 };
  int32_t dst_stride2 = dst_stride << 1;
  int32_t dst_stride3 = dst_stride2 + dst_stride;
  int32_t dst_stride4 = dst_stride2 << 1;

  filt0 = __lsx_vldrepl_h(filter, 0);

  for (; loop_cnt--;) {
    DUP2_ARG2(__lsx_vld, src, 0, src, 8, src0, src1);
    src += src_stride;
    DUP2_ARG2(__lsx_vld, src, 0, src, 8, src2, src3);
    src += src_stride;
    DUP2_ARG2(__lsx_vld, src, 0, src, 8, src4, src5);
    src += src_stride;
    DUP2_ARG2(__lsx_vld, src, 0, src, 8, src6, src7);
    src += src_stride;

    dst0 = __lsx_vld(dst, 0);
    DUP2_ARG2(__lsx_vldx, dst, dst_stride, dst, dst_stride2, dst1, dst2);
    dst3 = __lsx_vldx(dst, dst_stride3);
    dst += dst_stride4;

    DUP2_ARG3(__lsx_vshuf_b, src0, src0, mask, src1, src1, mask, vec0, vec1);
    DUP2_ARG3(__lsx_vshuf_b, src2, src2, mask, src3, src3, mask, vec2, vec3);
    DUP2_ARG3(__lsx_vshuf_b, src4, src4, mask, src5, src5, mask, vec4, vec5);
    DUP2_ARG3(__lsx_vshuf_b, src6, src6, mask, src7, src7, mask, vec6, vec7);

    DUP4_ARG2(__lsx_vdp2_h_bu, vec0, filt0, vec1, filt0, vec2, filt0, vec3,
              filt0, out0, out1, out2, out3);
    DUP4_ARG2(__lsx_vdp2_h_bu, vec4, filt0, vec5, filt0, vec6, filt0, vec7,
              filt0, out4, out5, out6, out7);
    DUP4_ARG3(__lsx_vssrarni_bu_h, out1, out0, FILTER_BITS, out3, out2,
              FILTER_BITS, out5, out4, FILTER_BITS, out7, out6, FILTER_BITS,
              src0, src1, src2, src3);
    CALC_MSE_AVG_B(src0, dst0, var, avg);
    CALC_MSE_AVG_B(src1, dst1, var, avg);
    CALC_MSE_AVG_B(src2, dst2, var, avg);
    CALC_MSE_AVG_B(src3, dst3, var, avg);
  }
  vec = __lsx_vhaddw_w_h(avg, avg);
  HADD_SW_S32(vec, *diff);
  HADD_SW_S32(var, res);

  return res;
}

static uint32_t sub_pixel_sse_diff_32width_h_lsx(
    const uint8_t *src, int32_t src_stride, const uint8_t *dst,
    int32_t dst_stride, const uint8_t *filter, int32_t height, int32_t *diff) {
  uint32_t sse = 0;
  int32_t diff0[2];

  sse += sub_pixel_sse_diff_16width_h_lsx(src, src_stride, dst, dst_stride,
                                          filter, height, &diff0[0]);
  src += 16;
  dst += 16;

  sse += sub_pixel_sse_diff_16width_h_lsx(src, src_stride, dst, dst_stride,
                                          filter, height, &diff0[1]);

  *diff = diff0[0] + diff0[1];

  return sse;
}

static uint32_t sub_pixel_sse_diff_16width_v_lsx(
    const uint8_t *src, int32_t src_stride, const uint8_t *dst,
    int32_t dst_stride, const uint8_t *filter, int32_t height, int32_t *diff) {
  uint32_t loop_cnt = (height >> 2);
  int32_t res;
  __m128i ref0, ref1, ref2, ref3, src0, src1, src2, src3, src4;
  __m128i out0, out1, out2, out3, tmp0, tmp1, filt0, vec;
  __m128i vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7;
  __m128i var = __lsx_vldi(0);
  __m128i avg = var;
  int32_t src_stride2 = src_stride << 1;
  int32_t src_stride3 = src_stride2 + src_stride;
  int32_t src_stride4 = src_stride2 << 1;
  int32_t dst_stride2 = dst_stride << 1;
  int32_t dst_stride3 = dst_stride2 + dst_stride;
  int32_t dst_stride4 = dst_stride2 << 1;

  filt0 = __lsx_vldrepl_h(filter, 0);

  src0 = __lsx_vld(src, 0);
  src += src_stride;

  for (; loop_cnt--;) {
    src1 = __lsx_vld(src, 0);
    DUP2_ARG2(__lsx_vldx, src, src_stride, src, src_stride2, src2, src3);
    src4 = __lsx_vldx(src, src_stride3);
    src += src_stride4;
    ref0 = __lsx_vld(dst, 0);
    DUP2_ARG2(__lsx_vldx, dst, dst_stride, dst, dst_stride2, ref1, ref2);
    ref3 = __lsx_vldx(dst, dst_stride3);
    dst += dst_stride4;

    DUP2_ARG2(__lsx_vilvl_b, src1, src0, src2, src1, vec0, vec2);
    DUP2_ARG2(__lsx_vilvh_b, src1, src0, src2, src1, vec1, vec3);
    DUP2_ARG2(__lsx_vdp2_h_bu, vec0, filt0, vec1, filt0, tmp0, tmp1);
    out0 = __lsx_vssrarni_bu_h(tmp1, tmp0, FILTER_BITS);

    DUP2_ARG2(__lsx_vilvl_b, src3, src2, src4, src3, vec4, vec6);
    DUP2_ARG2(__lsx_vilvh_b, src3, src2, src4, src3, vec5, vec7);
    DUP2_ARG2(__lsx_vdp2_h_bu, vec2, filt0, vec3, filt0, tmp0, tmp1);
    out1 = __lsx_vssrarni_bu_h(tmp1, tmp0, FILTER_BITS);

    DUP2_ARG2(__lsx_vdp2_h_bu, vec4, filt0, vec5, filt0, tmp0, tmp1);
    out2 = __lsx_vssrarni_bu_h(tmp1, tmp0, FILTER_BITS);
    DUP2_ARG2(__lsx_vdp2_h_bu, vec6, filt0, vec7, filt0, tmp0, tmp1);
    out3 = __lsx_vssrarni_bu_h(tmp1, tmp0, FILTER_BITS);

    src0 = src4;

    CALC_MSE_AVG_B(out0, ref0, var, avg);
    CALC_MSE_AVG_B(out1, ref1, var, avg);
    CALC_MSE_AVG_B(out2, ref2, var, avg);
    CALC_MSE_AVG_B(out3, ref3, var, avg);
  }
  vec = __lsx_vhaddw_w_h(avg, avg);
  HADD_SW_S32(vec, *diff);
  HADD_SW_S32(var, res);

  return res;
}

static uint32_t sub_pixel_sse_diff_32width_v_lsx(
    const uint8_t *src, int32_t src_stride, const uint8_t *dst,
    int32_t dst_stride, const uint8_t *filter, int32_t height, int32_t *diff) {
  uint32_t sse = 0;
  int32_t diff0[2];

  sse += sub_pixel_sse_diff_16width_v_lsx(src, src_stride, dst, dst_stride,
                                          filter, height, &diff0[0]);
  src += 16;
  dst += 16;

  sse += sub_pixel_sse_diff_16width_v_lsx(src, src_stride, dst, dst_stride,
                                          filter, height, &diff0[1]);

  *diff = diff0[0] + diff0[1];

  return sse;
}

static uint32_t sub_pixel_sse_diff_16width_hv_lsx(
    const uint8_t *src, int32_t src_stride, const uint8_t *dst,
    int32_t dst_stride, const uint8_t *filter_horiz, const uint8_t *filter_vert,
    int32_t height, int32_t *diff) {
  uint32_t loop_cnt = (height >> 2);
  int32_t res;
  __m128i src0, src1, src2, src3, src4, src5, src6, src7;
  __m128i ref0, ref1, ref2, ref3, filt_hz, filt_vt, vec0, vec1;
  __m128i hz_out0, hz_out1, hz_out2, hz_out3, tmp0, tmp1, vec;
  __m128i var = __lsx_vldi(0);
  __m128i avg = var;
  __m128i mask = { 0x0403030202010100, 0x0807070606050504 };
  int32_t dst_stride2 = dst_stride << 1;
  int32_t dst_stride3 = dst_stride2 + dst_stride;
  int32_t dst_stride4 = dst_stride2 << 1;

  filt_hz = __lsx_vldrepl_h(filter_horiz, 0);
  filt_vt = __lsx_vldrepl_h(filter_vert, 0);

  DUP2_ARG2(__lsx_vld, src, 0, src, 8, src0, src1);
  src += src_stride;

  HORIZ_2TAP_FILT_UH(src0, src0, mask, filt_hz, FILTER_BITS, hz_out0);
  HORIZ_2TAP_FILT_UH(src1, src1, mask, filt_hz, FILTER_BITS, hz_out2);

  for (; loop_cnt--;) {
    DUP2_ARG2(__lsx_vld, src, 0, src, 8, src0, src1);
    src += src_stride;
    DUP2_ARG2(__lsx_vld, src, 0, src, 8, src2, src3);
    src += src_stride;
    DUP2_ARG2(__lsx_vld, src, 0, src, 8, src4, src5);
    src += src_stride;
    DUP2_ARG2(__lsx_vld, src, 0, src, 8, src6, src7);
    src += src_stride;

    ref0 = __lsx_vld(dst, 0);
    DUP2_ARG2(__lsx_vldx, dst, dst_stride, dst, dst_stride2, ref1, ref2);
    ref3 = __lsx_vldx(dst, dst_stride3);
    dst += dst_stride4;

    HORIZ_2TAP_FILT_UH(src0, src0, mask, filt_hz, FILTER_BITS, hz_out1);
    HORIZ_2TAP_FILT_UH(src1, src1, mask, filt_hz, FILTER_BITS, hz_out3);
    DUP2_ARG2(__lsx_vpackev_b, hz_out1, hz_out0, hz_out3, hz_out2, vec0, vec1);
    DUP2_ARG2(__lsx_vdp2_h_bu, vec0, filt_vt, vec1, filt_vt, tmp0, tmp1);
    src0 = __lsx_vssrarni_bu_h(tmp1, tmp0, FILTER_BITS);

    HORIZ_2TAP_FILT_UH(src2, src2, mask, filt_hz, FILTER_BITS, hz_out0);
    HORIZ_2TAP_FILT_UH(src3, src3, mask, filt_hz, FILTER_BITS, hz_out2);
    DUP2_ARG2(__lsx_vpackev_b, hz_out0, hz_out1, hz_out2, hz_out3, vec0, vec1);
    DUP2_ARG2(__lsx_vdp2_h_bu, vec0, filt_vt, vec1, filt_vt, tmp0, tmp1);
    src1 = __lsx_vssrarni_bu_h(tmp1, tmp0, FILTER_BITS);

    HORIZ_2TAP_FILT_UH(src4, src4, mask, filt_hz, FILTER_BITS, hz_out1);
    HORIZ_2TAP_FILT_UH(src5, src5, mask, filt_hz, FILTER_BITS, hz_out3);
    DUP2_ARG2(__lsx_vpackev_b, hz_out1, hz_out0, hz_out3, hz_out2, vec0, vec1);
    DUP2_ARG2(__lsx_vdp2_h_bu, vec0, filt_vt, vec1, filt_vt, tmp0, tmp1);
    src2 = __lsx_vssrarni_bu_h(tmp1, tmp0, FILTER_BITS);

    HORIZ_2TAP_FILT_UH(src6, src6, mask, filt_hz, FILTER_BITS, hz_out0);
    HORIZ_2TAP_FILT_UH(src7, src7, mask, filt_hz, FILTER_BITS, hz_out2);
    DUP2_ARG2(__lsx_vpackev_b, hz_out0, hz_out1, hz_out2, hz_out3, vec0, vec1);
    DUP2_ARG2(__lsx_vdp2_h_bu, vec0, filt_vt, vec1, filt_vt, tmp0, tmp1);
    src3 = __lsx_vssrarni_bu_h(tmp1, tmp0, FILTER_BITS);

    CALC_MSE_AVG_B(src0, ref0, var, avg);
    CALC_MSE_AVG_B(src1, ref1, var, avg);
    CALC_MSE_AVG_B(src2, ref2, var, avg);
    CALC_MSE_AVG_B(src3, ref3, var, avg);
  }
  vec = __lsx_vhaddw_w_h(avg, avg);
  HADD_SW_S32(vec, *diff);
  HADD_SW_S32(var, res);

  return res;
}

static uint32_t sub_pixel_sse_diff_32width_hv_lsx(
    const uint8_t *src, int32_t src_stride, const uint8_t *dst,
    int32_t dst_stride, const uint8_t *filter_horiz, const uint8_t *filter_vert,
    int32_t height, int32_t *diff) {
  uint32_t sse = 0;
  int32_t diff0[2];

  sse += sub_pixel_sse_diff_16width_hv_lsx(src, src_stride, dst, dst_stride,
                                           filter_horiz, filter_vert, height,
                                           &diff0[0]);
  src += 16;
  dst += 16;

  sse += sub_pixel_sse_diff_16width_hv_lsx(src, src_stride, dst, dst_stride,
                                           filter_horiz, filter_vert, height,
                                           &diff0[1]);

  *diff = diff0[0] + diff0[1];

  return sse;
}

#define VARIANCE_32Wx32H(sse, diff) VARIANCE_LARGE_WxH(sse, diff, 10)

#define VPX_SUB_PIXEL_VARIANCE_WDXHT_LSX(wd, ht)                              \
  uint32_t vpx_sub_pixel_variance##wd##x##ht##_lsx(                           \
      const uint8_t *src, int32_t src_stride, int32_t x_offset,               \
      int32_t y_offset, const uint8_t *ref, int32_t ref_stride,               \
      uint32_t *sse) {                                                        \
    int32_t diff;                                                             \
    uint32_t var;                                                             \
    const uint8_t *h_filter = bilinear_filters_lsx[x_offset];                 \
    const uint8_t *v_filter = bilinear_filters_lsx[y_offset];                 \
                                                                              \
    if (y_offset) {                                                           \
      if (x_offset) {                                                         \
        *sse = sub_pixel_sse_diff_##wd##width_hv_lsx(                         \
            src, src_stride, ref, ref_stride, h_filter, v_filter, ht, &diff); \
      } else {                                                                \
        *sse = sub_pixel_sse_diff_##wd##width_v_lsx(                          \
            src, src_stride, ref, ref_stride, v_filter, ht, &diff);           \
      }                                                                       \
                                                                              \
      var = VARIANCE_##wd##Wx##ht##H(*sse, diff);                             \
    } else {                                                                  \
      if (x_offset) {                                                         \
        *sse = sub_pixel_sse_diff_##wd##width_h_lsx(                          \
            src, src_stride, ref, ref_stride, h_filter, ht, &diff);           \
                                                                              \
        var = VARIANCE_##wd##Wx##ht##H(*sse, diff);                           \
      } else {                                                                \
        var = vpx_variance##wd##x##ht##_lsx(src, src_stride, ref, ref_stride, \
                                            sse);                             \
      }                                                                       \
    }                                                                         \
                                                                              \
    return var;                                                               \
  }

VPX_SUB_PIXEL_VARIANCE_WDXHT_LSX(32, 32)
