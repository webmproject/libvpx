/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_dsp_rtcd.h"
#include "vpx_ports/mem.h"
#include "vpx_dsp/mips/macros_msa.h"
#include "vpx_dsp/variance.h"

static const uint8_t bilinear_filters_msa[8][2] = {
  { 128,   0, },
  { 112,  16, },
  {  96,  32, },
  {  80,  48, },
  {  64,  64, },
  {  48,  80, },
  {  32,  96, },
  {  16, 112, },
};

#define CALC_MSE_AVG_B(src, ref, var, sub) {                       \
  v16u8 src_l0_m, src_l1_m;                                        \
  v8i16 res_l0_m, res_l1_m;                                        \
                                                                   \
  ILVRL_B2_UB(src, ref, src_l0_m, src_l1_m);                       \
  HSUB_UB2_SH(src_l0_m, src_l1_m, res_l0_m, res_l1_m);             \
  DPADD_SH2_SW(res_l0_m, res_l1_m, res_l0_m, res_l1_m, var, var);  \
                                                                   \
  sub += res_l0_m + res_l1_m;                                      \
}

#define VARIANCE_WxH(sse, diff, shift) \
  sse - (((uint32_t)diff * diff) >> shift)

#define VARIANCE_LARGE_WxH(sse, diff, shift) \
  sse - (((int64_t)diff * diff) >> shift)

static uint32_t sub_pixel_sse_diff_4width_h_msa(const uint8_t *src,
                                                int32_t src_stride,
                                                const uint8_t *dst,
                                                int32_t dst_stride,
                                                const uint8_t *filter,
                                                int32_t height,
                                                int32_t *diff) {
  int16_t filtval;
  uint32_t loop_cnt;
  uint32_t ref0, ref1, ref2, ref3;
  v16u8 filt0, ref = { 0 };
  v16i8 src0, src1, src2, src3;
  v16i8 mask = { 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8 };
  v8u16 vec0, vec1, vec2, vec3;
  v8u16 const255;
  v8i16 avg = { 0 };
  v4i32 vec, var = { 0 };

  filtval = LH(filter);
  filt0 = (v16u8)__msa_fill_h(filtval);

  const255 = (v8u16)__msa_ldi_h(255);

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LD_SB4(src, src_stride, src0, src1, src2, src3);
    src += (4 * src_stride);
    LW4(dst, dst_stride, ref0, ref1, ref2, ref3);
    dst += (4 * dst_stride);
    INSERT_W4_UB(ref0, ref1, ref2, ref3, ref);
    VSHF_B2_UH(src0, src0, src1, src1, mask, mask, vec0, vec1);
    VSHF_B2_UH(src2, src2, src3, src3, mask, mask, vec2, vec3);
    DOTP_UB4_UH(vec0, vec1, vec2, vec3, filt0, filt0, filt0, filt0,
                vec0, vec1, vec2, vec3);
    SRARI_H4_UH(vec0, vec1, vec2, vec3, FILTER_BITS);
    MIN_UH4_UH(vec0, vec1, vec2, vec3, const255);
    PCKEV_B4_SB(vec0, vec0, vec1, vec1, vec2, vec2, vec3, vec3,
                src0, src1, src2, src3);
    ILVEV_W2_SB(src0, src1, src2, src3, src0, src2);
    src0 = (v16i8)__msa_ilvev_d((v2i64)src2, (v2i64)src0);
    CALC_MSE_AVG_B(src0, ref, var, avg);
  }

  vec = __msa_hadd_s_w(avg, avg);
  *diff = HADD_SW_S32(vec);

  return HADD_SW_S32(var);
}

static uint32_t sub_pixel_sse_diff_8width_h_msa(const uint8_t *src,
                                                int32_t src_stride,
                                                const uint8_t *dst,
                                                int32_t dst_stride,
                                                const uint8_t *filter,
                                                int32_t height,
                                                int32_t *diff) {
  int16_t filtval;
  uint32_t loop_cnt;
  v16u8 filt0, out, ref0, ref1, ref2, ref3;
  v16i8 src0, src1, src2, src3;
  v16i8 mask = { 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8 };
  v8u16 vec0, vec1, vec2, vec3, const255;
  v8i16 avg = { 0 };
  v4i32 vec, var = { 0 };

  filtval = LH(filter);
  filt0 = (v16u8)__msa_fill_h(filtval);

  const255 = (v8u16)__msa_ldi_h(255);

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LD_SB4(src, src_stride, src0, src1, src2, src3);
    src += (4 * src_stride);
    LD_UB4(dst, dst_stride, ref0, ref1, ref2, ref3);
    dst += (4 * dst_stride);

    PCKEV_D2_UB(ref1, ref0, ref3, ref2, ref0, ref1);
    VSHF_B2_UH(src0, src0, src1, src1, mask, mask, vec0, vec1);
    VSHF_B2_UH(src2, src2, src3, src3, mask, mask, vec2, vec3);
    DOTP_UB4_UH(vec0, vec1, vec2, vec3, filt0, filt0, filt0, filt0,
                vec0, vec1, vec2, vec3);
    SRARI_H4_UH(vec0, vec1, vec2, vec3, FILTER_BITS);
    MIN_UH4_UH(vec0, vec1, vec2, vec3, const255);
    PCKEV_B4_SB(vec0, vec0, vec1, vec1, vec2, vec2, vec3, vec3,
                src0, src1, src2, src3);
    out = (v16u8)__msa_ilvev_d((v2i64)src1, (v2i64)src0);
    CALC_MSE_AVG_B(out, ref0, var, avg);
    out = (v16u8)__msa_ilvev_d((v2i64)src3, (v2i64)src2);
    CALC_MSE_AVG_B(out, ref1, var, avg);
  }

  vec = __msa_hadd_s_w(avg, avg);
  *diff = HADD_SW_S32(vec);

  return HADD_SW_S32(var);
}

static uint32_t sub_pixel_sse_diff_16width_h_msa(const uint8_t *src,
                                                 int32_t src_stride,
                                                 const uint8_t *dst,
                                                 int32_t dst_stride,
                                                 const uint8_t *filter,
                                                 int32_t height,
                                                 int32_t *diff) {
  int16_t filtval;
  uint32_t loop_cnt;
  v16i8 src0, src1, src2, src3, src4, src5, src6, src7;
  v16i8 mask = { 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8 };
  v16u8 dst0, dst1, dst2, dst3, filt0;
  v8u16 vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7;
  v8u16 out0, out1, out2, out3, out4, out5, out6, out7;
  v8u16 const255;
  v8i16 avg = { 0 };
  v4i32 vec, var = { 0 };

  filtval = LH(filter);
  filt0 = (v16u8)__msa_fill_h(filtval);

  const255 = (v8u16)__msa_ldi_h(255);

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LD_SB4(src, src_stride, src0, src2, src4, src6);
    LD_SB4(src + 8, src_stride, src1, src3, src5, src7);
    src += (4 * src_stride);
    LD_UB4(dst, dst_stride, dst0, dst1, dst2, dst3);
    dst += (4 * dst_stride);

    VSHF_B2_UH(src0, src0, src1, src1, mask, mask, vec0, vec1);
    VSHF_B2_UH(src2, src2, src3, src3, mask, mask, vec2, vec3);
    VSHF_B2_UH(src4, src4, src5, src5, mask, mask, vec4, vec5);
    VSHF_B2_UH(src6, src6, src7, src7, mask, mask, vec6, vec7);
    DOTP_UB4_UH(vec0, vec1, vec2, vec3, filt0, filt0, filt0, filt0,
                out0, out1, out2, out3);
    DOTP_UB4_UH(vec4, vec5, vec6, vec7, filt0, filt0, filt0, filt0,
                out4, out5, out6, out7);
    SRARI_H4_UH(out0, out1, out2, out3, FILTER_BITS);
    SRARI_H4_UH(out4, out5, out6, out7, FILTER_BITS);
    MIN_UH4_UH(out0, out1, out2, out3, const255);
    MIN_UH4_UH(out4, out5, out6, out7, const255);
    PCKEV_B4_SB(out1, out0, out3, out2, out5, out4, out7, out6,
                src0, src1, src2, src3);
    CALC_MSE_AVG_B(src0, dst0, var, avg);
    CALC_MSE_AVG_B(src1, dst1, var, avg);
    CALC_MSE_AVG_B(src2, dst2, var, avg);
    CALC_MSE_AVG_B(src3, dst3, var, avg);
  }

  vec = __msa_hadd_s_w(avg, avg);
  *diff = HADD_SW_S32(vec);

  return HADD_SW_S32(var);
}

static uint32_t sub_pixel_sse_diff_32width_h_msa(const uint8_t *src,
                                                 int32_t src_stride,
                                                 const uint8_t *dst,
                                                 int32_t dst_stride,
                                                 const uint8_t *filter,
                                                 int32_t height,
                                                 int32_t *diff) {
  uint32_t loop_cnt, sse = 0;
  int32_t diff0[2];

  for (loop_cnt = 0; loop_cnt < 2; ++loop_cnt) {
    sse += sub_pixel_sse_diff_16width_h_msa(src, src_stride, dst, dst_stride,
                                            filter, height, &diff0[loop_cnt]);
    src += 16;
    dst += 16;
  }

  *diff = diff0[0] + diff0[1];

  return sse;
}

static uint32_t sub_pixel_sse_diff_64width_h_msa(const uint8_t *src,
                                                 int32_t src_stride,
                                                 const uint8_t *dst,
                                                 int32_t dst_stride,
                                                 const uint8_t *filter,
                                                 int32_t height,
                                                 int32_t *diff) {
  uint32_t loop_cnt, sse = 0;
  int32_t diff0[4];

  for (loop_cnt = 0; loop_cnt < 4; ++loop_cnt) {
    sse += sub_pixel_sse_diff_16width_h_msa(src, src_stride, dst, dst_stride,
                                            filter, height, &diff0[loop_cnt]);
    src += 16;
    dst += 16;
  }

  *diff = diff0[0] + diff0[1] + diff0[2] + diff0[3];

  return sse;
}

static uint32_t sub_pixel_sse_diff_4width_v_msa(const uint8_t *src,
                                                int32_t src_stride,
                                                const uint8_t *dst,
                                                int32_t dst_stride,
                                                const uint8_t *filter,
                                                int32_t height,
                                                int32_t *diff) {
  int16_t filtval;
  uint32_t loop_cnt;
  uint32_t ref0, ref1, ref2, ref3;
  v16u8 src0, src1, src2, src3, src4, out;
  v16u8 src10_r, src32_r, src21_r, src43_r;
  v16u8 ref = { 0 };
  v16u8 src2110, src4332;
  v16u8 filt0;
  v8i16 avg = { 0 };
  v4i32 vec, var = { 0 };
  v8u16 tmp0, tmp1;

  filtval = LH(filter);
  filt0 = (v16u8)__msa_fill_h(filtval);

  src0 = LD_UB(src);
  src += src_stride;

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LD_UB4(src, src_stride, src1, src2, src3, src4);
    src += (4 * src_stride);
    LW4(dst, dst_stride, ref0, ref1, ref2, ref3);
    dst += (4 * dst_stride);

    INSERT_W4_UB(ref0, ref1, ref2, ref3, ref);
    ILVR_B4_UB(src1, src0, src2, src1, src3, src2, src4, src3,
               src10_r, src21_r, src32_r, src43_r);
    ILVR_D2_UB(src21_r, src10_r, src43_r, src32_r, src2110, src4332);
    DOTP_UB2_UH(src2110, src4332, filt0, filt0, tmp0, tmp1);
    SRARI_H2_UH(tmp0, tmp1, FILTER_BITS);
    SAT_UH2_UH(tmp0, tmp1, 7);
    out = (v16u8)__msa_pckev_b((v16i8)tmp1, (v16i8)tmp0);
    CALC_MSE_AVG_B(out, ref, var, avg);
    src0 = src4;
  }

  vec = __msa_hadd_s_w(avg, avg);
  *diff = HADD_SW_S32(vec);

  return HADD_SW_S32(var);
}

static uint32_t sub_pixel_sse_diff_8width_v_msa(const uint8_t *src,
                                                int32_t src_stride,
                                                const uint8_t *dst,
                                                int32_t dst_stride,
                                                const uint8_t *filter,
                                                int32_t height,
                                                int32_t *diff) {
  int16_t filtval;
  uint32_t loop_cnt;
  v16u8 src0, src1, src2, src3, src4;
  v16u8 ref0, ref1, ref2, ref3;
  v8u16 vec0, vec1, vec2, vec3;
  v8u16 tmp0, tmp1, tmp2, tmp3;
  v16u8 filt0;
  v8i16 avg = { 0 };
  v4i32 vec, var = { 0 };

  filtval = LH(filter);
  filt0 = (v16u8)__msa_fill_h(filtval);

  src0 = LD_UB(src);
  src += src_stride;

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LD_UB4(src, src_stride, src1, src2, src3, src4);
    src += (4 * src_stride);
    LD_UB4(dst, dst_stride, ref0, ref1, ref2, ref3);
    dst += (4 * dst_stride);

    PCKEV_D2_UB(ref1, ref0, ref3, ref2, ref0, ref1);
    ILVR_B4_UH(src1, src0, src2, src1, src3, src2, src4, src3,
               vec0, vec1, vec2, vec3);
    DOTP_UB4_UH(vec0, vec1, vec2, vec3, filt0, filt0, filt0, filt0,
                tmp0, tmp1, tmp2, tmp3);
    SRARI_H4_UH(tmp0, tmp1, tmp2, tmp3, FILTER_BITS);
    SAT_UH4_UH(tmp0, tmp1, tmp2, tmp3, 7);
    PCKEV_B2_UB(tmp1, tmp0, tmp3, tmp2, src0, src1);
    CALC_MSE_AVG_B(src0, ref0, var, avg);
    CALC_MSE_AVG_B(src1, ref1, var, avg);
    src0 = src4;
  }

  vec = __msa_hadd_s_w(avg, avg);
  *diff = HADD_SW_S32(vec);

  return HADD_SW_S32(var);
}

static uint32_t sub_pixel_sse_diff_16width_v_msa(const uint8_t *src,
                                                 int32_t src_stride,
                                                 const uint8_t *dst,
                                                 int32_t dst_stride,
                                                 const uint8_t *filter,
                                                 int32_t height,
                                                 int32_t *diff) {
  int16_t filtval;
  uint32_t loop_cnt;
  v16u8 ref0, ref1, ref2, ref3;
  v16u8 src0, src1, src2, src3, src4;
  v16u8 out0, out1, out2, out3;
  v16u8 vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7;
  v8u16 tmp0, tmp1, tmp2, tmp3;
  v16u8 filt0;
  v8i16 avg = { 0 };
  v4i32 vec, var = { 0 };

  filtval = LH(filter);
  filt0 = (v16u8)__msa_fill_h(filtval);

  src0 = LD_UB(src);
  src += src_stride;

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LD_UB4(src, src_stride, src1, src2, src3, src4);
    src += (4 * src_stride);
    LD_UB4(dst, dst_stride, ref0, ref1, ref2, ref3);
    dst += (4 * dst_stride);

    ILVR_B2_UB(src1, src0, src2, src1, vec0, vec2);
    ILVL_B2_UB(src1, src0, src2, src1, vec1, vec3);
    DOTP_UB2_UH(vec0, vec1, filt0, filt0, tmp0, tmp1);
    SRARI_H2_UH(tmp0, tmp1, FILTER_BITS);
    SAT_UH2_UH(tmp0, tmp1, 7);
    out0 = (v16u8)__msa_pckev_b((v16i8)tmp1, (v16i8)tmp0);

    ILVR_B2_UB(src3, src2, src4, src3, vec4, vec6);
    ILVL_B2_UB(src3, src2, src4, src3, vec5, vec7);
    DOTP_UB2_UH(vec2, vec3, filt0, filt0, tmp2, tmp3);
    SRARI_H2_UH(tmp2, tmp3, FILTER_BITS);
    SAT_UH2_UH(tmp2, tmp3, 7);
    out1 = (v16u8)__msa_pckev_b((v16i8)tmp3, (v16i8)tmp2);

    DOTP_UB2_UH(vec4, vec5, filt0, filt0, tmp0, tmp1);
    SRARI_H2_UH(tmp0, tmp1, FILTER_BITS);
    SAT_UH2_UH(tmp0, tmp1, 7);
    out2 = (v16u8)__msa_pckev_b((v16i8)tmp1, (v16i8)tmp0);
    DOTP_UB2_UH(vec6, vec7, filt0, filt0, tmp2, tmp3);
    SRARI_H2_UH(tmp2, tmp3, FILTER_BITS);
    SAT_UH2_UH(tmp2, tmp3, 7);
    out3 = (v16u8)__msa_pckev_b((v16i8)tmp3, (v16i8)tmp2);

    src0 = src4;

    CALC_MSE_AVG_B(out0, ref0, var, avg);
    CALC_MSE_AVG_B(out1, ref1, var, avg);
    CALC_MSE_AVG_B(out2, ref2, var, avg);
    CALC_MSE_AVG_B(out3, ref3, var, avg);
  }

  vec = __msa_hadd_s_w(avg, avg);
  *diff = HADD_SW_S32(vec);

  return HADD_SW_S32(var);
}

static uint32_t sub_pixel_sse_diff_32width_v_msa(const uint8_t *src,
                                                 int32_t src_stride,
                                                 const uint8_t *dst,
                                                 int32_t dst_stride,
                                                 const uint8_t *filter,
                                                 int32_t height,
                                                 int32_t *diff) {
  uint32_t loop_cnt, sse = 0;
  int32_t diff0[2];

  for (loop_cnt = 0; loop_cnt < 2; ++loop_cnt) {
    sse += sub_pixel_sse_diff_16width_v_msa(src, src_stride, dst, dst_stride,
                                            filter, height, &diff0[loop_cnt]);
    src += 16;
    dst += 16;
  }

  *diff = diff0[0] + diff0[1];

  return sse;
}

static uint32_t sub_pixel_sse_diff_64width_v_msa(const uint8_t *src,
                                                 int32_t src_stride,
                                                 const uint8_t *dst,
                                                 int32_t dst_stride,
                                                 const uint8_t *filter,
                                                 int32_t height,
                                                 int32_t *diff) {
  uint32_t loop_cnt, sse = 0;
  int32_t diff0[4];

  for (loop_cnt = 0; loop_cnt < 4; ++loop_cnt) {
    sse += sub_pixel_sse_diff_16width_v_msa(src, src_stride, dst, dst_stride,
                                            filter, height, &diff0[loop_cnt]);
    src += 16;
    dst += 16;
  }

  *diff = diff0[0] + diff0[1] + diff0[2] + diff0[3];

  return sse;
}

static uint32_t sub_pixel_sse_diff_4width_hv_msa(const uint8_t *src,
                                                 int32_t src_stride,
                                                 const uint8_t *dst,
                                                 int32_t dst_stride,
                                                 const uint8_t *filter_horiz,
                                                 const uint8_t *filter_vert,
                                                 int32_t height,
                                                 int32_t *diff) {
  int16_t filtval;
  uint32_t loop_cnt;
  uint32_t ref0, ref1, ref2, ref3;
  v16u8 src0, src1, src2, src3, src4;
  v16u8 out, ref = { 0 };
  v16u8 filt_vt, filt_hz, vec0, vec1;
  v16u8 mask = { 0, 1, 1, 2, 2, 3, 3, 4, 16, 17, 17, 18, 18, 19, 19, 20 };
  v8u16 hz_out0, hz_out1, hz_out2, hz_out3, hz_out4;
  v8u16 tmp0, tmp1;
  v8i16 avg = { 0 };
  v4i32 vec, var = { 0 };

  filtval = LH(filter_horiz);
  filt_hz = (v16u8)__msa_fill_h(filtval);
  filtval = LH(filter_vert);
  filt_vt = (v16u8)__msa_fill_h(filtval);

  src0 = LD_UB(src);
  src += src_stride;

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LD_UB4(src, src_stride, src1, src2, src3, src4);
    src += (4 * src_stride);
    LW4(dst, dst_stride, ref0, ref1, ref2, ref3);
    dst += (4 * dst_stride);
    INSERT_W4_UB(ref0, ref1, ref2, ref3, ref);
    hz_out0 = HORIZ_2TAP_FILT_UH(src0, src1, mask, filt_hz, FILTER_BITS);
    hz_out2 = HORIZ_2TAP_FILT_UH(src2, src3, mask, filt_hz, FILTER_BITS);
    hz_out4 = HORIZ_2TAP_FILT_UH(src4, src4, mask, filt_hz, FILTER_BITS);
    hz_out1 = (v8u16)__msa_sldi_b((v16i8)hz_out2, (v16i8)hz_out0, 8);
    hz_out3 = (v8u16)__msa_pckod_d((v2i64)hz_out4, (v2i64)hz_out2);
    ILVEV_B2_UB(hz_out0, hz_out1, hz_out2, hz_out3, vec0, vec1);
    DOTP_UB2_UH(vec0, vec1, filt_vt, filt_vt, tmp0, tmp1);
    SRARI_H2_UH(tmp0, tmp1, FILTER_BITS);
    SAT_UH2_UH(tmp0, tmp1, 7);
    out = (v16u8)__msa_pckev_b((v16i8)tmp1, (v16i8)tmp0);
    CALC_MSE_AVG_B(out, ref, var, avg);
    src0 = src4;
  }

  vec = __msa_hadd_s_w(avg, avg);
  *diff = HADD_SW_S32(vec);

  return HADD_SW_S32(var);
}

static uint32_t sub_pixel_sse_diff_8width_hv_msa(const uint8_t *src,
                                                 int32_t src_stride,
                                                 const uint8_t *dst,
                                                 int32_t dst_stride,
                                                 const uint8_t *filter_horiz,
                                                 const uint8_t *filter_vert,
                                                 int32_t height,
                                                 int32_t *diff) {
  int16_t filtval;
  uint32_t loop_cnt;
  v16u8 ref0, ref1, ref2, ref3;
  v16u8 src0, src1, src2, src3, src4;
  v16u8 out0, out1;
  v16u8 mask = { 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8 };
  v8u16 hz_out0, hz_out1;
  v8u16 tmp0, tmp1, tmp2, tmp3;
  v16u8 filt_vt, filt_hz, vec0;
  v8i16 avg = { 0 };
  v4i32 vec, var = { 0 };

  filtval = LH(filter_horiz);
  filt_hz = (v16u8)__msa_fill_h(filtval);
  filtval = LH(filter_vert);
  filt_vt = (v16u8)__msa_fill_h(filtval);

  src0 = LD_UB(src);
  src += src_stride;
  hz_out0 = HORIZ_2TAP_FILT_UH(src0, src0, mask, filt_hz, FILTER_BITS);

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LD_UB4(src, src_stride, src1, src2, src3, src4);
    src += (4 * src_stride);
    LD_UB4(dst, dst_stride, ref0, ref1, ref2, ref3);
    dst += (4 * dst_stride);

    PCKEV_D2_UB(ref1, ref0, ref3, ref2, ref0, ref1);
    hz_out1 = HORIZ_2TAP_FILT_UH(src1, src1, mask, filt_hz, FILTER_BITS);
    vec0 = (v16u8)__msa_ilvev_b((v16i8)hz_out1, (v16i8)hz_out0);
    tmp0 = __msa_dotp_u_h(vec0, filt_vt);
    hz_out0 = HORIZ_2TAP_FILT_UH(src2, src2, mask, filt_hz, FILTER_BITS);
    vec0 = (v16u8)__msa_ilvev_b((v16i8)hz_out0, (v16i8)hz_out1);
    tmp1 = __msa_dotp_u_h(vec0, filt_vt);
    SRARI_H2_UH(tmp0, tmp1, FILTER_BITS);
    SAT_UH2_UH(tmp0, tmp1, 7);
    hz_out1 = HORIZ_2TAP_FILT_UH(src3, src3, mask, filt_hz, FILTER_BITS);
    vec0 = (v16u8)__msa_ilvev_b((v16i8)hz_out1, (v16i8)hz_out0);
    tmp2 = __msa_dotp_u_h(vec0, filt_vt);
    hz_out0 = HORIZ_2TAP_FILT_UH(src4, src4, mask, filt_hz, FILTER_BITS);
    vec0 = (v16u8)__msa_ilvev_b((v16i8)hz_out0, (v16i8)hz_out1);
    tmp3 = __msa_dotp_u_h(vec0, filt_vt);
    SRARI_H2_UH(tmp2, tmp3, FILTER_BITS);
    SAT_UH2_UH(tmp2, tmp3, 7);
    PCKEV_B2_UB(tmp1, tmp0, tmp3, tmp2, out0, out1);
    CALC_MSE_AVG_B(out0, ref0, var, avg);
    CALC_MSE_AVG_B(out1, ref1, var, avg);
  }

  vec = __msa_hadd_s_w(avg, avg);
  *diff = HADD_SW_S32(vec);

  return HADD_SW_S32(var);
}

static uint32_t sub_pixel_sse_diff_16width_hv_msa(const uint8_t *src,
                                                  int32_t src_stride,
                                                  const uint8_t *dst,
                                                  int32_t dst_stride,
                                                  const uint8_t *filter_horiz,
                                                  const uint8_t *filter_vert,
                                                  int32_t height,
                                                  int32_t *diff) {
  int16_t filtval;
  uint32_t loop_cnt;
  v16u8 src0, src1, src2, src3, src4, src5, src6, src7;
  v16u8 ref0, ref1, ref2, ref3;
  v16u8 filt_hz, filt_vt, vec0, vec1;
  v16u8 mask = { 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8 };
  v8u16 hz_out0, hz_out1, hz_out2, hz_out3;
  v8u16 tmp0, tmp1;
  v8i16 avg = { 0 };
  v4i32 vec, var = { 0 };

  filtval = LH(filter_horiz);
  filt_hz = (v16u8)__msa_fill_h(filtval);
  filtval = LH(filter_vert);
  filt_vt = (v16u8)__msa_fill_h(filtval);

  LD_UB2(src, 8, src0, src1);
  src += src_stride;

  hz_out0 = HORIZ_2TAP_FILT_UH(src0, src0, mask, filt_hz, FILTER_BITS);
  hz_out2 = HORIZ_2TAP_FILT_UH(src1, src1, mask, filt_hz, FILTER_BITS);

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LD_UB4(src, src_stride, src0, src2, src4, src6);
    LD_UB4(src + 8, src_stride, src1, src3, src5, src7);
    src += (4 * src_stride);
    LD_UB4(dst, dst_stride, ref0, ref1, ref2, ref3);
    dst += (4 * dst_stride);

    hz_out1 = HORIZ_2TAP_FILT_UH(src0, src0, mask, filt_hz, FILTER_BITS);
    hz_out3 = HORIZ_2TAP_FILT_UH(src1, src1, mask, filt_hz, FILTER_BITS);
    ILVEV_B2_UB(hz_out0, hz_out1, hz_out2, hz_out3, vec0, vec1);
    DOTP_UB2_UH(vec0, vec1, filt_vt, filt_vt, tmp0, tmp1);
    SRARI_H2_UH(tmp0, tmp1, FILTER_BITS);
    SAT_UH2_UH(tmp0, tmp1, 7);
    src0 = (v16u8)__msa_pckev_b((v16i8)tmp1, (v16i8)tmp0);

    hz_out0 = HORIZ_2TAP_FILT_UH(src2, src2, mask, filt_hz, FILTER_BITS);
    hz_out2 = HORIZ_2TAP_FILT_UH(src3, src3, mask, filt_hz, FILTER_BITS);
    ILVEV_B2_UB(hz_out1, hz_out0, hz_out3, hz_out2, vec0, vec1);
    DOTP_UB2_UH(vec0, vec1, filt_vt, filt_vt, tmp0, tmp1);
    SRARI_H2_UH(tmp0, tmp1, FILTER_BITS);
    SAT_UH2_UH(tmp0, tmp1, 7);
    src1 = (v16u8)__msa_pckev_b((v16i8)tmp1, (v16i8)tmp0);

    hz_out1 = HORIZ_2TAP_FILT_UH(src4, src4, mask, filt_hz, FILTER_BITS);
    hz_out3 = HORIZ_2TAP_FILT_UH(src5, src5, mask, filt_hz, FILTER_BITS);
    ILVEV_B2_UB(hz_out0, hz_out1, hz_out2, hz_out3, vec0, vec1);
    DOTP_UB2_UH(vec0, vec1, filt_vt, filt_vt, tmp0, tmp1);
    SRARI_H2_UH(tmp0, tmp1, FILTER_BITS);
    SAT_UH2_UH(tmp0, tmp1, 7);
    src2 = (v16u8)__msa_pckev_b((v16i8)tmp1, (v16i8)tmp0);

    hz_out0 = HORIZ_2TAP_FILT_UH(src6, src6, mask, filt_hz, FILTER_BITS);
    hz_out2 = HORIZ_2TAP_FILT_UH(src7, src7, mask, filt_hz, FILTER_BITS);
    ILVEV_B2_UB(hz_out1, hz_out0, hz_out3, hz_out2, vec0, vec1);
    DOTP_UB2_UH(vec0, vec1, filt_vt, filt_vt, tmp0, tmp1);
    SRARI_H2_UH(tmp0, tmp1, FILTER_BITS);
    SAT_UH2_UH(tmp0, tmp1, 7);
    src3 = (v16u8)__msa_pckev_b((v16i8)tmp1, (v16i8)tmp0);

    CALC_MSE_AVG_B(src0, ref0, var, avg);
    CALC_MSE_AVG_B(src1, ref1, var, avg);
    CALC_MSE_AVG_B(src2, ref2, var, avg);
    CALC_MSE_AVG_B(src3, ref3, var, avg);
  }

  vec = __msa_hadd_s_w(avg, avg);
  *diff = HADD_SW_S32(vec);

  return HADD_SW_S32(var);
}

static uint32_t sub_pixel_sse_diff_32width_hv_msa(const uint8_t *src,
                                                  int32_t src_stride,
                                                  const uint8_t *dst,
                                                  int32_t dst_stride,
                                                  const uint8_t *filter_horiz,
                                                  const uint8_t *filter_vert,
                                                  int32_t height,
                                                  int32_t *diff) {
  uint32_t loop_cnt, sse = 0;
  int32_t diff0[2];

  for (loop_cnt = 0; loop_cnt < 2; ++loop_cnt) {
    sse += sub_pixel_sse_diff_16width_hv_msa(src, src_stride, dst, dst_stride,
                                             filter_horiz, filter_vert, height,
                                             &diff0[loop_cnt]);
    src += 16;
    dst += 16;
  }

  *diff = diff0[0] + diff0[1];

  return sse;
}

static uint32_t sub_pixel_sse_diff_64width_hv_msa(const uint8_t *src,
                                                  int32_t src_stride,
                                                  const uint8_t *dst,
                                                  int32_t dst_stride,
                                                  const uint8_t *filter_horiz,
                                                  const uint8_t *filter_vert,
                                                  int32_t height,
                                                  int32_t *diff) {
  uint32_t loop_cnt, sse = 0;
  int32_t diff0[4];

  for (loop_cnt = 0; loop_cnt < 4; ++loop_cnt) {
    sse += sub_pixel_sse_diff_16width_hv_msa(src, src_stride, dst, dst_stride,
                                             filter_horiz, filter_vert, height,
                                             &diff0[loop_cnt]);
    src += 16;
    dst += 16;
  }

  *diff = diff0[0] + diff0[1] + diff0[2] + diff0[3];

  return sse;
}

#define VARIANCE_4Wx4H(sse, diff) VARIANCE_WxH(sse, diff, 4);
#define VARIANCE_4Wx8H(sse, diff) VARIANCE_WxH(sse, diff, 5);
#define VARIANCE_8Wx4H(sse, diff) VARIANCE_WxH(sse, diff, 5);
#define VARIANCE_8Wx8H(sse, diff) VARIANCE_WxH(sse, diff, 6);
#define VARIANCE_8Wx16H(sse, diff) VARIANCE_WxH(sse, diff, 7);
#define VARIANCE_16Wx8H(sse, diff) VARIANCE_WxH(sse, diff, 7);
#define VARIANCE_16Wx16H(sse, diff) VARIANCE_WxH(sse, diff, 8);

#define VARIANCE_16Wx32H(sse, diff) VARIANCE_LARGE_WxH(sse, diff, 9);
#define VARIANCE_32Wx16H(sse, diff) VARIANCE_LARGE_WxH(sse, diff, 9);
#define VARIANCE_32Wx32H(sse, diff) VARIANCE_LARGE_WxH(sse, diff, 10);
#define VARIANCE_32Wx64H(sse, diff) VARIANCE_LARGE_WxH(sse, diff, 11);
#define VARIANCE_64Wx32H(sse, diff) VARIANCE_LARGE_WxH(sse, diff, 11);
#define VARIANCE_64Wx64H(sse, diff) VARIANCE_LARGE_WxH(sse, diff, 12);

#define VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(wd, ht)                         \
uint32_t vpx_sub_pixel_variance##wd##x##ht##_msa(const uint8_t *src,     \
                                                 int32_t src_stride,     \
                                                 int32_t xoffset,        \
                                                 int32_t yoffset,        \
                                                 const uint8_t *ref,     \
                                                 int32_t ref_stride,     \
                                                 uint32_t *sse) {        \
  int32_t diff;                                                          \
  uint32_t var;                                                          \
  const uint8_t *h_filter = bilinear_filters_msa[xoffset];               \
  const uint8_t *v_filter = bilinear_filters_msa[yoffset];               \
                                                                         \
  if (yoffset) {                                                         \
    if (xoffset) {                                                       \
      *sse = sub_pixel_sse_diff_##wd##width_hv_msa(src, src_stride,      \
                                                   ref, ref_stride,      \
                                                   h_filter, v_filter,   \
                                                   ht, &diff);           \
    } else {                                                             \
      *sse = sub_pixel_sse_diff_##wd##width_v_msa(src, src_stride,       \
                                                  ref, ref_stride,       \
                                                  v_filter, ht, &diff);  \
    }                                                                    \
                                                                         \
    var = VARIANCE_##wd##Wx##ht##H(*sse, diff);                          \
  } else {                                                               \
    if (xoffset) {                                                       \
      *sse = sub_pixel_sse_diff_##wd##width_h_msa(src, src_stride,       \
                                                  ref, ref_stride,       \
                                                  h_filter, ht, &diff);  \
                                                                         \
      var = VARIANCE_##wd##Wx##ht##H(*sse, diff);                        \
    } else {                                                             \
      var = vpx_variance##wd##x##ht##_msa(src, src_stride,               \
                                          ref, ref_stride, sse);         \
    }                                                                    \
  }                                                                      \
                                                                         \
  return var;                                                            \
}

VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(4, 4);
VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(4, 8);

VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(8, 4);
VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(8, 8);
VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(8, 16);

VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(16, 8);
VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(16, 16);
VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(16, 32);

VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(32, 16);
VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(32, 32);
VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(32, 64);

VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(64, 32);
VPX_SUB_PIXEL_VARIANCE_WDXHT_MSA(64, 64);
