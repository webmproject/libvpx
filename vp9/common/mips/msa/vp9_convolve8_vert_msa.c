/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vp9_rtcd.h"
#include "vp9/common/mips/msa/vp9_convolve_msa.h"

static void common_vt_8t_4w_msa(const uint8_t *src, int32_t src_stride,
                                uint8_t *dst, int32_t dst_stride,
                                int8_t *filter, int32_t height) {
  uint32_t loop_cnt;
  v16i8 src0, src1, src2, src3, src4, src5, src6, src7, src8, src9, src10;
  v16i8 src10_r, src32_r, src54_r, src76_r, src98_r;
  v16i8 src21_r, src43_r, src65_r, src87_r, src109_r;
  v16i8 src2110, src4332, src6554, src8776, src10998;
  v8i16 filt, out10, out32;
  v16i8 filt0, filt1, filt2, filt3;

  src -= (3 * src_stride);

  filt = LOAD_SH(filter);
  filt0 = (v16i8)__msa_splati_h(filt, 0);
  filt1 = (v16i8)__msa_splati_h(filt, 1);
  filt2 = (v16i8)__msa_splati_h(filt, 2);
  filt3 = (v16i8)__msa_splati_h(filt, 3);

  LOAD_7VECS_SB(src, src_stride, src0, src1, src2, src3, src4, src5, src6);
  src += (7 * src_stride);

  ILVR_B_6VECS_SB(src0, src2, src4, src1, src3, src5,
                  src1, src3, src5, src2, src4, src6,
                  src10_r, src32_r, src54_r, src21_r, src43_r, src65_r);

  ILVR_D_3VECS_SB(src2110, src21_r, src10_r, src4332, src43_r, src32_r,
                  src6554, src65_r, src54_r);

  XORI_B_3VECS_SB(src2110, src4332, src6554, src2110, src4332, src6554, 128);

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LOAD_4VECS_SB(src, src_stride, src7, src8, src9, src10);
    src += (4 * src_stride);

    ILVR_B_4VECS_SB(src6, src7, src8, src9, src7, src8, src9, src10,
                    src76_r, src87_r, src98_r, src109_r);

    ILVR_D_2VECS_SB(src8776, src87_r, src76_r, src10998, src109_r, src98_r);

    XORI_B_2VECS_SB(src8776, src10998, src8776, src10998, 128);

    out10 = FILT_8TAP_DPADD_S_H(src2110, src4332, src6554, src8776,
                                filt0, filt1, filt2, filt3);
    out32 = FILT_8TAP_DPADD_S_H(src4332, src6554, src8776, src10998,
                                filt0, filt1, filt2, filt3);

    out10 = SRARI_SATURATE_SIGNED_H(out10, FILTER_BITS, 7);
    out32 = SRARI_SATURATE_SIGNED_H(out32, FILTER_BITS, 7);

    PCKEV_2B_XORI128_STORE_4_BYTES_4(out10, out32, dst, dst_stride);
    dst += (4 * dst_stride);

    src2110 = src6554;
    src4332 = src8776;
    src6554 = src10998;

    src6 = src10;
  }
}

static void common_vt_8t_8w_msa(const uint8_t *src, int32_t src_stride,
                                uint8_t *dst, int32_t dst_stride,
                                int8_t *filter, int32_t height) {
  uint32_t loop_cnt;
  v16i8 src0, src1, src2, src3, src4, src5, src6, src7, src8, src9, src10;
  v16i8 src10_r, src32_r, src54_r, src76_r, src98_r;
  v16i8 src21_r, src43_r, src65_r, src87_r, src109_r;
  v16i8 filt0, filt1, filt2, filt3;
  v8i16 filt, out0_r, out1_r, out2_r, out3_r;

  src -= (3 * src_stride);

  filt = LOAD_SH(filter);
  filt0 = (v16i8)__msa_splati_h(filt, 0);
  filt1 = (v16i8)__msa_splati_h(filt, 1);
  filt2 = (v16i8)__msa_splati_h(filt, 2);
  filt3 = (v16i8)__msa_splati_h(filt, 3);

  LOAD_7VECS_SB(src, src_stride, src0, src1, src2, src3, src4, src5, src6);
  src += (7 * src_stride);

  XORI_B_7VECS_SB(src0, src1, src2, src3, src4, src5, src6,
                  src0, src1, src2, src3, src4, src5, src6, 128);

  ILVR_B_6VECS_SB(src0, src2, src4, src1, src3, src5,
                  src1, src3, src5, src2, src4, src6,
                  src10_r, src32_r, src54_r, src21_r, src43_r, src65_r);

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LOAD_4VECS_SB(src, src_stride, src7, src8, src9, src10);
    src += (4 * src_stride);

    XORI_B_4VECS_SB(src7, src8, src9, src10, src7, src8, src9, src10, 128);

    ILVR_B_4VECS_SB(src6, src7, src8, src9, src7, src8, src9, src10,
                    src76_r, src87_r, src98_r, src109_r);

    out0_r = FILT_8TAP_DPADD_S_H(src10_r, src32_r, src54_r, src76_r,
                                 filt0, filt1, filt2, filt3);
    out1_r = FILT_8TAP_DPADD_S_H(src21_r, src43_r, src65_r, src87_r,
                                 filt0, filt1, filt2, filt3);
    out2_r = FILT_8TAP_DPADD_S_H(src32_r, src54_r, src76_r, src98_r,
                                 filt0, filt1, filt2, filt3);
    out3_r = FILT_8TAP_DPADD_S_H(src43_r, src65_r, src87_r, src109_r,
                                 filt0, filt1, filt2, filt3);

    out0_r = SRARI_SATURATE_SIGNED_H(out0_r, FILTER_BITS, 7);
    out1_r = SRARI_SATURATE_SIGNED_H(out1_r, FILTER_BITS, 7);
    out2_r = SRARI_SATURATE_SIGNED_H(out2_r, FILTER_BITS, 7);
    out3_r = SRARI_SATURATE_SIGNED_H(out3_r, FILTER_BITS, 7);

    PCKEV_B_4_XORI128_STORE_8_BYTES_4(out0_r, out1_r, out2_r, out3_r,
                                      dst, dst_stride);
    dst += (4 * dst_stride);

    src10_r = src54_r;
    src32_r = src76_r;
    src54_r = src98_r;
    src21_r = src65_r;
    src43_r = src87_r;
    src65_r = src109_r;

    src6 = src10;
  }
}

static void common_vt_8t_16w_mult_msa(const uint8_t *src, int32_t src_stride,
                                      uint8_t *dst, int32_t dst_stride,
                                      int8_t *filter, int32_t height,
                                      int32_t width) {
  const uint8_t *src_tmp;
  uint8_t *dst_tmp;
  uint32_t loop_cnt, cnt;
  v16i8 src0, src1, src2, src3, src4, src5, src6, src7, src8, src9, src10;
  v16i8 filt0, filt1, filt2, filt3;
  v16i8 src10_r, src32_r, src54_r, src76_r, src98_r;
  v16i8 src21_r, src43_r, src65_r, src87_r, src109_r;
  v16i8 src10_l, src32_l, src54_l, src76_l, src98_l;
  v16i8 src21_l, src43_l, src65_l, src87_l, src109_l;
  v8i16 out0_r, out1_r, out2_r, out3_r, out0_l, out1_l, out2_l, out3_l;
  v8i16 filt;
  v16u8 tmp0, tmp1, tmp2, tmp3;

  src -= (3 * src_stride);

  filt = LOAD_SH(filter);
  filt0 = (v16i8)__msa_splati_h(filt, 0);
  filt1 = (v16i8)__msa_splati_h(filt, 1);
  filt2 = (v16i8)__msa_splati_h(filt, 2);
  filt3 = (v16i8)__msa_splati_h(filt, 3);

  for (cnt = (width >> 4); cnt--;) {
    src_tmp = src;
    dst_tmp = dst;

    LOAD_7VECS_SB(src_tmp, src_stride,
                  src0, src1, src2, src3, src4, src5, src6);
    src_tmp += (7 * src_stride);

    XORI_B_7VECS_SB(src0, src1, src2, src3, src4, src5, src6,
                    src0, src1, src2, src3, src4, src5, src6, 128);

    ILVR_B_6VECS_SB(src0, src2, src4, src1, src3, src5,
                    src1, src3, src5, src2, src4, src6,
                    src10_r, src32_r, src54_r, src21_r, src43_r, src65_r);

    ILVL_B_6VECS_SB(src0, src2, src4, src1, src3, src5,
                    src1, src3, src5, src2, src4, src6,
                    src10_l, src32_l, src54_l, src21_l, src43_l, src65_l);

    for (loop_cnt = (height >> 2); loop_cnt--;) {
      LOAD_4VECS_SB(src_tmp, src_stride, src7, src8, src9, src10);
      src_tmp += (4 * src_stride);

      XORI_B_4VECS_SB(src7, src8, src9, src10, src7, src8, src9, src10, 128);

      ILVR_B_4VECS_SB(src6, src7, src8, src9, src7, src8, src9, src10,
                      src76_r, src87_r, src98_r, src109_r);

      ILVL_B_4VECS_SB(src6, src7, src8, src9, src7, src8, src9, src10,
                      src76_l, src87_l, src98_l, src109_l);

      out0_r = FILT_8TAP_DPADD_S_H(src10_r, src32_r, src54_r, src76_r,
                                   filt0, filt1, filt2, filt3);
      out1_r = FILT_8TAP_DPADD_S_H(src21_r, src43_r, src65_r, src87_r,
                                   filt0, filt1, filt2, filt3);
      out2_r = FILT_8TAP_DPADD_S_H(src32_r, src54_r, src76_r, src98_r,
                                   filt0, filt1, filt2, filt3);
      out3_r = FILT_8TAP_DPADD_S_H(src43_r, src65_r, src87_r, src109_r,
                                   filt0, filt1, filt2, filt3);

      out0_l = FILT_8TAP_DPADD_S_H(src10_l, src32_l, src54_l, src76_l,
                                   filt0, filt1, filt2, filt3);
      out1_l = FILT_8TAP_DPADD_S_H(src21_l, src43_l, src65_l, src87_l,
                                   filt0, filt1, filt2, filt3);
      out2_l = FILT_8TAP_DPADD_S_H(src32_l, src54_l, src76_l, src98_l,
                                   filt0, filt1, filt2, filt3);
      out3_l = FILT_8TAP_DPADD_S_H(src43_l, src65_l, src87_l, src109_l,
                                   filt0, filt1, filt2, filt3);

      out0_r = SRARI_SATURATE_SIGNED_H(out0_r, FILTER_BITS, 7);
      out1_r = SRARI_SATURATE_SIGNED_H(out1_r, FILTER_BITS, 7);
      out2_r = SRARI_SATURATE_SIGNED_H(out2_r, FILTER_BITS, 7);
      out3_r = SRARI_SATURATE_SIGNED_H(out3_r, FILTER_BITS, 7);
      out0_l = SRARI_SATURATE_SIGNED_H(out0_l, FILTER_BITS, 7);
      out1_l = SRARI_SATURATE_SIGNED_H(out1_l, FILTER_BITS, 7);
      out2_l = SRARI_SATURATE_SIGNED_H(out2_l, FILTER_BITS, 7);
      out3_l = SRARI_SATURATE_SIGNED_H(out3_l, FILTER_BITS, 7);

      out0_r = (v8i16)__msa_pckev_b((v16i8)out0_l, (v16i8)out0_r);
      out1_r = (v8i16)__msa_pckev_b((v16i8)out1_l, (v16i8)out1_r);
      out2_r = (v8i16)__msa_pckev_b((v16i8)out2_l, (v16i8)out2_r);
      out3_r = (v8i16)__msa_pckev_b((v16i8)out3_l, (v16i8)out3_r);

      XORI_B_4VECS_UB(out0_r, out1_r, out2_r, out3_r,
                      tmp0, tmp1, tmp2, tmp3, 128);

      STORE_4VECS_UB(dst_tmp, dst_stride, tmp0, tmp1, tmp2, tmp3);
      dst_tmp += (4 * dst_stride);

      src10_r = src54_r;
      src32_r = src76_r;
      src54_r = src98_r;
      src21_r = src65_r;
      src43_r = src87_r;
      src65_r = src109_r;

      src10_l = src54_l;
      src32_l = src76_l;
      src54_l = src98_l;
      src21_l = src65_l;
      src43_l = src87_l;
      src65_l = src109_l;

      src6 = src10;
    }

    src += 16;
    dst += 16;
  }
}

static void common_vt_8t_16w_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter, int32_t height) {
  common_vt_8t_16w_mult_msa(src, src_stride, dst, dst_stride,
                            filter, height, 16);
}

static void common_vt_8t_32w_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter, int32_t height) {
  common_vt_8t_16w_mult_msa(src, src_stride, dst, dst_stride,
                            filter, height, 32);
}

static void common_vt_8t_64w_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter, int32_t height) {
  common_vt_8t_16w_mult_msa(src, src_stride, dst, dst_stride,
                            filter, height, 64);
}

static void common_vt_2t_4x4_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter) {
  uint32_t out0, out1, out2, out3;
  v16i8 src0, src1, src2, src3, src4;
  v16i8 src10_r, src32_r, src21_r, src43_r, src2110, src4332;
  v16i8 filt0;
  v8u16 filt;

  filt = LOAD_UH(filter);
  filt0 = (v16i8)__msa_splati_h((v8i16)filt, 0);

  LOAD_5VECS_SB(src, src_stride, src0, src1, src2, src3, src4);
  src += (5 * src_stride);

  ILVR_B_4VECS_SB(src0, src1, src2, src3, src1, src2, src3, src4,
                  src10_r, src21_r, src32_r, src43_r);

  ILVR_D_2VECS_SB(src2110, src21_r, src10_r, src4332, src43_r, src32_r);

  src2110 = (v16i8)__msa_dotp_u_h((v16u8)src2110, (v16u8)filt0);
  src4332 = (v16i8)__msa_dotp_u_h((v16u8)src4332, (v16u8)filt0);

  src2110 = (v16i8)SRARI_SATURATE_UNSIGNED_H(src2110, FILTER_BITS, 7);
  src4332 = (v16i8)SRARI_SATURATE_UNSIGNED_H(src4332, FILTER_BITS, 7);

  src2110 = (v16i8)__msa_pckev_b((v16i8)src4332, (v16i8)src2110);

  out0 = __msa_copy_u_w((v4i32)src2110, 0);
  out1 = __msa_copy_u_w((v4i32)src2110, 1);
  out2 = __msa_copy_u_w((v4i32)src2110, 2);
  out3 = __msa_copy_u_w((v4i32)src2110, 3);

  STORE_WORD(dst, out0);
  dst += dst_stride;
  STORE_WORD(dst, out1);
  dst += dst_stride;
  STORE_WORD(dst, out2);
  dst += dst_stride;
  STORE_WORD(dst, out3);
}

static void common_vt_2t_4x8_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter) {
  uint32_t out0, out1, out2, out3, out4, out5, out6, out7;
  v16i8 src0, src1, src2, src3, src4, src5, src6, src7, src8;
  v16i8 src10_r, src32_r, src54_r, src76_r, src21_r, src43_r;
  v16i8 src65_r, src87_r, src2110, src4332, src6554, src8776;
  v16i8 filt0;
  v8u16 filt;

  filt = LOAD_UH(filter);
  filt0 = (v16i8)__msa_splati_h((v8i16)filt, 0);

  LOAD_8VECS_SB(src, src_stride,
                src0, src1, src2, src3, src4, src5, src6, src7);
  src += (8 * src_stride);

  src8 = LOAD_SB(src);
  src += src_stride;

  ILVR_B_8VECS_SB(src0, src1, src2, src3, src4, src5, src6, src7,
                  src1, src2, src3, src4, src5, src6, src7, src8,
                  src10_r, src21_r, src32_r, src43_r,
                  src54_r, src65_r, src76_r, src87_r);

  ILVR_D_4VECS_SB(src2110, src21_r, src10_r, src4332, src43_r, src32_r,
                  src6554, src65_r, src54_r, src8776, src87_r, src76_r);

  src2110 = (v16i8)__msa_dotp_u_h((v16u8)src2110, (v16u8)filt0);
  src4332 = (v16i8)__msa_dotp_u_h((v16u8)src4332, (v16u8)filt0);
  src6554 = (v16i8)__msa_dotp_u_h((v16u8)src6554, (v16u8)filt0);
  src8776 = (v16i8)__msa_dotp_u_h((v16u8)src8776, (v16u8)filt0);

  src2110 = (v16i8)SRARI_SATURATE_UNSIGNED_H(src2110, FILTER_BITS, 7);
  src4332 = (v16i8)SRARI_SATURATE_UNSIGNED_H(src4332, FILTER_BITS, 7);
  src6554 = (v16i8)SRARI_SATURATE_UNSIGNED_H(src6554, FILTER_BITS, 7);
  src8776 = (v16i8)SRARI_SATURATE_UNSIGNED_H(src8776, FILTER_BITS, 7);

  src2110 = (v16i8)__msa_pckev_b((v16i8)src4332, (v16i8)src2110);
  src4332 = (v16i8)__msa_pckev_b((v16i8)src8776, (v16i8)src6554);

  out0 = __msa_copy_u_w((v4i32)src2110, 0);
  out1 = __msa_copy_u_w((v4i32)src2110, 1);
  out2 = __msa_copy_u_w((v4i32)src2110, 2);
  out3 = __msa_copy_u_w((v4i32)src2110, 3);
  out4 = __msa_copy_u_w((v4i32)src4332, 0);
  out5 = __msa_copy_u_w((v4i32)src4332, 1);
  out6 = __msa_copy_u_w((v4i32)src4332, 2);
  out7 = __msa_copy_u_w((v4i32)src4332, 3);

  STORE_WORD(dst, out0);
  dst += dst_stride;
  STORE_WORD(dst, out1);
  dst += dst_stride;
  STORE_WORD(dst, out2);
  dst += dst_stride;
  STORE_WORD(dst, out3);
  dst += dst_stride;
  STORE_WORD(dst, out4);
  dst += dst_stride;
  STORE_WORD(dst, out5);
  dst += dst_stride;
  STORE_WORD(dst, out6);
  dst += dst_stride;
  STORE_WORD(dst, out7);
}

static void common_vt_2t_4w_msa(const uint8_t *src, int32_t src_stride,
                                uint8_t *dst, int32_t dst_stride,
                                int8_t *filter, int32_t height) {
  if (4 == height) {
    common_vt_2t_4x4_msa(src, src_stride, dst, dst_stride, filter);
  } else if (8 == height) {
    common_vt_2t_4x8_msa(src, src_stride, dst, dst_stride, filter);
  }
}

static void common_vt_2t_8x4_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter) {
  v16u8 src0, src1, src2, src3, src4;
  v16u8 vec0, vec1, vec2, vec3, filt0;
  v8u16 tmp0, tmp1, tmp2, tmp3;
  v8u16 filt;

  /* rearranging filter_y */
  filt = LOAD_UH(filter);
  filt0 = (v16u8)__msa_splati_h((v8i16)filt, 0);

  LOAD_5VECS_UB(src, src_stride, src0, src1, src2, src3, src4);

  ILVR_B_2VECS_UB(src0, src1, src1, src2, vec0, vec1);
  ILVR_B_2VECS_UB(src2, src3, src3, src4, vec2, vec3);

  /* filter calc */
  tmp0 = __msa_dotp_u_h(vec0, filt0);
  tmp1 = __msa_dotp_u_h(vec1, filt0);
  tmp2 = __msa_dotp_u_h(vec2, filt0);
  tmp3 = __msa_dotp_u_h(vec3, filt0);

  tmp0 = SRARI_SATURATE_UNSIGNED_H(tmp0, FILTER_BITS, 7);
  tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);
  tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);
  tmp3 = SRARI_SATURATE_UNSIGNED_H(tmp3, FILTER_BITS, 7);

  PCKEV_B_STORE_8_BYTES_4(tmp0, tmp1, tmp2, tmp3, dst, dst_stride);
}

static void common_vt_2t_8x8mult_msa(const uint8_t *src, int32_t src_stride,
                                     uint8_t *dst, int32_t dst_stride,
                                     int8_t *filter, int32_t height) {
  uint32_t loop_cnt;
  v16u8 src0, src1, src2, src3, src4, src5, src6, src7, src8;
  v16u8 vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7, filt0;
  v8u16 tmp0, tmp1, tmp2, tmp3;
  v8u16 filt;

  /* rearranging filter_y */
  filt = LOAD_UH(filter);
  filt0 = (v16u8)__msa_splati_h((v8i16)filt, 0);

  src0 = LOAD_UB(src);
  src += src_stride;

  for (loop_cnt = (height >> 3); loop_cnt--;) {
    LOAD_8VECS_UB(src, src_stride,
                  src1, src2, src3, src4, src5, src6, src7, src8);
    src += (8 * src_stride);

    ILVR_B_4VECS_UB(src0, src1, src2, src3, src1, src2, src3, src4,
                    vec0, vec1, vec2, vec3);

    ILVR_B_4VECS_UB(src4, src5, src6, src7, src5, src6, src7, src8,
                    vec4, vec5, vec6, vec7);

    tmp0 = __msa_dotp_u_h(vec0, filt0);
    tmp1 = __msa_dotp_u_h(vec1, filt0);
    tmp2 = __msa_dotp_u_h(vec2, filt0);
    tmp3 = __msa_dotp_u_h(vec3, filt0);

    tmp0 = SRARI_SATURATE_UNSIGNED_H(tmp0, FILTER_BITS, 7);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);
    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);
    tmp3 = SRARI_SATURATE_UNSIGNED_H(tmp3, FILTER_BITS, 7);

    PCKEV_B_STORE_8_BYTES_4(tmp0, tmp1, tmp2, tmp3, dst, dst_stride);
    dst += (4 * dst_stride);

    tmp0 = __msa_dotp_u_h(vec4, filt0);
    tmp1 = __msa_dotp_u_h(vec5, filt0);
    tmp2 = __msa_dotp_u_h(vec6, filt0);
    tmp3 = __msa_dotp_u_h(vec7, filt0);

    tmp0 = SRARI_SATURATE_UNSIGNED_H(tmp0, FILTER_BITS, 7);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);
    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);
    tmp3 = SRARI_SATURATE_UNSIGNED_H(tmp3, FILTER_BITS, 7);

    PCKEV_B_STORE_8_BYTES_4(tmp0, tmp1, tmp2, tmp3, dst, dst_stride);
    dst += (4 * dst_stride);

    src0 = src8;
  }
}

static void common_vt_2t_8w_msa(const uint8_t *src, int32_t src_stride,
                                uint8_t *dst, int32_t dst_stride,
                                int8_t *filter, int32_t height) {
  if (4 == height) {
    common_vt_2t_8x4_msa(src, src_stride, dst, dst_stride, filter);
  } else {
    common_vt_2t_8x8mult_msa(src, src_stride, dst, dst_stride, filter, height);
  }
}

static void common_vt_2t_16w_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter, int32_t height) {
  uint32_t loop_cnt;
  v16u8 src0, src1, src2, src3, src4;
  v16u8 vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7, filt0;
  v8u16 tmp0, tmp1, tmp2, tmp3;
  v8u16 filt;

  /* rearranging filter_y */
  filt = LOAD_UH(filter);
  filt0 = (v16u8)__msa_splati_h((v8i16)filt, 0);

  src0 = LOAD_UB(src);
  src += src_stride;

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LOAD_4VECS_UB(src, src_stride, src1, src2, src3, src4);
    src += (4 * src_stride);

    ILV_B_LRLR_UB(src0, src1, src1, src2, vec1, vec0, vec3, vec2);

    tmp0 = __msa_dotp_u_h(vec0, filt0);
    tmp1 = __msa_dotp_u_h(vec1, filt0);

    tmp0 = SRARI_SATURATE_UNSIGNED_H(tmp0, FILTER_BITS, 7);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp1, tmp0, dst);
    dst += dst_stride;

    ILV_B_LRLR_UB(src2, src3, src3, src4, vec5, vec4, vec7, vec6);

    tmp2 = __msa_dotp_u_h(vec2, filt0);
    tmp3 = __msa_dotp_u_h(vec3, filt0);

    tmp3 = SRARI_SATURATE_UNSIGNED_H(tmp3, FILTER_BITS, 7);
    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp3, tmp2, dst);
    dst += dst_stride;

    tmp0 = __msa_dotp_u_h(vec4, filt0);
    tmp1 = __msa_dotp_u_h(vec5, filt0);

    tmp0 = SRARI_SATURATE_UNSIGNED_H(tmp0, FILTER_BITS, 7);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp1, tmp0, dst);
    dst += dst_stride;

    tmp2 = __msa_dotp_u_h(vec6, filt0);
    tmp3 = __msa_dotp_u_h(vec7, filt0);

    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);
    tmp3 = SRARI_SATURATE_UNSIGNED_H(tmp3, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp3, tmp2, dst);
    dst += dst_stride;

    src0 = src4;
  }
}

static void common_vt_2t_32w_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter, int32_t height) {
  uint32_t loop_cnt;
  v16u8 src0, src1, src2, src3, src4, src5, src6, src7, src8, src9;
  v16u8 vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7, filt0;
  v8u16 tmp0, tmp1, tmp2, tmp3;
  v8u16 filt;

  /* rearranging filter_y */
  filt = LOAD_UH(filter);
  filt0 = (v16u8)__msa_splati_h((v8i16)filt, 0);

  src0 = LOAD_UB(src);
  src5 = LOAD_UB(src + 16);
  src += src_stride;

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LOAD_4VECS_UB(src, src_stride, src1, src2, src3, src4);

    ILV_B_LRLR_UB(src0, src1, src1, src2, vec1, vec0, vec3, vec2);

    LOAD_4VECS_UB(src + 16, src_stride, src6, src7, src8, src9);
    src += (4 * src_stride);

    tmp0 = __msa_dotp_u_h(vec0, filt0);
    tmp1 = __msa_dotp_u_h(vec1, filt0);

    tmp0 = SRARI_SATURATE_UNSIGNED_H(tmp0, FILTER_BITS, 7);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp1, tmp0, dst);

    tmp2 = __msa_dotp_u_h(vec2, filt0);
    tmp3 = __msa_dotp_u_h(vec3, filt0);

    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);
    tmp3 = SRARI_SATURATE_UNSIGNED_H(tmp3, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp3, tmp2, dst + dst_stride);

    ILV_B_LRLR_UB(src2, src3, src3, src4, vec5, vec4, vec7, vec6);

    tmp0 = __msa_dotp_u_h(vec4, filt0);
    tmp1 = __msa_dotp_u_h(vec5, filt0);

    tmp0 = SRARI_SATURATE_UNSIGNED_H(tmp0, FILTER_BITS, 7);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp1, tmp0, dst + 2 * dst_stride);

    tmp2 = __msa_dotp_u_h(vec6, filt0);
    tmp3 = __msa_dotp_u_h(vec7, filt0);

    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);
    tmp3 = SRARI_SATURATE_UNSIGNED_H(tmp3, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp3, tmp2, dst + 3 * dst_stride);

    ILV_B_LRLR_UB(src5, src6, src6, src7, vec1, vec0, vec3, vec2);

    tmp0 = __msa_dotp_u_h(vec0, filt0);
    tmp1 = __msa_dotp_u_h(vec1, filt0);

    tmp0 = SRARI_SATURATE_UNSIGNED_H(tmp0, FILTER_BITS, 7);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp1, tmp0, dst + 16);

    tmp2 = __msa_dotp_u_h(vec2, filt0);
    tmp3 = __msa_dotp_u_h(vec3, filt0);

    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);
    tmp3 = SRARI_SATURATE_UNSIGNED_H(tmp3, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp3, tmp2, dst + 16 + dst_stride);

    ILV_B_LRLR_UB(src7, src8, src8, src9, vec5, vec4, vec7, vec6);

    tmp0 = __msa_dotp_u_h(vec4, filt0);
    tmp1 = __msa_dotp_u_h(vec5, filt0);

    tmp0 = SRARI_SATURATE_UNSIGNED_H(tmp0, FILTER_BITS, 7);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp1, tmp0, dst + 16 + 2 * dst_stride);

    tmp2 = __msa_dotp_u_h(vec6, filt0);
    tmp3 = __msa_dotp_u_h(vec7, filt0);

    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);
    tmp3 = SRARI_SATURATE_UNSIGNED_H(tmp3, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp3, tmp2, dst + 16 + 3 * dst_stride);
    dst += (4 * dst_stride);

    src0 = src4;
    src5 = src9;
  }
}

static void common_vt_2t_64w_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter, int32_t height) {
  uint32_t loop_cnt;
  v16u8 src0, src1, src2, src3, src4, src5, src6, src7;
  v16u8 src8, src9, src10, src11;
  v16u8 vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7, filt0;
  v8u16 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  v8u16 filt;

  /* rearranging filter_y */
  filt = LOAD_UH(filter);
  filt0 = (v16u8)__msa_splati_h((v8i16)filt, 0);

  LOAD_4VECS_UB(src, 16, src0, src3, src6, src9);
  src += src_stride;

  for (loop_cnt = (height >> 1); loop_cnt--;) {
    LOAD_2VECS_UB(src, src_stride, src1, src2);
    LOAD_2VECS_UB(src + 16, src_stride, src4, src5);
    LOAD_2VECS_UB(src + 32, src_stride, src7, src8);
    LOAD_2VECS_UB(src + 48, src_stride, src10, src11);
    src += (2 * src_stride);

    ILV_B_LRLR_UB(src0, src1, src1, src2, vec1, vec0, vec3, vec2);

    tmp0 = __msa_dotp_u_h(vec0, filt0);
    tmp1 = __msa_dotp_u_h(vec1, filt0);

    tmp0 = SRARI_SATURATE_UNSIGNED_H(tmp0, FILTER_BITS, 7);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp1, tmp0, dst);

    tmp2 = __msa_dotp_u_h(vec2, filt0);
    tmp3 = __msa_dotp_u_h(vec3, filt0);

    tmp3 = SRARI_SATURATE_UNSIGNED_H(tmp3, FILTER_BITS, 7);
    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp3, tmp2, dst + dst_stride);

    ILV_B_LRLR_UB(src3, src4, src4, src5, vec5, vec4, vec7, vec6);

    tmp4 = __msa_dotp_u_h(vec4, filt0);
    tmp5 = __msa_dotp_u_h(vec5, filt0);

    tmp4 = SRARI_SATURATE_UNSIGNED_H(tmp4, FILTER_BITS, 7);
    tmp5 = SRARI_SATURATE_UNSIGNED_H(tmp5, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp5, tmp4, dst + 16);

    tmp6 = __msa_dotp_u_h(vec6, filt0);
    tmp7 = __msa_dotp_u_h(vec7, filt0);

    tmp6 = SRARI_SATURATE_UNSIGNED_H(tmp6, FILTER_BITS, 7);
    tmp7 = SRARI_SATURATE_UNSIGNED_H(tmp7, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp7, tmp6, dst + 16 + dst_stride);

    ILV_B_LRLR_UB(src6, src7, src7, src8, vec1, vec0, vec3, vec2);

    tmp0 = __msa_dotp_u_h(vec0, filt0);
    tmp1 = __msa_dotp_u_h(vec1, filt0);

    tmp0 = SRARI_SATURATE_UNSIGNED_H(tmp0, FILTER_BITS, 7);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp1, tmp0, dst + 32);

    tmp2 = __msa_dotp_u_h(vec2, filt0);
    tmp3 = __msa_dotp_u_h(vec3, filt0);

    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);
    tmp3 = SRARI_SATURATE_UNSIGNED_H(tmp3, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp3, tmp2, dst + 32 + dst_stride);

    ILV_B_LRLR_UB(src9, src10, src10, src11, vec5, vec4, vec7, vec6);

    tmp4 = __msa_dotp_u_h(vec4, filt0);
    tmp5 = __msa_dotp_u_h(vec5, filt0);

    tmp4 = SRARI_SATURATE_UNSIGNED_H(tmp4, FILTER_BITS, 7);
    tmp5 = SRARI_SATURATE_UNSIGNED_H(tmp5, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp5, tmp4, dst + 48);

    tmp6 = __msa_dotp_u_h(vec6, filt0);
    tmp7 = __msa_dotp_u_h(vec7, filt0);

    tmp6 = SRARI_SATURATE_UNSIGNED_H(tmp6, FILTER_BITS, 7);
    tmp7 = SRARI_SATURATE_UNSIGNED_H(tmp7, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp7, tmp6, dst + 48 + dst_stride);
    dst += (2 * dst_stride);

    src0 = src2;
    src3 = src5;
    src6 = src8;
    src9 = src11;
  }
}

void vp9_convolve8_vert_msa(const uint8_t *src, ptrdiff_t src_stride,
                            uint8_t *dst, ptrdiff_t dst_stride,
                            const int16_t *filter_x, int x_step_q4,
                            const int16_t *filter_y, int y_step_q4,
                            int w, int h) {
  int8_t cnt, filt_ver[8];

  if (16 != y_step_q4) {
    vp9_convolve8_vert_c(src, src_stride, dst, dst_stride,
                         filter_x, x_step_q4, filter_y, y_step_q4,
                         w, h);
    return;
  }

  if (((const int32_t *)filter_y)[1] == 0x800000) {
    vp9_convolve_copy(src, src_stride, dst, dst_stride,
                      filter_x, x_step_q4, filter_y, y_step_q4,
                      w, h);
    return;
  }

  for (cnt = 8; cnt--;) {
    filt_ver[cnt] = filter_y[cnt];
  }

  if (((const int32_t *)filter_y)[0] == 0) {
    switch (w) {
      case 4:
        common_vt_2t_4w_msa(src, (int32_t)src_stride,
                            dst, (int32_t)dst_stride,
                            &filt_ver[3], h);
        break;
      case 8:
        common_vt_2t_8w_msa(src, (int32_t)src_stride,
                            dst, (int32_t)dst_stride,
                            &filt_ver[3], h);
        break;
      case 16:
        common_vt_2t_16w_msa(src, (int32_t)src_stride,
                             dst, (int32_t)dst_stride,
                             &filt_ver[3], h);
        break;
      case 32:
        common_vt_2t_32w_msa(src, (int32_t)src_stride,
                             dst, (int32_t)dst_stride,
                             &filt_ver[3], h);
        break;
      case 64:
        common_vt_2t_64w_msa(src, (int32_t)src_stride,
                             dst, (int32_t)dst_stride,
                             &filt_ver[3], h);
        break;
      default:
        vp9_convolve8_vert_c(src, src_stride, dst, dst_stride,
                             filter_x, x_step_q4, filter_y, y_step_q4,
                             w, h);
        break;
    }
  } else {
    switch (w) {
      case 4:
        common_vt_8t_4w_msa(src, (int32_t)src_stride,
                            dst, (int32_t)dst_stride,
                            filt_ver, h);
        break;
      case 8:
        common_vt_8t_8w_msa(src, (int32_t)src_stride,
                            dst, (int32_t)dst_stride,
                            filt_ver, h);
        break;
      case 16:
        common_vt_8t_16w_msa(src, (int32_t)src_stride,
                             dst, (int32_t)dst_stride,
                             filt_ver, h);
        break;
      case 32:
        common_vt_8t_32w_msa(src, (int32_t)src_stride,
                             dst, (int32_t)dst_stride,
                             filt_ver, h);
        break;
      case 64:
        common_vt_8t_64w_msa(src, (int32_t)src_stride,
                             dst, (int32_t)dst_stride,
                             filt_ver, h);
        break;
      default:
        vp9_convolve8_vert_c(src, src_stride, dst, dst_stride,
                             filter_x, x_step_q4, filter_y, y_step_q4,
                             w, h);
        break;
    }
  }
}
