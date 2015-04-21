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

static void common_hz_8t_4x4_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter) {
  v16i8 filt0, filt1, filt2, filt3;
  v16i8 src0, src1, src2, src3;
  v16u8 mask0, mask1, mask2, mask3;
  v8i16 filt, out0, out1;

  mask0 = LOAD_UB(&mc_filt_mask_arr[16]);

  src -= 3;

  /* rearranging filter */
  filt = LOAD_SH(filter);
  filt0 = (v16i8)__msa_splati_h(filt, 0);
  filt1 = (v16i8)__msa_splati_h(filt, 1);
  filt2 = (v16i8)__msa_splati_h(filt, 2);
  filt3 = (v16i8)__msa_splati_h(filt, 3);

  mask1 = mask0 + 2;
  mask2 = mask0 + 4;
  mask3 = mask0 + 6;

  LOAD_4VECS_SB(src, src_stride, src0, src1, src2, src3);

  XORI_B_4VECS_SB(src0, src1, src2, src3, src0, src1, src2, src3, 128);

  HORIZ_8TAP_4WID_4VECS_FILT(src0, src1, src2, src3, mask0, mask1, mask2, mask3,
                             filt0, filt1, filt2, filt3, out0, out1);

  out0 = SRARI_SATURATE_SIGNED_H(out0, FILTER_BITS, 7);
  out1 = SRARI_SATURATE_SIGNED_H(out1, FILTER_BITS, 7);

  PCKEV_2B_XORI128_STORE_4_BYTES_4(out0, out1, dst, dst_stride);
}

static void common_hz_8t_4x8_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter) {
  v16i8 filt0, filt1, filt2, filt3;
  v16i8 src0, src1, src2, src3;
  v16u8 mask0, mask1, mask2, mask3;
  v8i16 filt, out0, out1, out2, out3;

  mask0 = LOAD_UB(&mc_filt_mask_arr[16]);

  src -= 3;

  /* rearranging filter */
  filt = LOAD_SH(filter);
  filt0 = (v16i8)__msa_splati_h(filt, 0);
  filt1 = (v16i8)__msa_splati_h(filt, 1);
  filt2 = (v16i8)__msa_splati_h(filt, 2);
  filt3 = (v16i8)__msa_splati_h(filt, 3);

  mask1 = mask0 + 2;
  mask2 = mask0 + 4;
  mask3 = mask0 + 6;

  LOAD_4VECS_SB(src, src_stride, src0, src1, src2, src3);
  src += (4 * src_stride);

  XORI_B_4VECS_SB(src0, src1, src2, src3, src0, src1, src2, src3, 128);

  HORIZ_8TAP_4WID_4VECS_FILT(src0, src1, src2, src3, mask0, mask1, mask2, mask3,
                             filt0, filt1, filt2, filt3, out0, out1);

  LOAD_4VECS_SB(src, src_stride, src0, src1, src2, src3);

  XORI_B_4VECS_SB(src0, src1, src2, src3, src0, src1, src2, src3, 128);

  HORIZ_8TAP_4WID_4VECS_FILT(src0, src1, src2, src3, mask0, mask1, mask2, mask3,
                             filt0, filt1, filt2, filt3, out2, out3);

  out0 = SRARI_SATURATE_SIGNED_H(out0, FILTER_BITS, 7);
  out1 = SRARI_SATURATE_SIGNED_H(out1, FILTER_BITS, 7);
  out2 = SRARI_SATURATE_SIGNED_H(out2, FILTER_BITS, 7);
  out3 = SRARI_SATURATE_SIGNED_H(out3, FILTER_BITS, 7);

  PCKEV_2B_XORI128_STORE_4_BYTES_4(out0, out1, dst, dst_stride);
  dst += (4 * dst_stride);
  PCKEV_2B_XORI128_STORE_4_BYTES_4(out2, out3, dst, dst_stride);
}

static void common_hz_8t_4w_msa(const uint8_t *src, int32_t src_stride,
                                uint8_t *dst, int32_t dst_stride,
                                int8_t *filter, int32_t height) {
  if (4 == height) {
    common_hz_8t_4x4_msa(src, src_stride, dst, dst_stride, filter);
  } else if (8 == height) {
    common_hz_8t_4x8_msa(src, src_stride, dst, dst_stride, filter);
  }
}

static void common_hz_8t_8x4_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter) {
  v16i8 filt0, filt1, filt2, filt3;
  v16i8 src0, src1, src2, src3;
  v16u8 mask0, mask1, mask2, mask3;
  v8i16 filt, out0, out1, out2, out3;

  mask0 = LOAD_UB(&mc_filt_mask_arr[0]);

  src -= 3;

  /* rearranging filter */
  filt = LOAD_SH(filter);
  filt0 = (v16i8)__msa_splati_h(filt, 0);
  filt1 = (v16i8)__msa_splati_h(filt, 1);
  filt2 = (v16i8)__msa_splati_h(filt, 2);
  filt3 = (v16i8)__msa_splati_h(filt, 3);

  mask1 = mask0 + 2;
  mask2 = mask0 + 4;
  mask3 = mask0 + 6;

  LOAD_4VECS_SB(src, src_stride, src0, src1, src2, src3);

  XORI_B_4VECS_SB(src0, src1, src2, src3, src0, src1, src2, src3, 128);

  HORIZ_8TAP_8WID_4VECS_FILT(src0, src1, src2, src3, mask0, mask1, mask2, mask3,
                             filt0, filt1, filt2, filt3, out0, out1, out2,
                             out3);

  out0 = SRARI_SATURATE_SIGNED_H(out0, FILTER_BITS, 7);
  out1 = SRARI_SATURATE_SIGNED_H(out1, FILTER_BITS, 7);
  out2 = SRARI_SATURATE_SIGNED_H(out2, FILTER_BITS, 7);
  out3 = SRARI_SATURATE_SIGNED_H(out3, FILTER_BITS, 7);

  PCKEV_B_4_XORI128_STORE_8_BYTES_4(out0, out1, out2, out3, dst, dst_stride);
}

static void common_hz_8t_8x8mult_msa(const uint8_t *src, int32_t src_stride,
                                     uint8_t *dst, int32_t dst_stride,
                                     int8_t *filter, int32_t height) {
  uint32_t loop_cnt;
  v16i8 filt0, filt1, filt2, filt3;
  v16i8 src0, src1, src2, src3;
  v16u8 mask0, mask1, mask2, mask3;
  v8i16 filt, out0, out1, out2, out3;

  mask0 = LOAD_UB(&mc_filt_mask_arr[0]);

  src -= 3;

  /* rearranging filter */
  filt = LOAD_SH(filter);
  filt0 = (v16i8)__msa_splati_h(filt, 0);
  filt1 = (v16i8)__msa_splati_h(filt, 1);
  filt2 = (v16i8)__msa_splati_h(filt, 2);
  filt3 = (v16i8)__msa_splati_h(filt, 3);

  mask1 = mask0 + 2;
  mask2 = mask0 + 4;
  mask3 = mask0 + 6;

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LOAD_4VECS_SB(src, src_stride, src0, src1, src2, src3);
    src += (4 * src_stride);

    XORI_B_4VECS_SB(src0, src1, src2, src3, src0, src1, src2, src3, 128);

    HORIZ_8TAP_8WID_4VECS_FILT(src0, src1, src2, src3, mask0, mask1, mask2,
                               mask3, filt0, filt1, filt2, filt3, out0, out1,
                               out2, out3);

    out0 = SRARI_SATURATE_SIGNED_H(out0, FILTER_BITS, 7);
    out1 = SRARI_SATURATE_SIGNED_H(out1, FILTER_BITS, 7);
    out2 = SRARI_SATURATE_SIGNED_H(out2, FILTER_BITS, 7);
    out3 = SRARI_SATURATE_SIGNED_H(out3, FILTER_BITS, 7);

    PCKEV_B_4_XORI128_STORE_8_BYTES_4(out0, out1, out2, out3, dst, dst_stride);
    dst += (4 * dst_stride);
  }
}

static void common_hz_8t_8w_msa(const uint8_t *src, int32_t src_stride,
                                uint8_t *dst, int32_t dst_stride,
                                int8_t *filter, int32_t height) {
  if (4 == height) {
    common_hz_8t_8x4_msa(src, src_stride, dst, dst_stride, filter);
  } else {
    common_hz_8t_8x8mult_msa(src, src_stride, dst, dst_stride, filter, height);
  }
}

static void common_hz_8t_16w_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter, int32_t height) {
  uint32_t loop_cnt;
  v16i8 src0, src1, src2, src3;
  v16i8 filt0, filt1, filt2, filt3;
  v16u8 mask0, mask1, mask2, mask3;
  v8i16 filt, out0, out1, out2, out3;

  mask0 = LOAD_UB(&mc_filt_mask_arr[0]);

  src -= 3;

  /* rearranging filter */
  filt = LOAD_SH(filter);
  filt0 = (v16i8)__msa_splati_h(filt, 0);
  filt1 = (v16i8)__msa_splati_h(filt, 1);
  filt2 = (v16i8)__msa_splati_h(filt, 2);
  filt3 = (v16i8)__msa_splati_h(filt, 3);

  mask1 = mask0 + 2;
  mask2 = mask0 + 4;
  mask3 = mask0 + 6;

  for (loop_cnt = (height >> 1); loop_cnt--;) {
    src0 = LOAD_SB(src);
    src1 = LOAD_SB(src + 8);
    src += src_stride;
    src2 = LOAD_SB(src);
    src3 = LOAD_SB(src + 8);
    src += src_stride;

    XORI_B_4VECS_SB(src0, src1, src2, src3, src0, src1, src2, src3, 128);

    HORIZ_8TAP_8WID_4VECS_FILT(src0, src1, src2, src3, mask0, mask1, mask2,
                               mask3, filt0, filt1, filt2, filt3, out0, out1,
                               out2, out3);

    out0 = SRARI_SATURATE_SIGNED_H(out0, FILTER_BITS, 7);
    out1 = SRARI_SATURATE_SIGNED_H(out1, FILTER_BITS, 7);
    out2 = SRARI_SATURATE_SIGNED_H(out2, FILTER_BITS, 7);
    out3 = SRARI_SATURATE_SIGNED_H(out3, FILTER_BITS, 7);

    PCKEV_B_XORI128_STORE_VEC(out1, out0, dst);
    dst += dst_stride;
    PCKEV_B_XORI128_STORE_VEC(out3, out2, dst);
    dst += dst_stride;
  }
}

static void common_hz_8t_32w_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter, int32_t height) {
  uint32_t loop_cnt;
  v16i8 src0, src1, src2, src3;
  v16i8 filt0, filt1, filt2, filt3;
  v16u8 mask0, mask1, mask2, mask3;
  v8i16 filt, out0, out1, out2, out3;

  mask0 = LOAD_UB(&mc_filt_mask_arr[0]);

  src -= 3;

  /* rearranging filter */
  filt = LOAD_SH(filter);
  filt0 = (v16i8)__msa_splati_h(filt, 0);
  filt1 = (v16i8)__msa_splati_h(filt, 1);
  filt2 = (v16i8)__msa_splati_h(filt, 2);
  filt3 = (v16i8)__msa_splati_h(filt, 3);

  mask1 = mask0 + 2;
  mask2 = mask0 + 4;
  mask3 = mask0 + 6;

  for (loop_cnt = (height >> 1); loop_cnt--;) {
    src0 = LOAD_SB(src);
    src2 = LOAD_SB(src + 16);
    src3 = LOAD_SB(src + 24);
    src1 = __msa_sld_b((v16i8)src2, (v16i8)src0, 8);
    src += src_stride;

    XORI_B_4VECS_SB(src0, src1, src2, src3, src0, src1, src2, src3, 128);

    HORIZ_8TAP_8WID_4VECS_FILT(src0, src1, src2, src3, mask0, mask1, mask2,
                               mask3, filt0, filt1, filt2, filt3, out0, out1,
                               out2, out3);

    out0 = SRARI_SATURATE_SIGNED_H(out0, FILTER_BITS, 7);
    out1 = SRARI_SATURATE_SIGNED_H(out1, FILTER_BITS, 7);
    out2 = SRARI_SATURATE_SIGNED_H(out2, FILTER_BITS, 7);
    out3 = SRARI_SATURATE_SIGNED_H(out3, FILTER_BITS, 7);

    src0 = LOAD_SB(src);
    src2 = LOAD_SB(src + 16);
    src3 = LOAD_SB(src + 24);
    src1 = __msa_sld_b((v16i8)src2, (v16i8)src0, 8);

    PCKEV_B_XORI128_STORE_VEC(out1, out0, dst);
    PCKEV_B_XORI128_STORE_VEC(out3, out2, (dst + 16));
    dst += dst_stride;

    XORI_B_4VECS_SB(src0, src1, src2, src3, src0, src1, src2, src3, 128);

    HORIZ_8TAP_8WID_4VECS_FILT(src0, src1, src2, src3, mask0, mask1, mask2,
                               mask3, filt0, filt1, filt2, filt3, out0, out1,
                               out2, out3);

    out0 = SRARI_SATURATE_SIGNED_H(out0, FILTER_BITS, 7);
    out1 = SRARI_SATURATE_SIGNED_H(out1, FILTER_BITS, 7);
    out2 = SRARI_SATURATE_SIGNED_H(out2, FILTER_BITS, 7);
    out3 = SRARI_SATURATE_SIGNED_H(out3, FILTER_BITS, 7);

    PCKEV_B_XORI128_STORE_VEC(out1, out0, dst);
    PCKEV_B_XORI128_STORE_VEC(out3, out2, (dst + 16));

    src += src_stride;
    dst += dst_stride;
  }
}

static void common_hz_8t_64w_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter, int32_t height) {
  uint32_t loop_cnt, cnt;
  v16i8 src0, src1, src2, src3;
  v16i8 filt0, filt1, filt2, filt3;
  v16u8 mask0, mask1, mask2, mask3;
  v8i16 filt, out0, out1, out2, out3;

  mask0 = LOAD_UB(&mc_filt_mask_arr[0]);

  src -= 3;

  /* rearranging filter */
  filt = LOAD_SH(filter);
  filt0 = (v16i8)__msa_splati_h(filt, 0);
  filt1 = (v16i8)__msa_splati_h(filt, 1);
  filt2 = (v16i8)__msa_splati_h(filt, 2);
  filt3 = (v16i8)__msa_splati_h(filt, 3);

  mask1 = mask0 + 2;
  mask2 = mask0 + 4;
  mask3 = mask0 + 6;

  for (loop_cnt = height; loop_cnt--;) {
    for (cnt = 0; cnt < 2; ++cnt) {
      src0 = LOAD_SB(&src[cnt << 5]);
      src2 = LOAD_SB(&src[16 + (cnt << 5)]);
      src3 = LOAD_SB(&src[24 + (cnt << 5)]);
      src1 = __msa_sld_b((v16i8)src2, (v16i8)src0, 8);

      XORI_B_4VECS_SB(src0, src1, src2, src3, src0, src1, src2, src3, 128);

      HORIZ_8TAP_8WID_4VECS_FILT(src0, src1, src2, src3, mask0, mask1, mask2,
                                 mask3, filt0, filt1, filt2, filt3, out0, out1,
                                 out2, out3);

      out0 = SRARI_SATURATE_SIGNED_H(out0, FILTER_BITS, 7);
      out1 = SRARI_SATURATE_SIGNED_H(out1, FILTER_BITS, 7);
      out2 = SRARI_SATURATE_SIGNED_H(out2, FILTER_BITS, 7);
      out3 = SRARI_SATURATE_SIGNED_H(out3, FILTER_BITS, 7);

      PCKEV_B_XORI128_STORE_VEC(out1, out0, &dst[cnt << 5]);
      PCKEV_B_XORI128_STORE_VEC(out3, out2, &dst[16 + (cnt << 5)]);
    }

    src += src_stride;
    dst += dst_stride;
  }
}

static void common_hz_2t_4x4_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter) {
  uint32_t out0, out1, out2, out3;
  v16i8 src0, src1, src2, src3, mask;
  v16u8 vec0, vec1, filt0;
  v16i8 res0, res1;
  v8u16 vec2, vec3, filt, const255;

  mask = LOAD_SB(&mc_filt_mask_arr[16]);

  /* rearranging filter */
  filt = LOAD_UH(filter);
  filt0 = (v16u8)__msa_splati_h((v8i16)filt, 0);

  const255 = (v8u16)__msa_ldi_h(255);

  LOAD_4VECS_SB(src, src_stride, src0, src1, src2, src3);

  vec0 = (v16u8)__msa_vshf_b(mask, src1, src0);
  vec1 = (v16u8)__msa_vshf_b(mask, src3, src2);

  vec2 = __msa_dotp_u_h(vec0, filt0);
  vec3 = __msa_dotp_u_h(vec1, filt0);

  vec2 = (v8u16)__msa_srari_h((v8i16)vec2, FILTER_BITS);
  vec3 = (v8u16)__msa_srari_h((v8i16)vec3, FILTER_BITS);

  vec2 = __msa_min_u_h(vec2, const255);
  vec3 = __msa_min_u_h(vec3, const255);

  res0 = __msa_pckev_b((v16i8)vec2, (v16i8)vec2);
  res1 = __msa_pckev_b((v16i8)vec3, (v16i8)vec3);

  out0 = __msa_copy_u_w((v4i32)res0, 0);
  out1 = __msa_copy_u_w((v4i32)res0, 1);
  out2 = __msa_copy_u_w((v4i32)res1, 0);
  out3 = __msa_copy_u_w((v4i32)res1, 1);

  STORE_WORD(dst, out0);
  dst += dst_stride;
  STORE_WORD(dst, out1);
  dst += dst_stride;
  STORE_WORD(dst, out2);
  dst += dst_stride;
  STORE_WORD(dst, out3);
}

static void common_hz_2t_4x8_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter) {
  uint32_t out0, out1, out2, out3;
  v16u8 filt0;
  v16i8 src0, src1, src2, src3, src4, src5, src6, src7, mask;
  v16u8 vec0, vec1, vec2, vec3;
  v8u16 vec4, vec5, vec6, vec7;
  v16i8 res0, res1, res2, res3;
  v8u16 filt, const255;

  mask = LOAD_SB(&mc_filt_mask_arr[16]);

  /* rearranging filter */
  filt = LOAD_UH(filter);
  filt0 = (v16u8)__msa_splati_h((v8i16)filt, 0);

  const255 = (v8u16)__msa_ldi_h(255);

  LOAD_8VECS_SB(src, src_stride,
                src0, src1, src2, src3, src4, src5, src6, src7);

  vec0 = (v16u8)__msa_vshf_b(mask, src1, src0);
  vec1 = (v16u8)__msa_vshf_b(mask, src3, src2);
  vec2 = (v16u8)__msa_vshf_b(mask, src5, src4);
  vec3 = (v16u8)__msa_vshf_b(mask, src7, src6);

  vec4 = __msa_dotp_u_h(vec0, filt0);
  vec5 = __msa_dotp_u_h(vec1, filt0);
  vec6 = __msa_dotp_u_h(vec2, filt0);
  vec7 = __msa_dotp_u_h(vec3, filt0);

  vec4 = (v8u16)__msa_srari_h((v8i16)vec4, FILTER_BITS);
  vec5 = (v8u16)__msa_srari_h((v8i16)vec5, FILTER_BITS);
  vec6 = (v8u16)__msa_srari_h((v8i16)vec6, FILTER_BITS);
  vec7 = (v8u16)__msa_srari_h((v8i16)vec7, FILTER_BITS);

  vec4 = __msa_min_u_h(vec4, const255);
  vec5 = __msa_min_u_h(vec5, const255);
  vec6 = __msa_min_u_h(vec6, const255);
  vec7 = __msa_min_u_h(vec7, const255);

  res0 = __msa_pckev_b((v16i8)vec4, (v16i8)vec4);
  res1 = __msa_pckev_b((v16i8)vec5, (v16i8)vec5);
  res2 = __msa_pckev_b((v16i8)vec6, (v16i8)vec6);
  res3 = __msa_pckev_b((v16i8)vec7, (v16i8)vec7);

  out0 = __msa_copy_u_w((v4i32)res0, 0);
  out1 = __msa_copy_u_w((v4i32)res0, 1);
  out2 = __msa_copy_u_w((v4i32)res1, 0);
  out3 = __msa_copy_u_w((v4i32)res1, 1);

  STORE_WORD(dst, out0);
  dst += dst_stride;
  STORE_WORD(dst, out1);
  dst += dst_stride;
  STORE_WORD(dst, out2);
  dst += dst_stride;
  STORE_WORD(dst, out3);
  dst += dst_stride;

  out0 = __msa_copy_u_w((v4i32)res2, 0);
  out1 = __msa_copy_u_w((v4i32)res2, 1);
  out2 = __msa_copy_u_w((v4i32)res3, 0);
  out3 = __msa_copy_u_w((v4i32)res3, 1);

  STORE_WORD(dst, out0);
  dst += dst_stride;
  STORE_WORD(dst, out1);
  dst += dst_stride;
  STORE_WORD(dst, out2);
  dst += dst_stride;
  STORE_WORD(dst, out3);
}

static void common_hz_2t_4w_msa(const uint8_t *src, int32_t src_stride,
                                uint8_t *dst, int32_t dst_stride,
                                int8_t *filter, int32_t height) {
  if (4 == height) {
    common_hz_2t_4x4_msa(src, src_stride, dst, dst_stride, filter);
  } else if (8 == height) {
    common_hz_2t_4x8_msa(src, src_stride, dst, dst_stride, filter);
  }
}

static void common_hz_2t_8x4_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter) {
  v16u8 filt0;
  v16i8 src0, src1, src2, src3, mask;
  v8u16 vec0, vec1, vec2, vec3;
  v8u16 out0, out1, out2, out3;
  v8u16 const255, filt;

  mask = LOAD_SB(&mc_filt_mask_arr[0]);

  /* rearranging filter */
  filt = LOAD_UH(filter);
  filt0 = (v16u8)__msa_splati_h((v8i16)filt, 0);

  const255 = (v8u16)__msa_ldi_h(255);

  LOAD_4VECS_SB(src, src_stride, src0, src1, src2, src3);

  vec0 = (v8u16)__msa_vshf_b(mask, src0, src0);
  vec1 = (v8u16)__msa_vshf_b(mask, src1, src1);
  vec2 = (v8u16)__msa_vshf_b(mask, src2, src2);
  vec3 = (v8u16)__msa_vshf_b(mask, src3, src3);

  vec0 = __msa_dotp_u_h((v16u8)vec0, filt0);
  vec1 = __msa_dotp_u_h((v16u8)vec1, filt0);
  vec2 = __msa_dotp_u_h((v16u8)vec2, filt0);
  vec3 = __msa_dotp_u_h((v16u8)vec3, filt0);

  SRARI_H_4VECS_UH(vec0, vec1, vec2, vec3, vec0, vec1, vec2, vec3, FILTER_BITS);

  out0 = __msa_min_u_h(vec0, const255);
  out1 = __msa_min_u_h(vec1, const255);
  out2 = __msa_min_u_h(vec2, const255);
  out3 = __msa_min_u_h(vec3, const255);

  PCKEV_B_STORE_8_BYTES_4(out0, out1, out2, out3, dst, dst_stride);
}

static void common_hz_2t_8x8mult_msa(const uint8_t *src, int32_t src_stride,
                                     uint8_t *dst, int32_t dst_stride,
                                     int8_t *filter, int32_t height) {
  v16u8 filt0;
  v16i8 src0, src1, src2, src3, mask;
  v8u16 vec0, vec1, vec2, vec3;
  v8u16 filt, const255;

  mask = LOAD_SB(&mc_filt_mask_arr[0]);

  /* rearranging filter */
  filt = LOAD_UH(filter);
  filt0 = (v16u8)__msa_splati_h((v8i16)filt, 0);

  const255 = (v8u16)__msa_ldi_h(255);

  LOAD_4VECS_SB(src, src_stride, src0, src1, src2, src3);
  src += (4 * src_stride);

  vec0 = (v8u16)__msa_vshf_b(mask, src0, src0);
  vec1 = (v8u16)__msa_vshf_b(mask, src1, src1);
  vec2 = (v8u16)__msa_vshf_b(mask, src2, src2);
  vec3 = (v8u16)__msa_vshf_b(mask, src3, src3);

  vec0 = __msa_dotp_u_h((v16u8)vec0, filt0);
  vec1 = __msa_dotp_u_h((v16u8)vec1, filt0);
  vec2 = __msa_dotp_u_h((v16u8)vec2, filt0);
  vec3 = __msa_dotp_u_h((v16u8)vec3, filt0);

  SRARI_H_4VECS_UH(vec0, vec1, vec2, vec3, vec0, vec1, vec2, vec3, FILTER_BITS);

  vec0 = __msa_min_u_h(vec0, const255);
  vec1 = __msa_min_u_h(vec1, const255);
  vec2 = __msa_min_u_h(vec2, const255);
  vec3 = __msa_min_u_h(vec3, const255);

  LOAD_4VECS_SB(src, src_stride, src0, src1, src2, src3);
  src += (4 * src_stride);

  PCKEV_B_STORE_8_BYTES_4(vec0, vec1, vec2, vec3, dst, dst_stride);
  dst += (4 * dst_stride);

  vec0 = (v8u16)__msa_vshf_b(mask, src0, src0);
  vec1 = (v8u16)__msa_vshf_b(mask, src1, src1);
  vec2 = (v8u16)__msa_vshf_b(mask, src2, src2);
  vec3 = (v8u16)__msa_vshf_b(mask, src3, src3);

  vec0 = __msa_dotp_u_h((v16u8)vec0, filt0);
  vec1 = __msa_dotp_u_h((v16u8)vec1, filt0);
  vec2 = __msa_dotp_u_h((v16u8)vec2, filt0);
  vec3 = __msa_dotp_u_h((v16u8)vec3, filt0);

  SRARI_H_4VECS_UH(vec0, vec1, vec2, vec3, vec0, vec1, vec2, vec3, FILTER_BITS);

  vec0 = __msa_min_u_h(vec0, const255);
  vec1 = __msa_min_u_h(vec1, const255);
  vec2 = __msa_min_u_h(vec2, const255);
  vec3 = __msa_min_u_h(vec3, const255);

  PCKEV_B_STORE_8_BYTES_4(vec0, vec1, vec2, vec3, dst, dst_stride);
  dst += (4 * dst_stride);

  if (16 == height) {
    LOAD_4VECS_SB(src, src_stride, src0, src1, src2, src3);
    src += (4 * src_stride);

    vec0 = (v8u16)__msa_vshf_b(mask, src0, src0);
    vec1 = (v8u16)__msa_vshf_b(mask, src1, src1);
    vec2 = (v8u16)__msa_vshf_b(mask, src2, src2);
    vec3 = (v8u16)__msa_vshf_b(mask, src3, src3);

    vec0 = __msa_dotp_u_h((v16u8)vec0, filt0);
    vec1 = __msa_dotp_u_h((v16u8)vec1, filt0);
    vec2 = __msa_dotp_u_h((v16u8)vec2, filt0);
    vec3 = __msa_dotp_u_h((v16u8)vec3, filt0);

    SRARI_H_4VECS_UH(vec0, vec1, vec2, vec3,
                     vec0, vec1, vec2, vec3, FILTER_BITS);

    vec0 = __msa_min_u_h(vec0, const255);
    vec1 = __msa_min_u_h(vec1, const255);
    vec2 = __msa_min_u_h(vec2, const255);
    vec3 = __msa_min_u_h(vec3, const255);

    LOAD_4VECS_SB(src, src_stride, src0, src1, src2, src3);
    src += (4 * src_stride);

    PCKEV_B_STORE_8_BYTES_4(vec0, vec1, vec2, vec3, dst, dst_stride);
    dst += (4 * dst_stride);

    vec0 = (v8u16)__msa_vshf_b(mask, src0, src0);
    vec1 = (v8u16)__msa_vshf_b(mask, src1, src1);
    vec2 = (v8u16)__msa_vshf_b(mask, src2, src2);
    vec3 = (v8u16)__msa_vshf_b(mask, src3, src3);

    vec0 = __msa_dotp_u_h((v16u8)vec0, filt0);
    vec1 = __msa_dotp_u_h((v16u8)vec1, filt0);
    vec2 = __msa_dotp_u_h((v16u8)vec2, filt0);
    vec3 = __msa_dotp_u_h((v16u8)vec3, filt0);

    SRARI_H_4VECS_UH(vec0, vec1, vec2, vec3,
                     vec0, vec1, vec2, vec3, FILTER_BITS);

    vec0 = __msa_min_u_h(vec0, const255);
    vec1 = __msa_min_u_h(vec1, const255);
    vec2 = __msa_min_u_h(vec2, const255);
    vec3 = __msa_min_u_h(vec3, const255);

    PCKEV_B_STORE_8_BYTES_4(vec0, vec1, vec2, vec3, dst, dst_stride);
  }
}

static void common_hz_2t_8w_msa(const uint8_t *src, int32_t src_stride,
                                uint8_t *dst, int32_t dst_stride,
                                int8_t *filter, int32_t height) {
  if (4 == height) {
    common_hz_2t_8x4_msa(src, src_stride, dst, dst_stride, filter);
  } else {
    common_hz_2t_8x8mult_msa(src, src_stride, dst, dst_stride, filter, height);
  }
}

static void common_hz_2t_16w_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter, int32_t height) {
  uint32_t loop_cnt;
  v16i8 src0, src1, src2, src3, src4, src5, src6, src7, mask;
  v16u8 filt0;
  v16u8 vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7;
  v8u16 out0, out1, out2, out3, out4, out5, out6, out7;
  v8u16 filt, const255;

  mask = LOAD_SB(&mc_filt_mask_arr[0]);

  loop_cnt = (height >> 2) - 1;

  /* rearranging filter */
  filt = LOAD_UH(filter);
  filt0 = (v16u8)__msa_splati_h((v8i16)filt, 0);

  const255 = (v8u16)__msa_ldi_h(255);

  src0 = LOAD_SB(src);
  src1 = LOAD_SB(src + 8);
  src += src_stride;
  src2 = LOAD_SB(src);
  src3 = LOAD_SB(src + 8);
  src += src_stride;
  src4 = LOAD_SB(src);
  src5 = LOAD_SB(src + 8);
  src += src_stride;
  src6 = LOAD_SB(src);
  src7 = LOAD_SB(src + 8);
  src += src_stride;

  vec0 = (v16u8)__msa_vshf_b(mask, src0, src0);
  vec1 = (v16u8)__msa_vshf_b(mask, src1, src1);
  vec2 = (v16u8)__msa_vshf_b(mask, src2, src2);
  vec3 = (v16u8)__msa_vshf_b(mask, src3, src3);
  vec4 = (v16u8)__msa_vshf_b(mask, src4, src4);
  vec5 = (v16u8)__msa_vshf_b(mask, src5, src5);
  vec6 = (v16u8)__msa_vshf_b(mask, src6, src6);
  vec7 = (v16u8)__msa_vshf_b(mask, src7, src7);

  out0 = __msa_dotp_u_h(vec0, filt0);
  out1 = __msa_dotp_u_h(vec1, filt0);
  out2 = __msa_dotp_u_h(vec2, filt0);
  out3 = __msa_dotp_u_h(vec3, filt0);
  out4 = __msa_dotp_u_h(vec4, filt0);
  out5 = __msa_dotp_u_h(vec5, filt0);
  out6 = __msa_dotp_u_h(vec6, filt0);
  out7 = __msa_dotp_u_h(vec7, filt0);

  out0 = (v8u16)__msa_srari_h((v8i16)out0, FILTER_BITS);
  out1 = (v8u16)__msa_srari_h((v8i16)out1, FILTER_BITS);
  out2 = (v8u16)__msa_srari_h((v8i16)out2, FILTER_BITS);
  out3 = (v8u16)__msa_srari_h((v8i16)out3, FILTER_BITS);
  out4 = (v8u16)__msa_srari_h((v8i16)out4, FILTER_BITS);
  out5 = (v8u16)__msa_srari_h((v8i16)out5, FILTER_BITS);
  out6 = (v8u16)__msa_srari_h((v8i16)out6, FILTER_BITS);
  out7 = (v8u16)__msa_srari_h((v8i16)out7, FILTER_BITS);

  out0 = __msa_min_u_h(out0, const255);
  out1 = __msa_min_u_h(out1, const255);
  out2 = __msa_min_u_h(out2, const255);
  out3 = __msa_min_u_h(out3, const255);
  out4 = __msa_min_u_h(out4, const255);
  out5 = __msa_min_u_h(out5, const255);
  out6 = __msa_min_u_h(out6, const255);
  out7 = __msa_min_u_h(out7, const255);

  PCKEV_B_STORE_VEC(out1, out0, dst);
  dst += dst_stride;
  PCKEV_B_STORE_VEC(out3, out2, dst);
  dst += dst_stride;
  PCKEV_B_STORE_VEC(out5, out4, dst);
  dst += dst_stride;
  PCKEV_B_STORE_VEC(out7, out6, dst);
  dst += dst_stride;

  for (; loop_cnt--;) {
    src0 = LOAD_SB(src);
    src1 = LOAD_SB(src + 8);
    src += src_stride;
    src2 = LOAD_SB(src);
    src3 = LOAD_SB(src + 8);
    src += src_stride;
    src4 = LOAD_SB(src);
    src5 = LOAD_SB(src + 8);
    src += src_stride;
    src6 = LOAD_SB(src);
    src7 = LOAD_SB(src + 8);
    src += src_stride;

    vec0 = (v16u8)__msa_vshf_b(mask, src0, src0);
    vec1 = (v16u8)__msa_vshf_b(mask, src1, src1);
    vec2 = (v16u8)__msa_vshf_b(mask, src2, src2);
    vec3 = (v16u8)__msa_vshf_b(mask, src3, src3);
    vec4 = (v16u8)__msa_vshf_b(mask, src4, src4);
    vec5 = (v16u8)__msa_vshf_b(mask, src5, src5);
    vec6 = (v16u8)__msa_vshf_b(mask, src6, src6);
    vec7 = (v16u8)__msa_vshf_b(mask, src7, src7);

    out0 = __msa_dotp_u_h(vec0, filt0);
    out1 = __msa_dotp_u_h(vec1, filt0);
    out2 = __msa_dotp_u_h(vec2, filt0);
    out3 = __msa_dotp_u_h(vec3, filt0);
    out4 = __msa_dotp_u_h(vec4, filt0);
    out5 = __msa_dotp_u_h(vec5, filt0);
    out6 = __msa_dotp_u_h(vec6, filt0);
    out7 = __msa_dotp_u_h(vec7, filt0);

    out0 = (v8u16)__msa_srari_h((v8i16)out0, FILTER_BITS);
    out1 = (v8u16)__msa_srari_h((v8i16)out1, FILTER_BITS);
    out2 = (v8u16)__msa_srari_h((v8i16)out2, FILTER_BITS);
    out3 = (v8u16)__msa_srari_h((v8i16)out3, FILTER_BITS);
    out4 = (v8u16)__msa_srari_h((v8i16)out4, FILTER_BITS);
    out5 = (v8u16)__msa_srari_h((v8i16)out5, FILTER_BITS);
    out6 = (v8u16)__msa_srari_h((v8i16)out6, FILTER_BITS);
    out7 = (v8u16)__msa_srari_h((v8i16)out7, FILTER_BITS);

    out0 = __msa_min_u_h(out0, const255);
    out1 = __msa_min_u_h(out1, const255);
    out2 = __msa_min_u_h(out2, const255);
    out3 = __msa_min_u_h(out3, const255);
    out4 = __msa_min_u_h(out4, const255);
    out5 = __msa_min_u_h(out5, const255);
    out6 = __msa_min_u_h(out6, const255);
    out7 = __msa_min_u_h(out7, const255);

    PCKEV_B_STORE_VEC(out1, out0, dst);
    dst += dst_stride;
    PCKEV_B_STORE_VEC(out3, out2, dst);
    dst += dst_stride;
    PCKEV_B_STORE_VEC(out5, out4, dst);
    dst += dst_stride;
    PCKEV_B_STORE_VEC(out7, out6, dst);
    dst += dst_stride;
  }
}

static void common_hz_2t_32w_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter, int32_t height) {
  uint32_t loop_cnt;
  v16i8 src0, src1, src2, src3, src4, src5, src6, src7, mask;
  v16u8 filt0;
  v16u8 vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7;
  v8u16 out0, out1, out2, out3, out4, out5, out6, out7;
  v8u16 filt, const255;

  mask = LOAD_SB(&mc_filt_mask_arr[0]);

  /* rearranging filter */
  filt = LOAD_UH(filter);
  filt0 = (v16u8)__msa_splati_h((v8i16)filt, 0);

  const255 = (v8u16)__msa_ldi_h(255);

  for (loop_cnt = height >> 1; loop_cnt--;) {
    src0 = LOAD_SB(src);
    src2 = LOAD_SB(src + 16);
    src3 = LOAD_SB(src + 24);
    src1 = __msa_sld_b(src2, src0, 8);
    src += src_stride;
    src4 = LOAD_SB(src);
    src6 = LOAD_SB(src + 16);
    src7 = LOAD_SB(src + 24);
    src5 = __msa_sld_b(src6, src4, 8);
    src += src_stride;

    vec0 = (v16u8)__msa_vshf_b(mask, src0, src0);
    vec1 = (v16u8)__msa_vshf_b(mask, src1, src1);
    vec2 = (v16u8)__msa_vshf_b(mask, src2, src2);
    vec3 = (v16u8)__msa_vshf_b(mask, src3, src3);
    vec4 = (v16u8)__msa_vshf_b(mask, src4, src4);
    vec5 = (v16u8)__msa_vshf_b(mask, src5, src5);
    vec6 = (v16u8)__msa_vshf_b(mask, src6, src6);
    vec7 = (v16u8)__msa_vshf_b(mask, src7, src7);

    out0 = __msa_dotp_u_h(vec0, filt0);
    out1 = __msa_dotp_u_h(vec1, filt0);
    out2 = __msa_dotp_u_h(vec2, filt0);
    out3 = __msa_dotp_u_h(vec3, filt0);
    out4 = __msa_dotp_u_h(vec4, filt0);
    out5 = __msa_dotp_u_h(vec5, filt0);
    out6 = __msa_dotp_u_h(vec6, filt0);
    out7 = __msa_dotp_u_h(vec7, filt0);

    out0 = (v8u16)__msa_srari_h((v8i16)out0, FILTER_BITS);
    out1 = (v8u16)__msa_srari_h((v8i16)out1, FILTER_BITS);
    out2 = (v8u16)__msa_srari_h((v8i16)out2, FILTER_BITS);
    out3 = (v8u16)__msa_srari_h((v8i16)out3, FILTER_BITS);
    out4 = (v8u16)__msa_srari_h((v8i16)out4, FILTER_BITS);
    out5 = (v8u16)__msa_srari_h((v8i16)out5, FILTER_BITS);
    out6 = (v8u16)__msa_srari_h((v8i16)out6, FILTER_BITS);
    out7 = (v8u16)__msa_srari_h((v8i16)out7, FILTER_BITS);

    out0 = __msa_min_u_h(out0, const255);
    out1 = __msa_min_u_h(out1, const255);
    out2 = __msa_min_u_h(out2, const255);
    out3 = __msa_min_u_h(out3, const255);
    out4 = __msa_min_u_h(out4, const255);
    out5 = __msa_min_u_h(out5, const255);
    out6 = __msa_min_u_h(out6, const255);
    out7 = __msa_min_u_h(out7, const255);

    PCKEV_B_STORE_VEC(out1, out0, dst);
    PCKEV_B_STORE_VEC(out3, out2, dst + 16);
    dst += dst_stride;
    PCKEV_B_STORE_VEC(out5, out4, dst);
    PCKEV_B_STORE_VEC(out7, out6, dst + 16);
    dst += dst_stride;
  }
}

static void common_hz_2t_64w_msa(const uint8_t *src, int32_t src_stride,
                                 uint8_t *dst, int32_t dst_stride,
                                 int8_t *filter, int32_t height) {
  uint32_t loop_cnt;
  v16i8 src0, src1, src2, src3, src4, src5, src6, src7, mask;
  v16u8 filt0;
  v16u8 vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7;
  v8u16 out0, out1, out2, out3, out4, out5, out6, out7;
  v8u16 filt, const255;

  mask = LOAD_SB(&mc_filt_mask_arr[0]);

  /* rearranging filter */
  filt = LOAD_UH(filter);
  filt0 = (v16u8)__msa_splati_h((v8i16)filt, 0);

  const255 = (v8u16)__msa_ldi_h(255);

  for (loop_cnt = height; loop_cnt--;) {
    src0 = LOAD_SB(src);
    src2 = LOAD_SB(src + 16);
    src4 = LOAD_SB(src + 32);
    src6 = LOAD_SB(src + 48);
    src7 = LOAD_SB(src + 56);
    src1 = __msa_sld_b(src2, src0, 8);
    src3 = __msa_sld_b(src4, src2, 8);
    src5 = __msa_sld_b(src6, src4, 8);
    src += src_stride;

    vec0 = (v16u8)__msa_vshf_b(mask, src0, src0);
    vec1 = (v16u8)__msa_vshf_b(mask, src1, src1);
    vec2 = (v16u8)__msa_vshf_b(mask, src2, src2);
    vec3 = (v16u8)__msa_vshf_b(mask, src3, src3);
    vec4 = (v16u8)__msa_vshf_b(mask, src4, src4);
    vec5 = (v16u8)__msa_vshf_b(mask, src5, src5);
    vec6 = (v16u8)__msa_vshf_b(mask, src6, src6);
    vec7 = (v16u8)__msa_vshf_b(mask, src7, src7);

    out0 = __msa_dotp_u_h(vec0, filt0);
    out1 = __msa_dotp_u_h(vec1, filt0);
    out2 = __msa_dotp_u_h(vec2, filt0);
    out3 = __msa_dotp_u_h(vec3, filt0);
    out4 = __msa_dotp_u_h(vec4, filt0);
    out5 = __msa_dotp_u_h(vec5, filt0);
    out6 = __msa_dotp_u_h(vec6, filt0);
    out7 = __msa_dotp_u_h(vec7, filt0);

    out0 = (v8u16)__msa_srari_h((v8i16)out0, FILTER_BITS);
    out1 = (v8u16)__msa_srari_h((v8i16)out1, FILTER_BITS);
    out2 = (v8u16)__msa_srari_h((v8i16)out2, FILTER_BITS);
    out3 = (v8u16)__msa_srari_h((v8i16)out3, FILTER_BITS);
    out4 = (v8u16)__msa_srari_h((v8i16)out4, FILTER_BITS);
    out5 = (v8u16)__msa_srari_h((v8i16)out5, FILTER_BITS);
    out6 = (v8u16)__msa_srari_h((v8i16)out6, FILTER_BITS);
    out7 = (v8u16)__msa_srari_h((v8i16)out7, FILTER_BITS);

    out0 = __msa_min_u_h(out0, const255);
    out1 = __msa_min_u_h(out1, const255);
    out2 = __msa_min_u_h(out2, const255);
    out3 = __msa_min_u_h(out3, const255);
    out4 = __msa_min_u_h(out4, const255);
    out5 = __msa_min_u_h(out5, const255);
    out6 = __msa_min_u_h(out6, const255);
    out7 = __msa_min_u_h(out7, const255);

    PCKEV_B_STORE_VEC(out1, out0, dst);
    PCKEV_B_STORE_VEC(out3, out2, dst + 16);
    PCKEV_B_STORE_VEC(out5, out4, dst + 32);
    PCKEV_B_STORE_VEC(out7, out6, dst + 48);
    dst += dst_stride;
  }
}

void vp9_convolve8_horiz_msa(const uint8_t *src, ptrdiff_t src_stride,
                             uint8_t *dst, ptrdiff_t dst_stride,
                             const int16_t *filter_x, int x_step_q4,
                             const int16_t *filter_y, int y_step_q4,
                             int w, int h) {
  int8_t cnt, filt_hor[8];

  if (16 != x_step_q4) {
    vp9_convolve8_horiz_c(src, src_stride, dst, dst_stride,
                          filter_x, x_step_q4, filter_y, y_step_q4,
                          w, h);
    return;
  }

  if (((const int32_t *)filter_x)[1] == 0x800000) {
    vp9_convolve_copy(src, src_stride, dst, dst_stride,
                      filter_x, x_step_q4, filter_y, y_step_q4,
                      w, h);
    return;
  }

  for (cnt = 0; cnt < 8; ++cnt) {
    filt_hor[cnt] = filter_x[cnt];
  }

  if (((const int32_t *)filter_x)[0] == 0) {
    switch (w) {
      case 4:
        common_hz_2t_4w_msa(src, (int32_t)src_stride,
                            dst, (int32_t)dst_stride,
                            &filt_hor[3], h);
        break;
      case 8:
        common_hz_2t_8w_msa(src, (int32_t)src_stride,
                            dst, (int32_t)dst_stride,
                            &filt_hor[3], h);
        break;
      case 16:
        common_hz_2t_16w_msa(src, (int32_t)src_stride,
                             dst, (int32_t)dst_stride,
                             &filt_hor[3], h);
        break;
      case 32:
        common_hz_2t_32w_msa(src, (int32_t)src_stride,
                             dst, (int32_t)dst_stride,
                             &filt_hor[3], h);
        break;
      case 64:
        common_hz_2t_64w_msa(src, (int32_t)src_stride,
                             dst, (int32_t)dst_stride,
                             &filt_hor[3], h);
        break;
      default:
        vp9_convolve8_horiz_c(src, src_stride, dst, dst_stride,
                              filter_x, x_step_q4, filter_y, y_step_q4,
                              w, h);
        break;
    }
  } else {
    switch (w) {
      case 4:
        common_hz_8t_4w_msa(src, (int32_t)src_stride,
                            dst, (int32_t)dst_stride,
                            filt_hor, h);
        break;
      case 8:
        common_hz_8t_8w_msa(src, (int32_t)src_stride,
                            dst, (int32_t)dst_stride,
                            filt_hor, h);
        break;
      case 16:
        common_hz_8t_16w_msa(src, (int32_t)src_stride,
                             dst, (int32_t)dst_stride,
                             filt_hor, h);
        break;
      case 32:
        common_hz_8t_32w_msa(src, (int32_t)src_stride,
                             dst, (int32_t)dst_stride,
                             filt_hor, h);
        break;
      case 64:
        common_hz_8t_64w_msa(src, (int32_t)src_stride,
                             dst, (int32_t)dst_stride,
                             filt_hor, h);
        break;
      default:
        vp9_convolve8_horiz_c(src, src_stride, dst, dst_stride,
                              filter_x, x_step_q4, filter_y, y_step_q4,
                              w, h);
        break;
    }
  }
}
