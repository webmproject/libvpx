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

const uint8_t mc_filt_mask_arr[16 * 3] = {
  /* 8 width cases */
  0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8,
  /* 4 width cases */
  0, 1, 1, 2, 2, 3, 3, 4, 16, 17, 17, 18, 18, 19, 19, 20,
  /* 4 width cases */
  8, 9, 9, 10, 10, 11, 11, 12, 24, 25, 25, 26, 26, 27, 27, 28
};

static void common_hv_8ht_8vt_4w_msa(const uint8_t *src, int32_t src_stride,
                                     uint8_t *dst, int32_t dst_stride,
                                     int8_t *filter_horiz, int8_t *filter_vert,
                                     int32_t height) {
  uint32_t loop_cnt;
  v16i8 src0, src1, src2, src3, src4, src5, src6, src7, src8, src9, src10;
  v16i8 filt_horiz0, filt_horiz1, filt_horiz2, filt_horiz3;
  v16u8 mask0, mask1, mask2, mask3;
  v8i16 filt_horiz;
  v8i16 horiz_out0, horiz_out1, horiz_out2, horiz_out3, horiz_out4;
  v8i16 horiz_out5, horiz_out6, horiz_out7, horiz_out8, horiz_out9;
  v8i16 tmp0, tmp1, out0, out1, out2, out3, out4;
  v8i16 filt, filt_vert0, filt_vert1, filt_vert2, filt_vert3;

  mask0 = LOAD_UB(&mc_filt_mask_arr[16]);

  src -= (3 + 3 * src_stride);

  /* rearranging filter */
  filt_horiz = LOAD_SH(filter_horiz);
  filt_horiz0 = (v16i8)__msa_splati_h(filt_horiz, 0);
  filt_horiz1 = (v16i8)__msa_splati_h(filt_horiz, 1);
  filt_horiz2 = (v16i8)__msa_splati_h(filt_horiz, 2);
  filt_horiz3 = (v16i8)__msa_splati_h(filt_horiz, 3);

  mask1 = mask0 + 2;
  mask2 = mask0 + 4;
  mask3 = mask0 + 6;

  LOAD_7VECS_SB(src, src_stride, src0, src1, src2, src3, src4, src5, src6);
  src += (7 * src_stride);

  XORI_B_7VECS_SB(src0, src1, src2, src3, src4, src5, src6,
                  src0, src1, src2, src3, src4, src5, src6, 128);

  horiz_out0 = HORIZ_8TAP_FILT_2VECS(src0, src1, mask0, mask1, mask2, mask3,
                                     filt_horiz0, filt_horiz1, filt_horiz2,
                                     filt_horiz3);
  horiz_out2 = HORIZ_8TAP_FILT_2VECS(src2, src3, mask0, mask1, mask2, mask3,
                                     filt_horiz0, filt_horiz1, filt_horiz2,
                                     filt_horiz3);
  horiz_out4 = HORIZ_8TAP_FILT_2VECS(src4, src5, mask0, mask1, mask2, mask3,
                                     filt_horiz0, filt_horiz1, filt_horiz2,
                                     filt_horiz3);
  horiz_out5 = HORIZ_8TAP_FILT_2VECS(src5, src6, mask0, mask1, mask2, mask3,
                                     filt_horiz0, filt_horiz1, filt_horiz2,
                                     filt_horiz3);
  horiz_out1 = (v8i16)__msa_sldi_b((v16i8)horiz_out2, (v16i8)horiz_out0, 8);
  horiz_out3 = (v8i16)__msa_sldi_b((v16i8)horiz_out4, (v16i8)horiz_out2, 8);

  filt = LOAD_SH(filter_vert);
  filt_vert0 = __msa_splati_h(filt, 0);
  filt_vert1 = __msa_splati_h(filt, 1);
  filt_vert2 = __msa_splati_h(filt, 2);
  filt_vert3 = __msa_splati_h(filt, 3);

  out0 = (v8i16)__msa_ilvev_b((v16i8)horiz_out1, (v16i8)horiz_out0);
  out1 = (v8i16)__msa_ilvev_b((v16i8)horiz_out3, (v16i8)horiz_out2);
  out2 = (v8i16)__msa_ilvev_b((v16i8)horiz_out5, (v16i8)horiz_out4);

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LOAD_4VECS_SB(src, src_stride, src7, src8, src9, src10);
    src += (4 * src_stride);

    XORI_B_4VECS_SB(src7, src8, src9, src10, src7, src8, src9, src10, 128);

    horiz_out7 = HORIZ_8TAP_FILT_2VECS(src7, src8, mask0, mask1, mask2, mask3,
                                       filt_horiz0, filt_horiz1, filt_horiz2,
                                       filt_horiz3);
    horiz_out6 = (v8i16)__msa_sldi_b((v16i8)horiz_out7, (v16i8)horiz_out5, 8);

    out3 = (v8i16)__msa_ilvev_b((v16i8)horiz_out7, (v16i8)horiz_out6);

    tmp0 = FILT_8TAP_DPADD_S_H(out0, out1, out2, out3, filt_vert0, filt_vert1,
                               filt_vert2, filt_vert3);

    horiz_out9 = HORIZ_8TAP_FILT_2VECS(src9, src10, mask0, mask1, mask2, mask3,
                                       filt_horiz0, filt_horiz1, filt_horiz2,
                                       filt_horiz3);
    horiz_out8 = (v8i16)__msa_sldi_b((v16i8)horiz_out9, (v16i8)horiz_out7, 8);

    out4 = (v8i16)__msa_ilvev_b((v16i8)horiz_out9, (v16i8)horiz_out8);

    tmp1 = FILT_8TAP_DPADD_S_H(out1, out2, out3, out4, filt_vert0, filt_vert1,
                               filt_vert2, filt_vert3);
    tmp0 = SRARI_SATURATE_SIGNED_H(tmp0, FILTER_BITS, 7);
    tmp1 = SRARI_SATURATE_SIGNED_H(tmp1, FILTER_BITS, 7);

    PCKEV_2B_XORI128_STORE_4_BYTES_4(tmp0, tmp1, dst, dst_stride);
    dst += (4 * dst_stride);

    horiz_out5 = horiz_out9;

    out0 = out2;
    out1 = out3;
    out2 = out4;
  }
}

static void common_hv_8ht_8vt_8w_msa(const uint8_t *src, int32_t src_stride,
                                     uint8_t *dst, int32_t dst_stride,
                                     int8_t *filter_horiz, int8_t *filter_vert,
                                     int32_t height) {
  uint32_t loop_cnt;
  v16i8 src0, src1, src2, src3, src4, src5, src6, src7, src8, src9, src10;
  v16i8 filt_horiz0, filt_horiz1, filt_horiz2, filt_horiz3;
  v8i16 filt_horiz, filt, filt_vert0, filt_vert1, filt_vert2, filt_vert3;
  v16u8 mask0, mask1, mask2, mask3;
  v8i16 horiz_out0, horiz_out1, horiz_out2, horiz_out3;
  v8i16 horiz_out4, horiz_out5, horiz_out6, horiz_out7;
  v8i16 horiz_out8, horiz_out9, horiz_out10;
  v8i16 out0, out1, out2, out3, out4, out5, out6, out7, out8, out9;
  v8i16 tmp0, tmp1, tmp2, tmp3;

  mask0 = LOAD_UB(&mc_filt_mask_arr[0]);

  src -= (3 + 3 * src_stride);

  /* rearranging filter */
  filt_horiz = LOAD_SH(filter_horiz);
  filt_horiz0 = (v16i8)__msa_splati_h(filt_horiz, 0);
  filt_horiz1 = (v16i8)__msa_splati_h(filt_horiz, 1);
  filt_horiz2 = (v16i8)__msa_splati_h(filt_horiz, 2);
  filt_horiz3 = (v16i8)__msa_splati_h(filt_horiz, 3);

  mask1 = mask0 + 2;
  mask2 = mask0 + 4;
  mask3 = mask0 + 6;

  LOAD_7VECS_SB(src, src_stride, src0, src1, src2, src3, src4, src5, src6);
  src += (7 * src_stride);

  XORI_B_7VECS_SB(src0, src1, src2, src3, src4, src5, src6,
                  src0, src1, src2, src3, src4, src5, src6, 128);

  horiz_out0 = HORIZ_8TAP_FILT(src0, mask0, mask1, mask2, mask3, filt_horiz0,
                               filt_horiz1, filt_horiz2, filt_horiz3);
  horiz_out1 = HORIZ_8TAP_FILT(src1, mask0, mask1, mask2, mask3, filt_horiz0,
                               filt_horiz1, filt_horiz2, filt_horiz3);
  horiz_out2 = HORIZ_8TAP_FILT(src2, mask0, mask1, mask2, mask3, filt_horiz0,
                               filt_horiz1, filt_horiz2, filt_horiz3);
  horiz_out3 = HORIZ_8TAP_FILT(src3, mask0, mask1, mask2, mask3, filt_horiz0,
                               filt_horiz1, filt_horiz2, filt_horiz3);
  horiz_out4 = HORIZ_8TAP_FILT(src4, mask0, mask1, mask2, mask3, filt_horiz0,
                               filt_horiz1, filt_horiz2, filt_horiz3);
  horiz_out5 = HORIZ_8TAP_FILT(src5, mask0, mask1, mask2, mask3, filt_horiz0,
                               filt_horiz1, filt_horiz2, filt_horiz3);
  horiz_out6 = HORIZ_8TAP_FILT(src6, mask0, mask1, mask2, mask3, filt_horiz0,
                               filt_horiz1, filt_horiz2, filt_horiz3);

  filt = LOAD_SH(filter_vert);
  filt_vert0 = __msa_splati_h(filt, 0);
  filt_vert1 = __msa_splati_h(filt, 1);
  filt_vert2 = __msa_splati_h(filt, 2);
  filt_vert3 = __msa_splati_h(filt, 3);

  out0 = (v8i16)__msa_ilvev_b((v16i8)horiz_out1, (v16i8)horiz_out0);
  out1 = (v8i16)__msa_ilvev_b((v16i8)horiz_out3, (v16i8)horiz_out2);
  out2 = (v8i16)__msa_ilvev_b((v16i8)horiz_out5, (v16i8)horiz_out4);
  out4 = (v8i16)__msa_ilvev_b((v16i8)horiz_out2, (v16i8)horiz_out1);
  out5 = (v8i16)__msa_ilvev_b((v16i8)horiz_out4, (v16i8)horiz_out3);
  out6 = (v8i16)__msa_ilvev_b((v16i8)horiz_out6, (v16i8)horiz_out5);

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LOAD_4VECS_SB(src, src_stride, src7, src8, src9, src10);
    src += (4 * src_stride);

    XORI_B_4VECS_SB(src7, src8, src9, src10, src7, src8, src9, src10, 128);

    horiz_out7 = HORIZ_8TAP_FILT(src7, mask0, mask1, mask2, mask3, filt_horiz0,
                                 filt_horiz1, filt_horiz2, filt_horiz3);

    out3 = (v8i16)__msa_ilvev_b((v16i8)horiz_out7, (v16i8)horiz_out6);
    tmp0 = FILT_8TAP_DPADD_S_H(out0, out1, out2, out3, filt_vert0, filt_vert1,
                               filt_vert2, filt_vert3);
    tmp0 = SRARI_SATURATE_SIGNED_H(tmp0, FILTER_BITS, 7);

    horiz_out8 = HORIZ_8TAP_FILT(src8, mask0, mask1, mask2, mask3, filt_horiz0,
                                 filt_horiz1, filt_horiz2, filt_horiz3);

    out7 = (v8i16)__msa_ilvev_b((v16i8)horiz_out8, (v16i8)horiz_out7);
    tmp1 = FILT_8TAP_DPADD_S_H(out4, out5, out6, out7, filt_vert0, filt_vert1,
                               filt_vert2, filt_vert3);
    tmp1 = SRARI_SATURATE_SIGNED_H(tmp1, FILTER_BITS, 7);

    horiz_out9 = HORIZ_8TAP_FILT(src9, mask0, mask1, mask2, mask3, filt_horiz0,
                                 filt_horiz1, filt_horiz2, filt_horiz3);

    out8 = (v8i16)__msa_ilvev_b((v16i8)horiz_out9, (v16i8)horiz_out8);
    tmp2 = FILT_8TAP_DPADD_S_H(out1, out2, out3, out8, filt_vert0, filt_vert1,
                               filt_vert2, filt_vert3);
    tmp2 = SRARI_SATURATE_SIGNED_H(tmp2, FILTER_BITS, 7);

    horiz_out10 = HORIZ_8TAP_FILT(src10, mask0, mask1, mask2, mask3,
                                  filt_horiz0, filt_horiz1, filt_horiz2,
                                  filt_horiz3);

    out9 = (v8i16)__msa_ilvev_b((v16i8)horiz_out10, (v16i8)horiz_out9);
    tmp3 = FILT_8TAP_DPADD_S_H(out5, out6, out7, out9, filt_vert0, filt_vert1,
                               filt_vert2, filt_vert3);
    tmp3 = SRARI_SATURATE_SIGNED_H(tmp3, FILTER_BITS, 7);

    PCKEV_B_4_XORI128_STORE_8_BYTES_4(tmp0, tmp1, tmp2, tmp3, dst, dst_stride);
    dst += (4 * dst_stride);

    horiz_out6 = horiz_out10;

    out0 = out2;
    out1 = out3;
    out2 = out8;
    out4 = out6;
    out5 = out7;
    out6 = out9;
  }
}

static void common_hv_8ht_8vt_16w_msa(const uint8_t *src, int32_t src_stride,
                                      uint8_t *dst, int32_t dst_stride,
                                      int8_t *filter_horiz, int8_t *filter_vert,
                                      int32_t height) {
  int32_t multiple8_cnt;
  for (multiple8_cnt = 2; multiple8_cnt--;) {
    common_hv_8ht_8vt_8w_msa(src, src_stride, dst, dst_stride, filter_horiz,
                             filter_vert, height);
    src += 8;
    dst += 8;
  }
}

static void common_hv_8ht_8vt_32w_msa(const uint8_t *src, int32_t src_stride,
                                      uint8_t *dst, int32_t dst_stride,
                                      int8_t *filter_horiz, int8_t *filter_vert,
                                      int32_t height) {
  int32_t multiple8_cnt;
  for (multiple8_cnt = 4; multiple8_cnt--;) {
    common_hv_8ht_8vt_8w_msa(src, src_stride, dst, dst_stride, filter_horiz,
                             filter_vert, height);
    src += 8;
    dst += 8;
  }
}

static void common_hv_8ht_8vt_64w_msa(const uint8_t *src, int32_t src_stride,
                                      uint8_t *dst, int32_t dst_stride,
                                      int8_t *filter_horiz, int8_t *filter_vert,
                                      int32_t height) {
  int32_t multiple8_cnt;
  for (multiple8_cnt = 8; multiple8_cnt--;) {
    common_hv_8ht_8vt_8w_msa(src, src_stride, dst, dst_stride, filter_horiz,
                             filter_vert, height);
    src += 8;
    dst += 8;
  }
}

static void common_hv_2ht_2vt_4x4_msa(const uint8_t *src, int32_t src_stride,
                                      uint8_t *dst, int32_t dst_stride,
                                      int8_t *filter_horiz,
                                      int8_t *filter_vert) {
  uint32_t out0, out1, out2, out3;
  v16i8 src0, src1, src2, src3, src4, mask;
  v16u8 res0, res1, horiz_vec;
  v16u8 filt_vert, filt_horiz, vec0, vec1;
  v8u16 filt, tmp0, tmp1;
  v8u16 horiz_out0, horiz_out1, horiz_out2, horiz_out3, horiz_out4;

  mask = LOAD_SB(&mc_filt_mask_arr[16]);

  /* rearranging filter */
  filt = LOAD_UH(filter_horiz);
  filt_horiz = (v16u8)__msa_splati_h((v8i16)filt, 0);

  filt = LOAD_UH(filter_vert);
  filt_vert = (v16u8)__msa_splati_h((v8i16)filt, 0);

  LOAD_5VECS_SB(src, src_stride, src0, src1, src2, src3, src4);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src1, src0);
  horiz_out0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out0 = SRARI_SATURATE_UNSIGNED_H(horiz_out0, FILTER_BITS, 7);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src3, src2);
  horiz_out2 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out2 = SRARI_SATURATE_UNSIGNED_H(horiz_out2, FILTER_BITS, 7);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src4, src4);
  horiz_out4 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out4 = SRARI_SATURATE_UNSIGNED_H(horiz_out4, FILTER_BITS, 7);

  horiz_out1 = (v8u16)__msa_sldi_b((v16i8)horiz_out2, (v16i8)horiz_out0, 8);
  horiz_out3 = (v8u16)__msa_pckod_d((v2i64)horiz_out4, (v2i64)horiz_out2);

  vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out1, (v16i8)horiz_out0);
  vec1 = (v16u8)__msa_ilvev_b((v16i8)horiz_out3, (v16i8)horiz_out2);

  tmp0 = __msa_dotp_u_h(vec0, filt_vert);
  tmp1 = __msa_dotp_u_h(vec1, filt_vert);
  tmp0 = SRARI_SATURATE_UNSIGNED_H(tmp0, FILTER_BITS, 7);
  tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);

  res0 = (v16u8)__msa_pckev_b((v16i8)tmp0, (v16i8)tmp0);
  res1 = (v16u8)__msa_pckev_b((v16i8)tmp1, (v16i8)tmp1);

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

static void common_hv_2ht_2vt_4x8_msa(const uint8_t *src, int32_t src_stride,
                                      uint8_t *dst, int32_t dst_stride,
                                      int8_t *filter_horiz,
                                      int8_t *filter_vert) {
  uint32_t out0, out1, out2, out3;
  v16i8 src0, src1, src2, src3, src4, src5, src6, src7, src8, mask;
  v16u8 filt_horiz, filt_vert, horiz_vec;
  v16u8 vec0, vec1, vec2, vec3;
  v8u16 horiz_out0, horiz_out1, horiz_out2, horiz_out3;
  v8u16 vec4, vec5, vec6, vec7, filt;
  v8u16 horiz_out4, horiz_out5, horiz_out6, horiz_out7, horiz_out8;
  v16i8 res0, res1, res2, res3;

  mask = LOAD_SB(&mc_filt_mask_arr[16]);

  /* rearranging filter */
  filt = LOAD_UH(filter_horiz);
  filt_horiz = (v16u8)__msa_splati_h((v8i16)filt, 0);

  filt = LOAD_UH(filter_vert);
  filt_vert = (v16u8)__msa_splati_h((v8i16)filt, 0);

  LOAD_8VECS_SB(src, src_stride,
                src0, src1, src2, src3, src4, src5, src6, src7);
  src += (8 * src_stride);
  src8 = LOAD_SB(src);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src1, src0);
  horiz_out0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out0 = SRARI_SATURATE_UNSIGNED_H(horiz_out0, FILTER_BITS, 7);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src3, src2);
  horiz_out2 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out2 = SRARI_SATURATE_UNSIGNED_H(horiz_out2, FILTER_BITS, 7);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src5, src4);
  horiz_out4 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out4 = SRARI_SATURATE_UNSIGNED_H(horiz_out4, FILTER_BITS, 7);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src7, src6);
  horiz_out6 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out6 = SRARI_SATURATE_UNSIGNED_H(horiz_out6, FILTER_BITS, 7);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src8, src8);
  horiz_out8 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out8 = SRARI_SATURATE_UNSIGNED_H(horiz_out8, FILTER_BITS, 7);

  horiz_out1 = (v8u16)__msa_sldi_b((v16i8)horiz_out2, (v16i8)horiz_out0, 8);
  horiz_out3 = (v8u16)__msa_sldi_b((v16i8)horiz_out4, (v16i8)horiz_out2, 8);
  horiz_out5 = (v8u16)__msa_sldi_b((v16i8)horiz_out6, (v16i8)horiz_out4, 8);
  horiz_out7 = (v8u16)__msa_pckod_d((v2i64)horiz_out8, (v2i64)horiz_out6);

  vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out1, (v16i8)horiz_out0);
  vec1 = (v16u8)__msa_ilvev_b((v16i8)horiz_out3, (v16i8)horiz_out2);
  vec2 = (v16u8)__msa_ilvev_b((v16i8)horiz_out5, (v16i8)horiz_out4);
  vec3 = (v16u8)__msa_ilvev_b((v16i8)horiz_out7, (v16i8)horiz_out6);

  vec4 = __msa_dotp_u_h(vec0, filt_vert);
  vec5 = __msa_dotp_u_h(vec1, filt_vert);
  vec6 = __msa_dotp_u_h(vec2, filt_vert);
  vec7 = __msa_dotp_u_h(vec3, filt_vert);

  vec4 = SRARI_SATURATE_UNSIGNED_H(vec4, FILTER_BITS, 7);
  vec5 = SRARI_SATURATE_UNSIGNED_H(vec5, FILTER_BITS, 7);
  vec6 = SRARI_SATURATE_UNSIGNED_H(vec6, FILTER_BITS, 7);
  vec7 = SRARI_SATURATE_UNSIGNED_H(vec7, FILTER_BITS, 7);

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

static void common_hv_2ht_2vt_4w_msa(const uint8_t *src, int32_t src_stride,
                                     uint8_t *dst, int32_t dst_stride,
                                     int8_t *filter_horiz,
                                     int8_t *filter_vert,
                                     int32_t height) {
  if (4 == height) {
    common_hv_2ht_2vt_4x4_msa(src, src_stride, dst, dst_stride,
                              filter_horiz, filter_vert);
  } else if (8 == height) {
    common_hv_2ht_2vt_4x8_msa(src, src_stride, dst, dst_stride,
                              filter_horiz, filter_vert);
  }
}

static void common_hv_2ht_2vt_8x4_msa(const uint8_t *src, int32_t src_stride,
                                      uint8_t *dst, int32_t dst_stride,
                                      int8_t *filter_horiz,
                                      int8_t *filter_vert) {
  v16i8 src0, src1, src2, src3, src4, mask;
  v16u8 filt_horiz, filt_vert, horiz_vec;
  v16u8 vec0, vec1, vec2, vec3;
  v8u16 horiz_out0, horiz_out1;
  v8u16 tmp0, tmp1, tmp2, tmp3;
  v8i16 filt;

  mask = LOAD_SB(&mc_filt_mask_arr[0]);

  /* rearranging filter */
  filt = LOAD_SH(filter_horiz);
  filt_horiz = (v16u8)__msa_splati_h(filt, 0);

  filt = LOAD_SH(filter_vert);
  filt_vert = (v16u8)__msa_splati_h(filt, 0);

  LOAD_5VECS_SB(src, src_stride, src0, src1, src2, src3, src4);
  src += (5 * src_stride);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src0, src0);
  horiz_out0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out0 = SRARI_SATURATE_UNSIGNED_H(horiz_out0, FILTER_BITS, 7);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src1, src1);
  horiz_out1 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out1 = SRARI_SATURATE_UNSIGNED_H(horiz_out1, FILTER_BITS, 7);

  vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out1, (v16i8)horiz_out0);
  tmp0 = __msa_dotp_u_h(vec0, filt_vert);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src2, src2);
  horiz_out0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out0 = SRARI_SATURATE_UNSIGNED_H(horiz_out0, FILTER_BITS, 7);

  vec1 = (v16u8)__msa_ilvev_b((v16i8)horiz_out0, (v16i8)horiz_out1);
  tmp1 = __msa_dotp_u_h(vec1, filt_vert);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src3, src3);
  horiz_out1 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out1 = SRARI_SATURATE_UNSIGNED_H(horiz_out1, FILTER_BITS, 7);

  vec2 = (v16u8)__msa_ilvev_b((v16i8)horiz_out1, (v16i8)horiz_out0);
  tmp2 = __msa_dotp_u_h(vec2, filt_vert);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src4, src4);
  horiz_out0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out0 = SRARI_SATURATE_UNSIGNED_H(horiz_out0, FILTER_BITS, 7);

  vec3 = (v16u8)__msa_ilvev_b((v16i8)horiz_out0, (v16i8)horiz_out1);
  tmp3 = __msa_dotp_u_h(vec3, filt_vert);

  tmp0 = SRARI_SATURATE_UNSIGNED_H(tmp0, FILTER_BITS, 7);
  tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);
  tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);
  tmp3 = SRARI_SATURATE_UNSIGNED_H(tmp3, FILTER_BITS, 7);

  PCKEV_B_STORE_8_BYTES_4(tmp0, tmp1, tmp2, tmp3, dst, dst_stride);
}

static void common_hv_2ht_2vt_8x8mult_msa(const uint8_t *src,
                                          int32_t src_stride,
                                          uint8_t *dst,
                                          int32_t dst_stride,
                                          int8_t *filter_horiz,
                                          int8_t *filter_vert,
                                          int32_t height) {
  uint32_t loop_cnt;
  v16i8 src0, src1, src2, src3, src4, mask;
  v16u8 filt_horiz, filt_vert, vec0, horiz_vec;
  v8u16 horiz_out0, horiz_out1;
  v8u16 tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  v8i16 filt;

  mask = LOAD_SB(&mc_filt_mask_arr[0]);

  /* rearranging filter */
  filt = LOAD_SH(filter_horiz);
  filt_horiz = (v16u8)__msa_splati_h(filt, 0);

  filt = LOAD_SH(filter_vert);
  filt_vert = (v16u8)__msa_splati_h(filt, 0);

  src0 = LOAD_SB(src);
  src += src_stride;

  horiz_vec = (v16u8)__msa_vshf_b(mask, src0, src0);
  horiz_out0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out0 = SRARI_SATURATE_UNSIGNED_H(horiz_out0, FILTER_BITS, 7);

  for (loop_cnt = (height >> 3); loop_cnt--;) {
    LOAD_4VECS_SB(src, src_stride, src1, src2, src3, src4);
    src += (4 * src_stride);

    horiz_vec = (v16u8)__msa_vshf_b(mask, src1, src1);
    horiz_out1 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out1 = SRARI_SATURATE_UNSIGNED_H(horiz_out1, FILTER_BITS, 7);

    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out1, (v16i8)horiz_out0);
    tmp1 = __msa_dotp_u_h(vec0, filt_vert);

    horiz_vec = (v16u8)__msa_vshf_b(mask, src2, src2);
    horiz_out0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out0 = SRARI_SATURATE_UNSIGNED_H(horiz_out0, FILTER_BITS, 7);

    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out0, (v16i8)horiz_out1);
    tmp2 = (v8u16)__msa_dotp_u_h(vec0, filt_vert);

    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);
    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);

    horiz_vec = (v16u8)__msa_vshf_b(mask, src3, src3);
    horiz_out1 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out1 = SRARI_SATURATE_UNSIGNED_H(horiz_out1, FILTER_BITS, 7);

    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out1, (v16i8)horiz_out0);
    tmp3 = __msa_dotp_u_h(vec0, filt_vert);

    horiz_vec = (v16u8)__msa_vshf_b(mask, src4, src4);
    horiz_out0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out0 = SRARI_SATURATE_UNSIGNED_H(horiz_out0, FILTER_BITS, 7);

    LOAD_4VECS_SB(src, src_stride, src1, src2, src3, src4);
    src += (4 * src_stride);

    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out0, (v16i8)horiz_out1);
    tmp4 = __msa_dotp_u_h(vec0, filt_vert);

    tmp3 = SRARI_SATURATE_UNSIGNED_H(tmp3, FILTER_BITS, 7);
    tmp4 = SRARI_SATURATE_UNSIGNED_H(tmp4, FILTER_BITS, 7);

    PCKEV_B_STORE_8_BYTES_4(tmp1, tmp2, tmp3, tmp4, dst, dst_stride);
    dst += (4 * dst_stride);

    horiz_vec = (v16u8)__msa_vshf_b(mask, src1, src1);
    horiz_out1 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out1 = SRARI_SATURATE_UNSIGNED_H(horiz_out1, FILTER_BITS, 7);

    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out1, (v16i8)horiz_out0);
    tmp5 = __msa_dotp_u_h(vec0, filt_vert);

    horiz_vec = (v16u8)__msa_vshf_b(mask, src2, src2);
    horiz_out0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out0 = SRARI_SATURATE_UNSIGNED_H(horiz_out0, FILTER_BITS, 7);

    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out0, (v16i8)horiz_out1);
    tmp6 = __msa_dotp_u_h(vec0, filt_vert);

    horiz_vec = (v16u8)__msa_vshf_b(mask, src3, src3);
    horiz_out1 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out1 = SRARI_SATURATE_UNSIGNED_H(horiz_out1, FILTER_BITS, 7);

    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out1, (v16i8)horiz_out0);
    tmp7 = __msa_dotp_u_h(vec0, filt_vert);

    horiz_vec = (v16u8)__msa_vshf_b(mask, src4, src4);
    horiz_out0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out0 = SRARI_SATURATE_UNSIGNED_H(horiz_out0, FILTER_BITS, 7);

    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out0, (v16i8)horiz_out1);
    tmp8 = __msa_dotp_u_h(vec0, filt_vert);

    tmp5 = SRARI_SATURATE_UNSIGNED_H(tmp5, FILTER_BITS, 7);
    tmp6 = SRARI_SATURATE_UNSIGNED_H(tmp6, FILTER_BITS, 7);
    tmp7 = SRARI_SATURATE_UNSIGNED_H(tmp7, FILTER_BITS, 7);
    tmp8 = SRARI_SATURATE_UNSIGNED_H(tmp8, FILTER_BITS, 7);

    PCKEV_B_STORE_8_BYTES_4(tmp5, tmp6, tmp7, tmp8, dst, dst_stride);
    dst += (4 * dst_stride);
  }
}

static void common_hv_2ht_2vt_8w_msa(const uint8_t *src, int32_t src_stride,
                                     uint8_t *dst, int32_t dst_stride,
                                     int8_t *filter_horiz, int8_t *filter_vert,
                                     int32_t height) {
  if (4 == height) {
    common_hv_2ht_2vt_8x4_msa(src, src_stride, dst, dst_stride, filter_horiz,
                              filter_vert);
  } else {
    common_hv_2ht_2vt_8x8mult_msa(src, src_stride, dst, dst_stride,
                                  filter_horiz, filter_vert, height);
  }
}

static void common_hv_2ht_2vt_16w_msa(const uint8_t *src, int32_t src_stride,
                                      uint8_t *dst, int32_t dst_stride,
                                      int8_t *filter_horiz, int8_t *filter_vert,
                                      int32_t height) {
  uint32_t loop_cnt;
  v16i8 src0, src1, src2, src3, src4, src5, src6, src7, mask;
  v16u8 filt_horiz, filt_vert, vec0, horiz_vec;
  v8u16 horiz_vec0, horiz_vec1, tmp1, tmp2;
  v8u16 horiz_out0, horiz_out1, horiz_out2, horiz_out3;
  v8i16 filt;

  mask = LOAD_SB(&mc_filt_mask_arr[0]);

  /* rearranging filter */
  filt = LOAD_SH(filter_horiz);
  filt_horiz = (v16u8)__msa_splati_h(filt, 0);

  filt = LOAD_SH(filter_vert);
  filt_vert = (v16u8)__msa_splati_h(filt, 0);

  src0 = LOAD_SB(src);
  src1 = LOAD_SB(src + 8);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src0, src0);
  horiz_vec0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out0 = SRARI_SATURATE_UNSIGNED_H(horiz_vec0, FILTER_BITS, 7);

  horiz_vec = (v16u8)__msa_vshf_b(mask, src1, src1);
  horiz_vec1 = __msa_dotp_u_h(horiz_vec, filt_horiz);
  horiz_out2 = SRARI_SATURATE_UNSIGNED_H(horiz_vec1, FILTER_BITS, 7);

  src += src_stride;

  for (loop_cnt = (height >> 2); loop_cnt--;) {
    LOAD_4VECS_SB(src, src_stride, src0, src2, src4, src6);
    LOAD_4VECS_SB(src + 8, src_stride, src1, src3, src5, src7);
    src += (4 * src_stride);

    horiz_vec = (v16u8)__msa_vshf_b(mask, src0, src0);
    horiz_vec0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out1 = SRARI_SATURATE_UNSIGNED_H(horiz_vec0, FILTER_BITS, 7);

    horiz_vec = (v16u8)__msa_vshf_b(mask, src1, src1);
    horiz_vec1 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out3 = SRARI_SATURATE_UNSIGNED_H(horiz_vec1, FILTER_BITS, 7);

    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out1, (v16i8)horiz_out0);
    tmp1 = __msa_dotp_u_h(vec0, filt_vert);
    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out3, (v16i8)horiz_out2);
    tmp2 = __msa_dotp_u_h(vec0, filt_vert);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);
    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp2, tmp1, dst);
    dst += dst_stride;

    horiz_vec = (v16u8)__msa_vshf_b(mask, src2, src2);
    horiz_vec0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out0 = SRARI_SATURATE_UNSIGNED_H(horiz_vec0, FILTER_BITS, 7);

    horiz_vec = (v16u8)__msa_vshf_b(mask, src3, src3);
    horiz_vec1 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out2 = SRARI_SATURATE_UNSIGNED_H(horiz_vec1, FILTER_BITS, 7);

    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out0, (v16i8)horiz_out1);
    tmp1 = __msa_dotp_u_h(vec0, filt_vert);
    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out2, (v16i8)horiz_out3);
    tmp2 = __msa_dotp_u_h(vec0, filt_vert);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);
    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp2, tmp1, dst);
    dst += dst_stride;

    horiz_vec = (v16u8)__msa_vshf_b(mask, src4, src4);
    horiz_vec0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out1 = SRARI_SATURATE_UNSIGNED_H(horiz_vec0, FILTER_BITS, 7);

    horiz_vec = (v16u8)__msa_vshf_b(mask, src5, src5);
    horiz_vec1 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out3 = SRARI_SATURATE_UNSIGNED_H(horiz_vec1, FILTER_BITS, 7);

    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out1, (v16i8)horiz_out0);
    tmp1 = __msa_dotp_u_h(vec0, filt_vert);
    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out3, (v16i8)horiz_out2);
    tmp2 = __msa_dotp_u_h(vec0, filt_vert);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);
    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp2, tmp1, dst);
    dst += dst_stride;

    horiz_vec = (v16u8)__msa_vshf_b(mask, src6, src6);
    horiz_vec0 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out0 = SRARI_SATURATE_UNSIGNED_H(horiz_vec0, FILTER_BITS, 7);

    horiz_vec = (v16u8)__msa_vshf_b(mask, src7, src7);
    horiz_vec1 = __msa_dotp_u_h(horiz_vec, filt_horiz);
    horiz_out2 = SRARI_SATURATE_UNSIGNED_H(horiz_vec1, FILTER_BITS, 7);

    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out0, (v16i8)horiz_out1);
    tmp1 = __msa_dotp_u_h(vec0, filt_vert);
    vec0 = (v16u8)__msa_ilvev_b((v16i8)horiz_out2, (v16i8)horiz_out3);
    tmp2 = __msa_dotp_u_h(vec0, filt_vert);
    tmp1 = SRARI_SATURATE_UNSIGNED_H(tmp1, FILTER_BITS, 7);
    tmp2 = SRARI_SATURATE_UNSIGNED_H(tmp2, FILTER_BITS, 7);

    PCKEV_B_STORE_VEC(tmp2, tmp1, dst);
    dst += dst_stride;
  }
}

static void common_hv_2ht_2vt_32w_msa(const uint8_t *src, int32_t src_stride,
                                      uint8_t *dst, int32_t dst_stride,
                                      int8_t *filter_horiz, int8_t *filter_vert,
                                      int32_t height) {
  int32_t multiple8_cnt;
  for (multiple8_cnt = 2; multiple8_cnt--;) {
    common_hv_2ht_2vt_16w_msa(src, src_stride, dst, dst_stride, filter_horiz,
                              filter_vert, height);
    src += 16;
    dst += 16;
  }
}

static void common_hv_2ht_2vt_64w_msa(const uint8_t *src, int32_t src_stride,
                                      uint8_t *dst, int32_t dst_stride,
                                      int8_t *filter_horiz, int8_t *filter_vert,
                                      int32_t height) {
  int32_t multiple8_cnt;
  for (multiple8_cnt = 4; multiple8_cnt--;) {
    common_hv_2ht_2vt_16w_msa(src, src_stride, dst, dst_stride, filter_horiz,
                              filter_vert, height);
    src += 16;
    dst += 16;
  }
}

void vp9_convolve8_msa(const uint8_t *src, ptrdiff_t src_stride,
                       uint8_t *dst, ptrdiff_t dst_stride,
                       const int16_t *filter_x, int32_t x_step_q4,
                       const int16_t *filter_y, int32_t y_step_q4,
                       int32_t w, int32_t h) {
  int8_t cnt, filt_hor[8], filt_ver[8];

  if (16 != x_step_q4 || 16 != y_step_q4) {
    vp9_convolve8_c(src, src_stride, dst, dst_stride,
                    filter_x, x_step_q4, filter_y, y_step_q4,
                    w, h);
    return;
  }

  if (((const int32_t *)filter_x)[1] == 0x800000 &&
      ((const int32_t *)filter_y)[1] == 0x800000) {
    vp9_convolve_copy(src, src_stride, dst, dst_stride,
                      filter_x, x_step_q4, filter_y, y_step_q4,
                      w, h);
    return;
  }

  for (cnt = 0; cnt < 8; ++cnt) {
    filt_hor[cnt] = filter_x[cnt];
    filt_ver[cnt] = filter_y[cnt];
  }

  if (((const int32_t *)filter_x)[0] == 0 &&
      ((const int32_t *)filter_y)[0] == 0) {
    switch (w) {
      case 4:
        common_hv_2ht_2vt_4w_msa(src, (int32_t)src_stride,
                                 dst, (int32_t)dst_stride,
                                 &filt_hor[3], &filt_ver[3], (int32_t)h);
        break;
      case 8:
        common_hv_2ht_2vt_8w_msa(src, (int32_t)src_stride,
                                 dst, (int32_t)dst_stride,
                                 &filt_hor[3], &filt_ver[3], (int32_t)h);
        break;
      case 16:
        common_hv_2ht_2vt_16w_msa(src, (int32_t)src_stride,
                                  dst, (int32_t)dst_stride,
                                  &filt_hor[3], &filt_ver[3], (int32_t)h);
        break;
      case 32:
        common_hv_2ht_2vt_32w_msa(src, (int32_t)src_stride,
                                  dst, (int32_t)dst_stride,
                                  &filt_hor[3], &filt_ver[3], (int32_t)h);
        break;
      case 64:
        common_hv_2ht_2vt_64w_msa(src, (int32_t)src_stride,
                                  dst, (int32_t)dst_stride,
                                  &filt_hor[3], &filt_ver[3], (int32_t)h);
        break;
      default:
        vp9_convolve8_c(src, src_stride, dst, dst_stride,
                        filter_x, x_step_q4, filter_y, y_step_q4,
                        w, h);
        break;
    }
  } else if (((const int32_t *)filter_x)[0] == 0 ||
             ((const int32_t *)filter_y)[0] == 0) {
    vp9_convolve8_c(src, src_stride, dst, dst_stride,
                    filter_x, x_step_q4, filter_y, y_step_q4,
                    w, h);
  } else {
    switch (w) {
      case 4:
        common_hv_8ht_8vt_4w_msa(src, (int32_t)src_stride,
                                 dst, (int32_t)dst_stride,
                                 filt_hor, filt_ver, (int32_t)h);
        break;
      case 8:
        common_hv_8ht_8vt_8w_msa(src, (int32_t)src_stride,
                                 dst, (int32_t)dst_stride,
                                 filt_hor, filt_ver, (int32_t)h);
        break;
      case 16:
        common_hv_8ht_8vt_16w_msa(src, (int32_t)src_stride,
                                  dst, (int32_t)dst_stride,
                                  filt_hor, filt_ver, (int32_t)h);
        break;
      case 32:
        common_hv_8ht_8vt_32w_msa(src, (int32_t)src_stride,
                                  dst, (int32_t)dst_stride,
                                  filt_hor, filt_ver, (int32_t)h);
        break;
      case 64:
        common_hv_8ht_8vt_64w_msa(src, (int32_t)src_stride,
                                  dst, (int32_t)dst_stride,
                                  filt_hor, filt_ver, (int32_t)h);
        break;
      default:
        vp9_convolve8_c(src, src_stride, dst, dst_stride,
                        filter_x, x_step_q4, filter_y, y_step_q4,
                        w, h);
        break;
    }
  }
}
