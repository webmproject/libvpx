/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/common/mips/msa/vp9_macros_msa.h"

static void avg_width4_msa(const uint8_t *src, int32_t src_stride,
                           uint8_t *dst, int32_t dst_stride, int32_t height) {
  int32_t cnt;
  uint32_t out0, out1, out2, out3;
  v16u8 src0, src1, src2, src3;
  v16u8 dst0, dst1, dst2, dst3;

  if (0 == (height % 4)) {
    for (cnt = (height / 4); cnt--;) {
      LOAD_4VECS_UB(src, src_stride, src0, src1, src2, src3);
      src += (4 * src_stride);

      LOAD_4VECS_UB(dst, dst_stride, dst0, dst1, dst2, dst3);

      dst0 = __msa_aver_u_b(src0, dst0);
      dst1 = __msa_aver_u_b(src1, dst1);
      dst2 = __msa_aver_u_b(src2, dst2);
      dst3 = __msa_aver_u_b(src3, dst3);

      out0 = __msa_copy_u_w((v4i32)dst0, 0);
      out1 = __msa_copy_u_w((v4i32)dst1, 0);
      out2 = __msa_copy_u_w((v4i32)dst2, 0);
      out3 = __msa_copy_u_w((v4i32)dst3, 0);

      STORE_WORD(dst, out0);
      dst += dst_stride;
      STORE_WORD(dst, out1);
      dst += dst_stride;
      STORE_WORD(dst, out2);
      dst += dst_stride;
      STORE_WORD(dst, out3);
      dst += dst_stride;
    }
  } else if (0 == (height % 2)) {
    for (cnt = (height / 2); cnt--;) {
      LOAD_2VECS_UB(src, src_stride, src0, src1);
      src += (2 * src_stride);

      LOAD_2VECS_UB(dst, dst_stride, dst0, dst1);

      dst0 = __msa_aver_u_b(src0, dst0);
      dst1 = __msa_aver_u_b(src1, dst1);

      out0 = __msa_copy_u_w((v4i32)dst0, 0);
      out1 = __msa_copy_u_w((v4i32)dst1, 0);

      STORE_WORD(dst, out0);
      dst += dst_stride;
      STORE_WORD(dst, out1);
      dst += dst_stride;
    }
  }
}

static void avg_width8_msa(const uint8_t *src, int32_t src_stride,
                           uint8_t *dst, int32_t dst_stride, int32_t height) {
  int32_t cnt;
  uint64_t out0, out1, out2, out3;
  v16u8 src0, src1, src2, src3;
  v16u8 dst0, dst1, dst2, dst3;

  for (cnt = (height / 4); cnt--;) {
    LOAD_4VECS_UB(src, src_stride, src0, src1, src2, src3);
    src += (4 * src_stride);

    LOAD_4VECS_UB(dst, dst_stride, dst0, dst1, dst2, dst3);

    dst0 = __msa_aver_u_b(src0, dst0);
    dst1 = __msa_aver_u_b(src1, dst1);
    dst2 = __msa_aver_u_b(src2, dst2);
    dst3 = __msa_aver_u_b(src3, dst3);

    out0 = __msa_copy_u_d((v2i64)dst0, 0);
    out1 = __msa_copy_u_d((v2i64)dst1, 0);
    out2 = __msa_copy_u_d((v2i64)dst2, 0);
    out3 = __msa_copy_u_d((v2i64)dst3, 0);

    STORE_DWORD(dst, out0);
    dst += dst_stride;
    STORE_DWORD(dst, out1);
    dst += dst_stride;
    STORE_DWORD(dst, out2);
    dst += dst_stride;
    STORE_DWORD(dst, out3);
    dst += dst_stride;
  }
}

static void avg_width16_msa(const uint8_t *src, int32_t src_stride,
                            uint8_t *dst, int32_t dst_stride, int32_t height) {
  int32_t cnt;
  v16u8 src0, src1, src2, src3, src4, src5, src6, src7;
  v16u8 dst0, dst1, dst2, dst3, dst4, dst5, dst6, dst7;

  for (cnt = (height / 8); cnt--;) {
    LOAD_8VECS_UB(src, src_stride,
                  src0, src1, src2, src3, src4, src5, src6, src7);
    src += (8 * src_stride);

    LOAD_8VECS_UB(dst, dst_stride,
                  dst0, dst1, dst2, dst3, dst4, dst5, dst6, dst7);

    dst0 = __msa_aver_u_b(src0, dst0);
    dst1 = __msa_aver_u_b(src1, dst1);
    dst2 = __msa_aver_u_b(src2, dst2);
    dst3 = __msa_aver_u_b(src3, dst3);
    dst4 = __msa_aver_u_b(src4, dst4);
    dst5 = __msa_aver_u_b(src5, dst5);
    dst6 = __msa_aver_u_b(src6, dst6);
    dst7 = __msa_aver_u_b(src7, dst7);

    STORE_8VECS_UB(dst, dst_stride,
                   dst0, dst1, dst2, dst3, dst4, dst5, dst6, dst7);
    dst += (8 * dst_stride);
  }
}

static void avg_width32_msa(const uint8_t *src, int32_t src_stride,
                            uint8_t *dst, int32_t dst_stride, int32_t height) {
  int32_t cnt;
  uint8_t *dst_dup = dst;
  v16u8 src0, src1, src2, src3, src4, src5, src6, src7;
  v16u8 src8, src9, src10, src11, src12, src13, src14, src15;
  v16u8 dst0, dst1, dst2, dst3, dst4, dst5, dst6, dst7;
  v16u8 dst8, dst9, dst10, dst11, dst12, dst13, dst14, dst15;

  for (cnt = (height / 8); cnt--;) {
    src0 = LOAD_UB(src);
    src1 = LOAD_UB(src + 16);
    src += src_stride;
    src2 = LOAD_UB(src);
    src3 = LOAD_UB(src + 16);
    src += src_stride;
    src4 = LOAD_UB(src);
    src5 = LOAD_UB(src + 16);
    src += src_stride;
    src6 = LOAD_UB(src);
    src7 = LOAD_UB(src + 16);
    src += src_stride;

    dst0 = LOAD_UB(dst_dup);
    dst1 = LOAD_UB(dst_dup + 16);
    dst_dup += dst_stride;
    dst2 = LOAD_UB(dst_dup);
    dst3 = LOAD_UB(dst_dup + 16);
    dst_dup += dst_stride;
    dst4 = LOAD_UB(dst_dup);
    dst5 = LOAD_UB(dst_dup + 16);
    dst_dup += dst_stride;
    dst6 = LOAD_UB(dst_dup);
    dst7 = LOAD_UB(dst_dup + 16);
    dst_dup += dst_stride;

    src8 = LOAD_UB(src);
    src9 = LOAD_UB(src + 16);
    src += src_stride;
    src10 = LOAD_UB(src);
    src11 = LOAD_UB(src + 16);
    src += src_stride;
    src12 = LOAD_UB(src);
    src13 = LOAD_UB(src + 16);
    src += src_stride;
    src14 = LOAD_UB(src);
    src15 = LOAD_UB(src + 16);
    src += src_stride;

    dst8 = LOAD_UB(dst_dup);
    dst9 = LOAD_UB(dst_dup + 16);
    dst_dup += dst_stride;
    dst10 = LOAD_UB(dst_dup);
    dst11 = LOAD_UB(dst_dup + 16);
    dst_dup += dst_stride;
    dst12 = LOAD_UB(dst_dup);
    dst13 = LOAD_UB(dst_dup + 16);
    dst_dup += dst_stride;
    dst14 = LOAD_UB(dst_dup);
    dst15 = LOAD_UB(dst_dup + 16);
    dst_dup += dst_stride;

    dst0 = __msa_aver_u_b(src0, dst0);
    dst1 = __msa_aver_u_b(src1, dst1);
    dst2 = __msa_aver_u_b(src2, dst2);
    dst3 = __msa_aver_u_b(src3, dst3);
    dst4 = __msa_aver_u_b(src4, dst4);
    dst5 = __msa_aver_u_b(src5, dst5);
    dst6 = __msa_aver_u_b(src6, dst6);
    dst7 = __msa_aver_u_b(src7, dst7);
    dst8 = __msa_aver_u_b(src8, dst8);
    dst9 = __msa_aver_u_b(src9, dst9);
    dst10 = __msa_aver_u_b(src10, dst10);
    dst11 = __msa_aver_u_b(src11, dst11);
    dst12 = __msa_aver_u_b(src12, dst12);
    dst13 = __msa_aver_u_b(src13, dst13);
    dst14 = __msa_aver_u_b(src14, dst14);
    dst15 = __msa_aver_u_b(src15, dst15);

    STORE_UB(dst0, dst);
    STORE_UB(dst1, dst + 16);
    dst += dst_stride;
    STORE_UB(dst2, dst);
    STORE_UB(dst3, dst + 16);
    dst += dst_stride;
    STORE_UB(dst4, dst);
    STORE_UB(dst5, dst + 16);
    dst += dst_stride;
    STORE_UB(dst6, dst);
    STORE_UB(dst7, dst + 16);
    dst += dst_stride;
    STORE_UB(dst8, dst);
    STORE_UB(dst9, dst + 16);
    dst += dst_stride;
    STORE_UB(dst10, dst);
    STORE_UB(dst11, dst + 16);
    dst += dst_stride;
    STORE_UB(dst12, dst);
    STORE_UB(dst13, dst + 16);
    dst += dst_stride;
    STORE_UB(dst14, dst);
    STORE_UB(dst15, dst + 16);
    dst += dst_stride;
  }
}

static void avg_width64_msa(const uint8_t *src, int32_t src_stride,
                            uint8_t *dst, int32_t dst_stride, int32_t height) {
  int32_t cnt;
  uint8_t *dst_dup = dst;
  v16u8 src0, src1, src2, src3, src4, src5, src6, src7;
  v16u8 src8, src9, src10, src11, src12, src13, src14, src15;
  v16u8 dst0, dst1, dst2, dst3, dst4, dst5, dst6, dst7;
  v16u8 dst8, dst9, dst10, dst11, dst12, dst13, dst14, dst15;

  for (cnt = (height / 4); cnt--;) {
    LOAD_4VECS_UB(src, 16, src0, src1, src2, src3);
    src += src_stride;
    LOAD_4VECS_UB(src, 16, src4, src5, src6, src7);
    src += src_stride;
    LOAD_4VECS_UB(src, 16, src8, src9, src10, src11);
    src += src_stride;
    LOAD_4VECS_UB(src, 16, src12, src13, src14, src15);
    src += src_stride;

    LOAD_4VECS_UB(dst_dup, 16, dst0, dst1, dst2, dst3);
    dst_dup += dst_stride;
    LOAD_4VECS_UB(dst_dup, 16, dst4, dst5, dst6, dst7);
    dst_dup += dst_stride;
    LOAD_4VECS_UB(dst_dup, 16, dst8, dst9, dst10, dst11);
    dst_dup += dst_stride;
    LOAD_4VECS_UB(dst_dup, 16, dst12, dst13, dst14, dst15);
    dst_dup += dst_stride;

    dst0 = __msa_aver_u_b(src0, dst0);
    dst1 = __msa_aver_u_b(src1, dst1);
    dst2 = __msa_aver_u_b(src2, dst2);
    dst3 = __msa_aver_u_b(src3, dst3);
    dst4 = __msa_aver_u_b(src4, dst4);
    dst5 = __msa_aver_u_b(src5, dst5);
    dst6 = __msa_aver_u_b(src6, dst6);
    dst7 = __msa_aver_u_b(src7, dst7);
    dst8 = __msa_aver_u_b(src8, dst8);
    dst9 = __msa_aver_u_b(src9, dst9);
    dst10 = __msa_aver_u_b(src10, dst10);
    dst11 = __msa_aver_u_b(src11, dst11);
    dst12 = __msa_aver_u_b(src12, dst12);
    dst13 = __msa_aver_u_b(src13, dst13);
    dst14 = __msa_aver_u_b(src14, dst14);
    dst15 = __msa_aver_u_b(src15, dst15);

    STORE_4VECS_UB(dst, 16, dst0, dst1, dst2, dst3);
    dst += dst_stride;
    STORE_4VECS_UB(dst, 16, dst4, dst5, dst6, dst7);
    dst += dst_stride;
    STORE_4VECS_UB(dst, 16, dst8, dst9, dst10, dst11);
    dst += dst_stride;
    STORE_4VECS_UB(dst, 16, dst12, dst13, dst14, dst15);
    dst += dst_stride;
  }
}

void vp9_convolve_avg_msa(const uint8_t *src, ptrdiff_t src_stride,
                          uint8_t *dst, ptrdiff_t dst_stride,
                          const int16_t *filter_x, int32_t filter_x_stride,
                          const int16_t *filter_y, int32_t filter_y_stride,
                          int32_t w, int32_t h) {
  (void)filter_x;
  (void)filter_y;
  (void)filter_x_stride;
  (void)filter_y_stride;

  switch (w) {
    case 4: {
      avg_width4_msa(src, src_stride, dst, dst_stride, h);
      break;
    }
    case 8: {
      avg_width8_msa(src, src_stride, dst, dst_stride, h);
      break;
    }
    case 16: {
      avg_width16_msa(src, src_stride, dst, dst_stride, h);
      break;
    }
    case 32: {
      avg_width32_msa(src, src_stride, dst, dst_stride, h);
      break;
    }
    case 64: {
      avg_width64_msa(src, src_stride, dst, dst_stride, h);
      break;
    }
    default: {
      int32_t lp, cnt;
      for (cnt = h; cnt--;) {
        for (lp = 0; lp < w; ++lp) {
          dst[lp] = (((dst[lp] + src[lp]) + 1) >> 1);
        }
        src += src_stride;
        dst += dst_stride;
      }
      break;
    }
  }
}
