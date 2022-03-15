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

#define SAD_UB2_UH(in0, in1, ref0, ref1)           \
  ({                                               \
    __m128i diff0_m, diff1_m, sad_m0;              \
    __m128i sad_m = __lsx_vldi(0);                 \
                                                   \
    diff0_m = __lsx_vabsd_bu(in0, ref0);           \
    diff1_m = __lsx_vabsd_bu(in1, ref1);           \
                                                   \
    sad_m0 = __lsx_vhaddw_hu_bu(diff0_m, diff0_m); \
    sad_m = __lsx_vadd_h(sad_m, sad_m0);           \
    sad_m0 = __lsx_vhaddw_hu_bu(diff1_m, diff1_m); \
    sad_m = __lsx_vadd_h(sad_m, sad_m0);           \
                                                   \
    sad_m;                                         \
  })

#define HADD_UW_U32(in)                          \
  ({                                             \
    __m128i res0_m;                              \
    uint32_t sum_m;                              \
    res0_m = __lsx_vhaddw_du_wu(in, in);         \
    res0_m = __lsx_vhaddw_qu_du(res0_m, res0_m); \
    sum_m = __lsx_vpickve2gr_w(res0_m, 0);       \
    sum_m;                                       \
  })

#define HADD_UH_U32(in)                 \
  ({                                    \
    __m128i res_m;                      \
    uint32_t sum_m;                     \
    res_m = __lsx_vhaddw_wu_hu(in, in); \
    sum_m = HADD_UW_U32(res_m);         \
    sum_m;                              \
  })

static uint32_t sad_32width_lsx(const uint8_t *src, int32_t src_stride,
                                const uint8_t *ref, int32_t ref_stride,
                                int32_t height) {
  int32_t ht_cnt = (height >> 2);
  __m128i src0, src1, ref0, ref1;
  __m128i sad_tmp;
  __m128i sad = __lsx_vldi(0);

  for (; ht_cnt--;) {
    DUP2_ARG2(__lsx_vld, src, 0, src, 16, src0, src1);
    src += src_stride;
    DUP2_ARG2(__lsx_vld, ref, 0, ref, 16, ref0, ref1);
    ref += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad = __lsx_vadd_h(sad, sad_tmp);

    DUP2_ARG2(__lsx_vld, src, 0, src, 16, src0, src1);
    src += src_stride;
    DUP2_ARG2(__lsx_vld, ref, 0, ref, 16, ref0, ref1);
    ref += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad = __lsx_vadd_h(sad, sad_tmp);

    DUP2_ARG2(__lsx_vld, src, 0, src, 16, src0, src1);
    src += src_stride;
    DUP2_ARG2(__lsx_vld, ref, 0, ref, 16, ref0, ref1);
    ref += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad = __lsx_vadd_h(sad, sad_tmp);

    DUP2_ARG2(__lsx_vld, src, 0, src, 16, src0, src1);
    src += src_stride;
    DUP2_ARG2(__lsx_vld, ref, 0, ref, 16, ref0, ref1);
    ref += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad = __lsx_vadd_h(sad, sad_tmp);
  }
  return HADD_UH_U32(sad);
}

static uint32_t sad_64width_lsx(const uint8_t *src, int32_t src_stride,
                                const uint8_t *ref, int32_t ref_stride,
                                int32_t height) {
  int32_t ht_cnt = (height >> 1);
  uint32_t sad = 0;
  __m128i src0, src1, src2, src3;
  __m128i ref0, ref1, ref2, ref3;
  __m128i sad_tmp;
  __m128i sad0 = __lsx_vldi(0);
  __m128i sad1 = sad0;

  for (; ht_cnt--;) {
    DUP4_ARG2(__lsx_vld, src, 0, src, 16, src, 32, src, 48, src0, src1, src2,
              src3);
    src += src_stride;
    DUP4_ARG2(__lsx_vld, ref, 0, ref, 16, ref, 32, ref, 48, ref0, ref1, ref2,
              ref3);
    ref += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad0 = __lsx_vadd_h(sad0, sad_tmp);
    sad_tmp = SAD_UB2_UH(src2, src3, ref2, ref3);
    sad1 = __lsx_vadd_h(sad1, sad_tmp);

    DUP4_ARG2(__lsx_vld, src, 0, src, 16, src, 32, src, 48, src0, src1, src2,
              src3);
    src += src_stride;
    DUP4_ARG2(__lsx_vld, ref, 0, ref, 16, ref, 32, ref, 48, ref0, ref1, ref2,
              ref3);
    ref += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad0 = __lsx_vadd_h(sad0, sad_tmp);
    sad_tmp = SAD_UB2_UH(src2, src3, ref2, ref3);
    sad1 = __lsx_vadd_h(sad1, sad_tmp);
  }

  sad = HADD_UH_U32(sad0);
  sad += HADD_UH_U32(sad1);

  return sad;
}

static void sad_16width_x4d_lsx(const uint8_t *src_ptr, int32_t src_stride,
                                const uint8_t *const aref_ptr[],
                                int32_t ref_stride, int32_t height,
                                uint32_t *sad_array) {
  int32_t ht_cnt = (height >> 1);
  const uint8_t *ref0_ptr, *ref1_ptr, *ref2_ptr, *ref3_ptr;
  __m128i src, ref0, ref1, ref2, ref3, diff, sad_tmp;
  __m128i sad0 = __lsx_vldi(0);
  __m128i sad1 = sad0;
  __m128i sad2 = sad0;
  __m128i sad3 = sad0;

  ref0_ptr = aref_ptr[0];
  ref1_ptr = aref_ptr[1];
  ref2_ptr = aref_ptr[2];
  ref3_ptr = aref_ptr[3];

  for (; ht_cnt--;) {
    src = __lsx_vld(src_ptr, 0);
    src_ptr += src_stride;
    ref0 = __lsx_vld(ref0_ptr, 0);
    ref0_ptr += ref_stride;
    ref1 = __lsx_vld(ref1_ptr, 0);
    ref1_ptr += ref_stride;
    ref2 = __lsx_vld(ref2_ptr, 0);
    ref2_ptr += ref_stride;
    ref3 = __lsx_vld(ref3_ptr, 0);
    ref3_ptr += ref_stride;

    diff = __lsx_vabsd_bu(src, ref0);
    sad_tmp = __lsx_vhaddw_hu_bu(diff, diff);
    sad0 = __lsx_vadd_h(sad0, sad_tmp);
    diff = __lsx_vabsd_bu(src, ref1);
    sad_tmp = __lsx_vhaddw_hu_bu(diff, diff);
    sad1 = __lsx_vadd_h(sad1, sad_tmp);
    diff = __lsx_vabsd_bu(src, ref2);
    sad_tmp = __lsx_vhaddw_hu_bu(diff, diff);
    sad2 = __lsx_vadd_h(sad2, sad_tmp);
    diff = __lsx_vabsd_bu(src, ref3);
    sad_tmp = __lsx_vhaddw_hu_bu(diff, diff);
    sad3 = __lsx_vadd_h(sad3, sad_tmp);

    src = __lsx_vld(src_ptr, 0);
    src_ptr += src_stride;
    ref0 = __lsx_vld(ref0_ptr, 0);
    ref0_ptr += ref_stride;
    ref1 = __lsx_vld(ref1_ptr, 0);
    ref1_ptr += ref_stride;
    ref2 = __lsx_vld(ref2_ptr, 0);
    ref2_ptr += ref_stride;
    ref3 = __lsx_vld(ref3_ptr, 0);
    ref3_ptr += ref_stride;

    diff = __lsx_vabsd_bu(src, ref0);
    sad_tmp = __lsx_vhaddw_hu_bu(diff, diff);
    sad0 = __lsx_vadd_h(sad0, sad_tmp);
    diff = __lsx_vabsd_bu(src, ref1);
    sad_tmp = __lsx_vhaddw_hu_bu(diff, diff);
    sad1 = __lsx_vadd_h(sad1, sad_tmp);
    diff = __lsx_vabsd_bu(src, ref2);
    sad_tmp = __lsx_vhaddw_hu_bu(diff, diff);
    sad2 = __lsx_vadd_h(sad2, sad_tmp);
    diff = __lsx_vabsd_bu(src, ref3);
    sad_tmp = __lsx_vhaddw_hu_bu(diff, diff);
    sad3 = __lsx_vadd_h(sad3, sad_tmp);
  }
  sad_array[0] = HADD_UH_U32(sad0);
  sad_array[1] = HADD_UH_U32(sad1);
  sad_array[2] = HADD_UH_U32(sad2);
  sad_array[3] = HADD_UH_U32(sad3);
}

static void sad_32width_x4d_lsx(const uint8_t *src, int32_t src_stride,
                                const uint8_t *const aref_ptr[],
                                int32_t ref_stride, int32_t height,
                                uint32_t *sad_array) {
  const uint8_t *ref0_ptr, *ref1_ptr, *ref2_ptr, *ref3_ptr;
  int32_t ht_cnt = height;
  __m128i src0, src1, ref0, ref1, sad_tmp;
  __m128i sad0 = __lsx_vldi(0);
  __m128i sad1 = sad0;
  __m128i sad2 = sad0;
  __m128i sad3 = sad0;

  ref0_ptr = aref_ptr[0];
  ref1_ptr = aref_ptr[1];
  ref2_ptr = aref_ptr[2];
  ref3_ptr = aref_ptr[3];

  for (; ht_cnt--;) {
    DUP2_ARG2(__lsx_vld, src, 0, src, 16, src0, src1);
    src += src_stride;

    DUP2_ARG2(__lsx_vld, ref0_ptr, 0, ref0_ptr, 16, ref0, ref1);
    ref0_ptr += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad0 = __lsx_vadd_h(sad0, sad_tmp);

    DUP2_ARG2(__lsx_vld, ref1_ptr, 0, ref1_ptr, 16, ref0, ref1);
    ref1_ptr += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad1 = __lsx_vadd_h(sad1, sad_tmp);

    DUP2_ARG2(__lsx_vld, ref2_ptr, 0, ref2_ptr, 16, ref0, ref1);
    ref2_ptr += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad2 = __lsx_vadd_h(sad2, sad_tmp);

    DUP2_ARG2(__lsx_vld, ref3_ptr, 0, ref3_ptr, 16, ref0, ref1);
    ref3_ptr += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad3 = __lsx_vadd_h(sad3, sad_tmp);
  }
  sad_array[0] = HADD_UH_U32(sad0);
  sad_array[1] = HADD_UH_U32(sad1);
  sad_array[2] = HADD_UH_U32(sad2);
  sad_array[3] = HADD_UH_U32(sad3);
}

static void sad_64width_x4d_lsx(const uint8_t *src, int32_t src_stride,
                                const uint8_t *const aref_ptr[],
                                int32_t ref_stride, int32_t height,
                                uint32_t *sad_array) {
  const uint8_t *ref0_ptr, *ref1_ptr, *ref2_ptr, *ref3_ptr;
  int32_t ht_cnt = height;
  __m128i src0, src1, src2, src3;
  __m128i ref0, ref1, ref2, ref3;
  __m128i sad, sad_tmp;

  __m128i sad0_0 = __lsx_vldi(0);
  __m128i sad0_1 = sad0_0;
  __m128i sad1_0 = sad0_0;
  __m128i sad1_1 = sad0_0;
  __m128i sad2_0 = sad0_0;
  __m128i sad2_1 = sad0_0;
  __m128i sad3_0 = sad0_0;
  __m128i sad3_1 = sad0_0;

  ref0_ptr = aref_ptr[0];
  ref1_ptr = aref_ptr[1];
  ref2_ptr = aref_ptr[2];
  ref3_ptr = aref_ptr[3];

  for (; ht_cnt--;) {
    DUP4_ARG2(__lsx_vld, src, 0, src, 16, src, 32, src, 48, src0, src1, src2,
              src3);
    src += src_stride;

    DUP4_ARG2(__lsx_vld, ref0_ptr, 0, ref0_ptr, 16, ref0_ptr, 32, ref0_ptr, 48,
              ref0, ref1, ref2, ref3);
    ref0_ptr += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad0_0 = __lsx_vadd_h(sad0_0, sad_tmp);
    sad_tmp = SAD_UB2_UH(src2, src3, ref2, ref3);
    sad0_1 = __lsx_vadd_h(sad0_1, sad_tmp);

    DUP4_ARG2(__lsx_vld, ref1_ptr, 0, ref1_ptr, 16, ref1_ptr, 32, ref1_ptr, 48,
              ref0, ref1, ref2, ref3);
    ref1_ptr += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad1_0 = __lsx_vadd_h(sad1_0, sad_tmp);
    sad_tmp = SAD_UB2_UH(src2, src3, ref2, ref3);
    sad1_1 = __lsx_vadd_h(sad1_1, sad_tmp);

    DUP4_ARG2(__lsx_vld, ref2_ptr, 0, ref2_ptr, 16, ref2_ptr, 32, ref2_ptr, 48,
              ref0, ref1, ref2, ref3);
    ref2_ptr += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad2_0 = __lsx_vadd_h(sad2_0, sad_tmp);
    sad_tmp = SAD_UB2_UH(src2, src3, ref2, ref3);
    sad2_1 = __lsx_vadd_h(sad2_1, sad_tmp);

    DUP4_ARG2(__lsx_vld, ref3_ptr, 0, ref3_ptr, 16, ref3_ptr, 32, ref3_ptr, 48,
              ref0, ref1, ref2, ref3);
    ref3_ptr += ref_stride;
    sad_tmp = SAD_UB2_UH(src0, src1, ref0, ref1);
    sad3_0 = __lsx_vadd_h(sad3_0, sad_tmp);
    sad_tmp = SAD_UB2_UH(src2, src3, ref2, ref3);
    sad3_1 = __lsx_vadd_h(sad3_1, sad_tmp);
  }
  sad = __lsx_vhaddw_wu_hu(sad0_0, sad0_0);
  sad_tmp = __lsx_vhaddw_wu_hu(sad0_1, sad0_1);
  sad = __lsx_vadd_w(sad, sad_tmp);
  sad_array[0] = HADD_UW_U32(sad);

  sad = __lsx_vhaddw_wu_hu(sad1_0, sad1_0);
  sad_tmp = __lsx_vhaddw_wu_hu(sad1_1, sad1_1);
  sad = __lsx_vadd_w(sad, sad_tmp);
  sad_array[1] = HADD_UW_U32(sad);

  sad = __lsx_vhaddw_wu_hu(sad2_0, sad2_0);
  sad_tmp = __lsx_vhaddw_wu_hu(sad2_1, sad2_1);
  sad = __lsx_vadd_w(sad, sad_tmp);
  sad_array[2] = HADD_UW_U32(sad);

  sad = __lsx_vhaddw_wu_hu(sad3_0, sad3_0);
  sad_tmp = __lsx_vhaddw_wu_hu(sad3_1, sad3_1);
  sad = __lsx_vadd_w(sad, sad_tmp);
  sad_array[3] = HADD_UW_U32(sad);
}

#define VPX_SAD_32xHT_LSX(height)                                             \
  uint32_t vpx_sad32x##height##_lsx(const uint8_t *src, int32_t src_stride,   \
                                    const uint8_t *ref, int32_t ref_stride) { \
    return sad_32width_lsx(src, src_stride, ref, ref_stride, height);         \
  }

#define VPX_SAD_64xHT_LSX(height)                                             \
  uint32_t vpx_sad64x##height##_lsx(const uint8_t *src, int32_t src_stride,   \
                                    const uint8_t *ref, int32_t ref_stride) { \
    return sad_64width_lsx(src, src_stride, ref, ref_stride, height);         \
  }

#define VPX_SAD_16xHTx4D_LSX(height)                                       \
  void vpx_sad16x##height##x4d_lsx(const uint8_t *src, int32_t src_stride, \
                                   const uint8_t *const refs[],            \
                                   int32_t ref_stride, uint32_t *sads) {   \
    sad_16width_x4d_lsx(src, src_stride, refs, ref_stride, height, sads);  \
  }

#define VPX_SAD_32xHTx4D_LSX(height)                                       \
  void vpx_sad32x##height##x4d_lsx(const uint8_t *src, int32_t src_stride, \
                                   const uint8_t *const refs[],            \
                                   int32_t ref_stride, uint32_t *sads) {   \
    sad_32width_x4d_lsx(src, src_stride, refs, ref_stride, height, sads);  \
  }

#define VPX_SAD_64xHTx4D_LSX(height)                                       \
  void vpx_sad64x##height##x4d_lsx(const uint8_t *src, int32_t src_stride, \
                                   const uint8_t *const refs[],            \
                                   int32_t ref_stride, uint32_t *sads) {   \
    sad_64width_x4d_lsx(src, src_stride, refs, ref_stride, height, sads);  \
  }

#define SAD64 VPX_SAD_64xHT_LSX(64) VPX_SAD_64xHTx4D_LSX(64)

SAD64

#define SAD32 VPX_SAD_32xHT_LSX(32) VPX_SAD_32xHTx4D_LSX(32)

SAD32

#define SAD16 VPX_SAD_16xHTx4D_LSX(16)

SAD16

#undef SAD64
#undef SAD32
#undef SAD16
