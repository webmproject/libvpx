/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <stdlib.h>

#include "./vpx_dsp_rtcd.h"

#include "vpx_dsp/ppc/types_vsx.h"

#include "vpx/vpx_integer.h"
#include "vpx_ports/mem.h"

#define PROCESS16(offset)           \
  v_a = vec_vsx_ld(offset, a);      \
  v_b = vec_vsx_ld(offset, b);      \
  v_ah = unpack_to_s16_h(v_a);      \
  v_al = unpack_to_s16_l(v_a);      \
  v_bh = unpack_to_s16_h(v_b);      \
  v_bl = unpack_to_s16_l(v_b);      \
  v_subh = vec_sub(v_ah, v_bh);     \
  v_subl = vec_sub(v_al, v_bl);     \
  v_absh = vec_abs(v_subh);         \
  v_absl = vec_abs(v_subl);         \
  v_sad = vec_sum4s(v_absh, v_sad); \
  v_sad = vec_sum4s(v_absl, v_sad);

#define SAD16(height)                                                     \
  unsigned int vpx_sad16x##height##_vsx(const uint8_t *a, int a_stride,   \
                                        const uint8_t *b, int b_stride) { \
    int y;                                                                \
    unsigned int sad[4];                                                  \
    uint8x16_t v_a, v_b;                                                  \
    int16x8_t v_ah, v_al, v_bh, v_bl, v_absh, v_absl, v_subh, v_subl;     \
    int32x4_t v_sad = vec_splat_s32(0);                                   \
                                                                          \
    for (y = 0; y < height; y++) {                                        \
      PROCESS16(0);                                                       \
                                                                          \
      a += a_stride;                                                      \
      b += b_stride;                                                      \
    }                                                                     \
    vec_vsx_st((uint32x4_t)v_sad, 0, sad);                                \
                                                                          \
    return sad[3] + sad[2] + sad[1] + sad[0];                             \
  }

#define SAD32(height)                                                     \
  unsigned int vpx_sad32x##height##_vsx(const uint8_t *a, int a_stride,   \
                                        const uint8_t *b, int b_stride) { \
    int y;                                                                \
    unsigned int sad[4];                                                  \
    uint8x16_t v_a, v_b;                                                  \
    int16x8_t v_ah, v_al, v_bh, v_bl, v_absh, v_absl, v_subh, v_subl;     \
    int32x4_t v_sad = vec_splat_s32(0);                                   \
                                                                          \
    for (y = 0; y < height; y++) {                                        \
      PROCESS16(0);                                                       \
      PROCESS16(16);                                                      \
                                                                          \
      a += a_stride;                                                      \
      b += b_stride;                                                      \
    }                                                                     \
    vec_vsx_st((uint32x4_t)v_sad, 0, sad);                                \
                                                                          \
    return sad[3] + sad[2] + sad[1] + sad[0];                             \
  }

#define SAD64(height)                                                     \
  unsigned int vpx_sad64x##height##_vsx(const uint8_t *a, int a_stride,   \
                                        const uint8_t *b, int b_stride) { \
    int y;                                                                \
    unsigned int sad[4];                                                  \
    uint8x16_t v_a, v_b;                                                  \
    int16x8_t v_ah, v_al, v_bh, v_bl, v_absh, v_absl, v_subh, v_subl;     \
    int32x4_t v_sad = vec_splat_s32(0);                                   \
                                                                          \
    for (y = 0; y < height; y++) {                                        \
      PROCESS16(0);                                                       \
      PROCESS16(16);                                                      \
      PROCESS16(32);                                                      \
      PROCESS16(48);                                                      \
                                                                          \
      a += a_stride;                                                      \
      b += b_stride;                                                      \
    }                                                                     \
    vec_vsx_st((uint32x4_t)v_sad, 0, sad);                                \
                                                                          \
    return sad[3] + sad[2] + sad[1] + sad[0];                             \
  }

SAD16(8);
SAD16(16);
SAD16(32);
SAD32(16);
SAD32(32);
SAD32(64);
SAD64(32);
SAD64(64);

#define SAD16AVG(height)                                                      \
  unsigned int vpx_sad16x##height##_avg_vsx(                                  \
      const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, \
      const uint8_t *second_pred) {                                           \
    DECLARE_ALIGNED(16, uint8_t, comp_pred[16 * height]);                     \
    vpx_comp_avg_pred_vsx(comp_pred, second_pred, 16, height, ref,            \
                          ref_stride);                                        \
                                                                              \
    return vpx_sad16x##height##_vsx(src, src_stride, comp_pred, 16);          \
  }

#define SAD32AVG(height)                                                      \
  unsigned int vpx_sad32x##height##_avg_vsx(                                  \
      const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, \
      const uint8_t *second_pred) {                                           \
    DECLARE_ALIGNED(32, uint8_t, comp_pred[32 * height]);                     \
    vpx_comp_avg_pred_vsx(comp_pred, second_pred, 32, height, ref,            \
                          ref_stride);                                        \
                                                                              \
    return vpx_sad32x##height##_vsx(src, src_stride, comp_pred, 32);          \
  }

#define SAD64AVG(height)                                                      \
  unsigned int vpx_sad64x##height##_avg_vsx(                                  \
      const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, \
      const uint8_t *second_pred) {                                           \
    DECLARE_ALIGNED(64, uint8_t, comp_pred[64 * height]);                     \
    vpx_comp_avg_pred_vsx(comp_pred, second_pred, 64, height, ref,            \
                          ref_stride);                                        \
                                                                              \
    return vpx_sad64x##height##_vsx(src, src_stride, comp_pred, 64);          \
  }

SAD16AVG(8);
SAD16AVG(16);
SAD16AVG(32);
SAD32AVG(16);
SAD32AVG(32);
SAD32AVG(64);
SAD64AVG(32);
SAD64AVG(64);
