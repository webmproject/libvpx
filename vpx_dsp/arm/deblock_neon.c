/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <arm_neon.h>
#include <assert.h>

#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/arm/transpose_neon.h"

static uint8x8_t average_k_out(const uint8x8_t a2, const uint8x8_t a1,
                               const uint8x8_t v0, const uint8x8_t b1,
                               const uint8x8_t b2) {
  const uint8x8_t k1 = vrhadd_u8(a2, a1);
  const uint8x8_t k2 = vrhadd_u8(b2, b1);
  const uint8x8_t k3 = vrhadd_u8(k1, k2);
  return vrhadd_u8(k3, v0);
}

static uint8x8_t generate_mask(const uint8x8_t a2, const uint8x8_t a1,
                               const uint8x8_t v0, const uint8x8_t b1,
                               const uint8x8_t b2, const uint8x8_t filter) {
  const uint8x8_t a2_v0 = vabd_u8(a2, v0);
  const uint8x8_t a1_v0 = vabd_u8(a1, v0);
  const uint8x8_t b1_v0 = vabd_u8(b1, v0);
  const uint8x8_t b2_v0 = vabd_u8(b2, v0);

  uint8x8_t max = vmax_u8(a2_v0, a1_v0);
  max = vmax_u8(b1_v0, max);
  max = vmax_u8(b2_v0, max);
  return vclt_u8(max, filter);
}

static uint8x8_t generate_output(const uint8x8_t a2, const uint8x8_t a1,
                                 const uint8x8_t v0, const uint8x8_t b1,
                                 const uint8x8_t b2, const uint8x8_t filter) {
  const uint8x8_t k_out = average_k_out(a2, a1, v0, b1, b2);
  const uint8x8_t mask = generate_mask(a2, a1, v0, b1, b2, filter);

  return vbsl_u8(mask, k_out, v0);
}

// Same functions but for uint8x16_t.
static uint8x16_t average_k_outq(const uint8x16_t a2, const uint8x16_t a1,
                                 const uint8x16_t v0, const uint8x16_t b1,
                                 const uint8x16_t b2) {
  const uint8x16_t k1 = vrhaddq_u8(a2, a1);
  const uint8x16_t k2 = vrhaddq_u8(b2, b1);
  const uint8x16_t k3 = vrhaddq_u8(k1, k2);
  return vrhaddq_u8(k3, v0);
}

static uint8x16_t generate_maskq(const uint8x16_t a2, const uint8x16_t a1,
                                 const uint8x16_t v0, const uint8x16_t b1,
                                 const uint8x16_t b2, const uint8x16_t filter) {
  const uint8x16_t a2_v0 = vabdq_u8(a2, v0);
  const uint8x16_t a1_v0 = vabdq_u8(a1, v0);
  const uint8x16_t b1_v0 = vabdq_u8(b1, v0);
  const uint8x16_t b2_v0 = vabdq_u8(b2, v0);

  uint8x16_t max = vmaxq_u8(a2_v0, a1_v0);
  max = vmaxq_u8(b1_v0, max);
  max = vmaxq_u8(b2_v0, max);
  return vcltq_u8(max, filter);
}

static uint8x16_t generate_outputq(const uint8x16_t a2, const uint8x16_t a1,
                                   const uint8x16_t v0, const uint8x16_t b1,
                                   const uint8x16_t b2,
                                   const uint8x16_t filter) {
  const uint8x16_t k_out = average_k_outq(a2, a1, v0, b1, b2);
  const uint8x16_t mask = generate_maskq(a2, a1, v0, b1, b2, filter);

  return vbslq_u8(mask, k_out, v0);
}

void vpx_post_proc_down_and_across_mb_row_neon(uint8_t *src_ptr,
                                               uint8_t *dst_ptr, int src_stride,
                                               int dst_stride, int cols,
                                               uint8_t *f, int size) {
  uint8_t *src, *dst;
  int row;
  int col;

  // Process a stripe of macroblocks. The stripe will be a multiple of 16 (for
  // Y) or 8 (for U/V) wide (cols) and the height (size) will be 16 (for Y) or 8
  // (for U/V).
  assert((size == 8 || size == 16) && cols % 8 == 0);

  // While columns of length 16 can be processed, load them.
  for (col = 0; col < cols - 8; col += 16) {
    uint8x16_t a0, a1, a2, a3, a4, a5, a6, a7;
    src = src_ptr - 2 * src_stride;
    dst = dst_ptr;

    a0 = vld1q_u8(src);
    src += src_stride;
    a1 = vld1q_u8(src);
    src += src_stride;
    a2 = vld1q_u8(src);
    src += src_stride;
    a3 = vld1q_u8(src);
    src += src_stride;

    for (row = 0; row < size; row += 4) {
      uint8x16_t v_out_0, v_out_1, v_out_2, v_out_3;
      const uint8x16_t filterq = vld1q_u8(f + col);

      a4 = vld1q_u8(src);
      src += src_stride;
      a5 = vld1q_u8(src);
      src += src_stride;
      a6 = vld1q_u8(src);
      src += src_stride;
      a7 = vld1q_u8(src);
      src += src_stride;

      v_out_0 = generate_outputq(a0, a1, a2, a3, a4, filterq);
      v_out_1 = generate_outputq(a1, a2, a3, a4, a5, filterq);
      v_out_2 = generate_outputq(a2, a3, a4, a5, a6, filterq);
      v_out_3 = generate_outputq(a3, a4, a5, a6, a7, filterq);

      vst1q_u8(dst, v_out_0);
      dst += dst_stride;
      vst1q_u8(dst, v_out_1);
      dst += dst_stride;
      vst1q_u8(dst, v_out_2);
      dst += dst_stride;
      vst1q_u8(dst, v_out_3);
      dst += dst_stride;

      // Rotate over to the next slot.
      a0 = a4;
      a1 = a5;
      a2 = a6;
      a3 = a7;
    }

    src_ptr += 16;
    dst_ptr += 16;
  }

  // Clean up any left over column of length 8.
  if (col != cols) {
    uint8x8_t a0, a1, a2, a3, a4, a5, a6, a7;
    src = src_ptr - 2 * src_stride;
    dst = dst_ptr;

    a0 = vld1_u8(src);
    src += src_stride;
    a1 = vld1_u8(src);
    src += src_stride;
    a2 = vld1_u8(src);
    src += src_stride;
    a3 = vld1_u8(src);
    src += src_stride;

    for (row = 0; row < size; row += 4) {
      uint8x8_t v_out_0, v_out_1, v_out_2, v_out_3;
      const uint8x8_t filter = vld1_u8(f + col);

      a4 = vld1_u8(src);
      src += src_stride;
      a5 = vld1_u8(src);
      src += src_stride;
      a6 = vld1_u8(src);
      src += src_stride;
      a7 = vld1_u8(src);
      src += src_stride;

      v_out_0 = generate_output(a0, a1, a2, a3, a4, filter);
      v_out_1 = generate_output(a1, a2, a3, a4, a5, filter);
      v_out_2 = generate_output(a2, a3, a4, a5, a6, filter);
      v_out_3 = generate_output(a3, a4, a5, a6, a7, filter);

      vst1_u8(dst, v_out_0);
      dst += dst_stride;
      vst1_u8(dst, v_out_1);
      dst += dst_stride;
      vst1_u8(dst, v_out_2);
      dst += dst_stride;
      vst1_u8(dst, v_out_3);
      dst += dst_stride;

      // Rotate over to the next slot.
      a0 = a4;
      a1 = a5;
      a2 = a6;
      a3 = a7;
    }

    // Not strictly necessary but makes resetting dst_ptr easier.
    dst_ptr += 8;
  }

  dst_ptr -= cols;

  for (row = 0; row < size; row += 8) {
    uint8x8_t a0, a1, a2, a3;
    uint8x8_t b0, b1, b2, b3, b4, b5, b6, b7;

    src = dst_ptr;
    dst = dst_ptr;

    // Load 8 values, transpose 4 of them, and discard 2 because they will be
    // reloaded later.
    load_and_transpose_u8_4x8(src, dst_stride, &a0, &a1, &a2, &a3);
    a3 = a1;
    a2 = a1 = a0;  // Extend left border.

    src += 2;

    for (col = 0; col < cols; col += 8) {
      uint8x8_t v_out_0, v_out_1, v_out_2, v_out_3, v_out_4, v_out_5, v_out_6,
          v_out_7;
      // Although the filter is meant to be applied vertically and is instead
      // being applied horizontally here it's OK because it's set in blocks of 8
      // (or 16).
      const uint8x8_t filter = vld1_u8(f + col);

      load_and_transpose_u8_8x8(src, dst_stride, &b0, &b1, &b2, &b3, &b4, &b5,
                                &b6, &b7);

      if (col + 8 == cols) {
        // Last row. Extend border (b5).
        b6 = b7 = b5;
      }

      v_out_0 = generate_output(a0, a1, a2, a3, b0, filter);
      v_out_1 = generate_output(a1, a2, a3, b0, b1, filter);
      v_out_2 = generate_output(a2, a3, b0, b1, b2, filter);
      v_out_3 = generate_output(a3, b0, b1, b2, b3, filter);
      v_out_4 = generate_output(b0, b1, b2, b3, b4, filter);
      v_out_5 = generate_output(b1, b2, b3, b4, b5, filter);
      v_out_6 = generate_output(b2, b3, b4, b5, b6, filter);
      v_out_7 = generate_output(b3, b4, b5, b6, b7, filter);

      transpose_and_store_u8_8x8(dst, dst_stride, v_out_0, v_out_1, v_out_2,
                                 v_out_3, v_out_4, v_out_5, v_out_6, v_out_7);

      a0 = b4;
      a1 = b5;
      a2 = b6;
      a3 = b7;

      src += 8;
      dst += 8;
    }

    dst_ptr += 8 * dst_stride;
  }
}
