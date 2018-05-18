/*
 *  Copyright (c) 2018 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/ppc/types_vsx.h"

extern const int16_t vpx_rv[];

static const uint8x16_t load_merge = { 0x00, 0x02, 0x04, 0x06, 0x08, 0x0A,
                                       0x0C, 0x0E, 0x18, 0x19, 0x1A, 0x1B,
                                       0x1C, 0x1D, 0x1E, 0x1F };

static const uint8x16_t mask_merge = { 0x00, 0x01, 0x10, 0x11, 0x04, 0x05,
                                       0x14, 0x15, 0x08, 0x09, 0x18, 0x19,
                                       0x0C, 0x0D, 0x1C, 0x1D };

void vpx_mbpost_proc_down_vsx(uint8_t *dst, int pitch, int rows, int cols,
                              int flimit) {
  int col, row, i;
  int16x8_t window[16];
  const int32x4_t lim = vec_splats(flimit);

  // 8 columns are processed at a time.
  assert(cols % 8 == 0);
  // If rows is less than 8 the bottom border extension fails.
  assert(rows >= 8);

  for (col = 0; col < cols; col += 8) {
    // The sum is signed and requires at most 13 bits.
    // (8 bits + sign) * 15 (4 bits)
    int16x8_t r1, sum;
    // The sum of squares requires at most 20 bits.
    // (16 bits + sign) * 15 (4 bits)
    int32x4_t sumsq_even, sumsq_odd;

    r1 = unpack_to_s16_h(vec_vsx_ld(0, dst));
    // Fill sliding window with first row.
    for (i = 0; i <= 8; i++) {
      window[i] = r1;
    }
    // First 9 rows of the sliding window are the same.
    // sum = r1 * 9
    sum = vec_mladd(r1, vec_splats((int16_t)9), vec_zeros_s16);

    // sumsq = r1 * r1 * 9
    sumsq_even = vec_mule(sum, r1);
    sumsq_odd = vec_mulo(sum, r1);

    // Fill the next 6 rows of the sliding window with rows 2 to 7.
    for (i = 1; i <= 6; ++i) {
      const int16x8_t next_row = unpack_to_s16_h(vec_vsx_ld(i * pitch, dst));
      window[i + 8] = next_row;
      sum = vec_add(sum, next_row);
      sumsq_odd = vec_add(sumsq_odd, vec_mulo(next_row, next_row));
      sumsq_even = vec_add(sumsq_even, vec_mule(next_row, next_row));
    }

    for (row = 0; row < rows; row++) {
      int32x4_t d15_even, d15_odd, d0_even, d0_odd;
      int32x4_t sumsq_odd_scaled, sumsq_even_scaled;
      int32x4_t thres_odd, thres_even;
      bool32x4_t mask_odd, mask_even;
      bool16x8_t mask;
      int16x8_t filtered, masked;
      uint8x16_t out;

      const int16x8_t rv = vec_vsx_ld(0, vpx_rv + (row & 127));

      // Move the sliding window
      if (row + 7 < rows) {
        window[15] = unpack_to_s16_h(vec_vsx_ld((row + 7) * pitch, dst));
      } else {
        window[15] = window[14];
      }

      // C: sum += s[7 * pitch] - s[-8 * pitch];
      sum = vec_add(sum, vec_sub(window[15], window[0]));

      // C: sumsq += s[7 * pitch] * s[7 * pitch] - s[-8 * pitch] * s[-8 *
      // pitch];
      // Optimization Note: Caching a squared-window for odd and even is
      // slower than just repeating the multiplies.
      d15_odd = vec_mulo(window[15], window[15]);
      d15_even = vec_mule(window[15], window[15]);
      d0_odd = vec_mulo(window[0], window[0]);
      d0_even = vec_mule(window[0], window[0]);
      sumsq_odd = vec_add(sumsq_odd, vec_sub(d15_odd, d0_odd));
      sumsq_even = vec_add(sumsq_even, vec_sub(d15_even, d0_even));

      // C: sumsq * 15 - sum * sum
      sumsq_odd_scaled = vec_mul(sumsq_odd, vec_splats((int32_t)15));
      sumsq_even_scaled = vec_mul(sumsq_even, vec_splats((int32_t)15));
      thres_odd = vec_sub(sumsq_odd_scaled, vec_mulo(sum, sum));
      thres_even = vec_sub(sumsq_even_scaled, vec_mule(sum, sum));

      // C: (vpx_rv[(r & 127) + (c & 7)] + sum + s[0]) >> 4
      filtered = vec_add(vec_add(rv, sum), window[8]);
      filtered = vec_sra(filtered, vec_splats((uint16_t)4));

      mask_odd = vec_cmplt(thres_odd, lim);
      mask_even = vec_cmplt(thres_even, lim);
      mask = vec_perm((bool16x8_t)mask_even, (bool16x8_t)mask_odd, mask_merge);
      masked = vec_sel(window[8], filtered, mask);

      // TODO(ltrudeau) If cols % 16 == 0, we could just process 16 per
      // iteration
      out = vec_perm((uint8x16_t)masked, vec_vsx_ld(0, dst + row * pitch),
                     load_merge);
      vec_vsx_st(out, 0, dst + row * pitch);

      // Optimization Note: Turns out that the following loop is faster than
      // using pointers to manage the sliding window.
      for (i = 1; i < 16; i++) {
        window[i - 1] = window[i];
      }
    }
    dst += 8;
  }
}
