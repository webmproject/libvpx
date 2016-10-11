/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "./aom_dsp_rtcd.h"

SIMD_INLINE void calc_delta(v128 o, v128 x, v128 a, v128 b, v128 c, v128 d,
                            v128 e, v128 f, uint8_t *dst, v128 sp, v128 sm,
                            int dstride) {
  const v128 c8 = v128_dup_8(8);
  const v128 tmp =
      v128_add_8(v128_max_s8(v128_min_s8(v128_ssub_s8(c, x), sp), sm),
                 v128_max_s8(v128_min_s8(v128_ssub_s8(d, x), sp), sm));
  const v128 delta = v128_add_8(
      v128_add_8(
          v128_shl_8(
              v128_add_8(v128_max_s8(v128_min_s8(v128_ssub_s8(a, x), sp), sm),
                         v128_max_s8(v128_min_s8(v128_ssub_s8(f, x), sp), sm)),
              2),
          v128_add_8(v128_max_s8(v128_min_s8(v128_ssub_s8(b, x), sp), sm),
                     v128_max_s8(v128_min_s8(v128_ssub_s8(e, x), sp), sm))),
      v128_add_8(v128_add_8(tmp, tmp), tmp));
  o = v128_add_8(
      o,
      v128_shr_s8(
          v128_add_8(c8, v128_add_8(delta, v128_cmplt_s8(delta, v128_zero()))),
          4));
  v64_store_aligned(dst, v128_high_v64(o));
  v64_store_aligned(dst + dstride, v128_low_v64(o));
}

static void clpf_block(const uint8_t *src, uint8_t *dst, int sstride,
                       int dstride, int x0, int y0, int sizey, int width,
                       int height, unsigned int strength) {
  int bottom = height - 2 - y0;
  const v128 sp = v128_dup_8(strength);
  const v128 sm = v128_dup_8(-(int)strength);
  const v128 c128 = v128_dup_8(128);
  dst += x0 + y0 * dstride;
  src += x0 + y0 * sstride;

  if (!x0) {  // Clip left
    const v128 b_shuff = v128_from_v64(v64_from_64(0x0d0c0b0a09080808LL),
                                       v64_from_64(0x0504030201000000LL));
    const v128 c_shuff = v128_from_v64(v64_from_64(0x0e0d0c0b0a090808LL),
                                       v64_from_64(0x0605040302010000LL));
    int y;

    for (y = 0; y < sizey; y += 2) {
      const v64 l1 = v64_load_aligned(src);
      const v64 l2 = v64_load_aligned(src + sstride);
      v128 o = v128_from_v64(l1, l2);
      const v128 x = v128_add_8(c128, o);
      const v128 a = v128_add_8(
          c128,
          v128_from_v64(v64_load_aligned(src - (y != -y0) * sstride), l1));
      const v128 b = v128_shuffle_8(x, b_shuff);
      const v128 c = v128_shuffle_8(x, c_shuff);
      const v128 d = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src + 1),
                              v64_load_unaligned(src + 1 + sstride)));
      const v128 e = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src + 2),
                              v64_load_unaligned(src + 2 + sstride)));
      const v128 f = v128_add_8(
          c128, v128_from_v64(
                    l2, v64_load_aligned(src + ((y != bottom) + 1) * sstride)));
      calc_delta(o, x, a, b, c, d, e, f, dst, sp, sm, dstride);
      src += sstride * 2;
      dst += dstride * 2;
    }
  } else if (!(width - x0 - 8)) {  // Clip right
    const v128 d_shuff = v128_from_v64(v64_from_64(0x0f0f0e0d0c0b0a09LL),
                                       v64_from_64(0x0707060504030201LL));
    const v128 e_shuff = v128_from_v64(v64_from_64(0x0f0f0f0e0d0c0b0aLL),
                                       v64_from_64(0x0707070605040302LL));
    int y;

    for (y = 0; y < sizey; y += 2) {
      const v64 l1 = v64_load_aligned(src);
      const v64 l2 = v64_load_aligned(src + sstride);
      v128 o = v128_from_v64(l1, l2);
      const v128 x = v128_add_8(c128, o);
      const v128 a = v128_add_8(
          c128,
          v128_from_v64(v64_load_aligned(src - (y != -y0) * sstride), l1));
      const v128 b = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src - 2),
                              v64_load_unaligned(src - 2 + sstride)));
      const v128 c = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src - 1),
                              v64_load_unaligned(src - 1 + sstride)));
      const v128 d = v128_shuffle_8(x, d_shuff);
      const v128 e = v128_shuffle_8(x, e_shuff);
      const v128 f = v128_add_8(
          c128, v128_from_v64(
                    l2, v64_load_aligned(src + ((y != bottom) + 1) * sstride)));
      calc_delta(o, x, a, b, c, d, e, f, dst, sp, sm, dstride);
      src += sstride * 2;
      dst += dstride * 2;
    }
  } else {  // No left/right clipping
    int y;
    for (y = 0; y < sizey; y += 2) {
      const v64 l1 = v64_load_aligned(src);
      const v64 l2 = v64_load_aligned(src + sstride);
      v128 o = v128_from_v64(l1, l2);
      const v128 x = v128_add_8(c128, o);
      const v128 a = v128_add_8(
          c128,
          v128_from_v64(v64_load_aligned(src - (y != -y0) * sstride), l1));
      const v128 b = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src - 2),
                              v64_load_unaligned(src - 2 + sstride)));
      const v128 c = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src - 1),
                              v64_load_unaligned(src - 1 + sstride)));
      const v128 d = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src + 1),
                              v64_load_unaligned(src + 1 + sstride)));
      const v128 e = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src + 2),
                              v64_load_unaligned(src + 2 + sstride)));
      const v128 f = v128_add_8(
          c128, v128_from_v64(
                    l2, v64_load_aligned(src + ((y != bottom) + 1) * sstride)));
      calc_delta(o, x, a, b, c, d, e, f, dst, sp, sm, dstride);
      src += sstride * 2;
      dst += dstride * 2;
    }
  }
}

void SIMD_FUNC(aom_clpf_block)(const uint8_t *src, uint8_t *dst, int sstride,
                               int dstride, int x0, int y0, int sizex,
                               int sizey, int width, int height,
                               unsigned int strength) {
  // TODO(stemidts):
  // A sizex different from 8 will only be needed if CLPF is extended to chroma.
  // This will only be used if 4:2:0 and width not a multiple of 16 and along
  // the right edge only, so we can fall back to the plain C implementation in
  // this case.  If not extended to chroma, this test will be redundant.
  if (sizex != 8 || width < 16 || y0 + 8 > height || x0 + 8 > width) {
    // Fallback to C for odd sizes
    aom_clpf_block_c(src, dst, sstride, dstride, x0, y0, sizex, sizey, width,
                     height, strength);
  } else {
    clpf_block(src, dst, sstride, dstride, x0, y0, sizey, width, height,
               strength);
  }
}

#if CONFIG_AOM_HIGHBITDEPTH
static void calc_delta_hbd(v128 o, v128 a, v128 b, v128 c, v128 d, v128 e,
                           v128 f, uint16_t *dst, v128 sp, v128 sm) {
  const v128 c8 = v128_dup_16(8);
  const v128 tmp =
      v128_add_16(v128_max_s16(v128_min_s16(v128_sub_16(c, o), sp), sm),
                  v128_max_s16(v128_min_s16(v128_sub_16(d, o), sp), sm));
  const v128 delta = v128_add_16(
      v128_add_16(
          v128_shl_16(
              v128_add_16(
                  v128_max_s16(v128_min_s16(v128_sub_16(a, o), sp), sm),
                  v128_max_s16(v128_min_s16(v128_sub_16(f, o), sp), sm)),
              2),
          v128_add_16(v128_max_s16(v128_min_s16(v128_sub_16(b, o), sp), sm),
                      v128_max_s16(v128_min_s16(v128_sub_16(e, o), sp), sm))),
      v128_add_16(v128_add_16(tmp, tmp), tmp));
  v128_store_aligned(
      dst,
      v128_add_16(
          o, v128_shr_s16(
                 v128_add_16(c8, v128_add_16(delta, v128_cmplt_s16(
                                                        delta, v128_zero()))),
                 4)));
}

SIMD_INLINE void clpf_block_hbd(const uint16_t *src, uint16_t *dst, int sstride,
                                int dstride, int x0, int y0, int sizey,
                                int width, int height, unsigned int strength) {
  int y;
  int bottom = height - 2 - y0;
  const v128 sp = v128_dup_16(strength);
  const v128 sm = v128_dup_16(-(int)strength);

  dst += x0 + y0 * dstride;
  src += x0 + y0 * sstride;

  if (!x0) {  // Clip left
    const v128 b_shuff = v128_from_v64(v64_from_64(0x0b0a090807060504LL),
                                       v64_from_64(0x0302010001000100LL));
    const v128 c_shuff = v128_from_v64(v64_from_64(0x0d0c0b0a09080706LL),
                                       v64_from_64(0x0504030201000100LL));
    for (y = 0; y < sizey; y++) {
      const v128 o = v128_load_aligned(src);
      const v128 a = v128_load_aligned(src - (y != -y0) * sstride);
      const v128 b = v128_shuffle_8(o, b_shuff);
      const v128 c = v128_shuffle_8(o, c_shuff);
      const v128 d = v128_load_unaligned(src + 1);
      const v128 e = v128_load_unaligned(src + 2);
      const v128 f = v128_load_aligned(src + (y - 1 != bottom) * sstride);
      calc_delta_hbd(o, a, b, c, d, e, f, dst, sp, sm);
      src += sstride;
      dst += dstride;
    }
  } else if (!(width - x0 - 8)) {  // Clip right
    const v128 d_shuff = v128_from_v64(v64_from_64(0x0f0e0f0e0d0c0b0aLL),
                                       v64_from_64(0x0908070605040302LL));
    const v128 e_shuff = v128_from_v64(v64_from_64(0x0f0e0f0e0f0e0d0cLL),
                                       v64_from_64(0x0b0a090807060504LL));
    for (y = 0; y < sizey; y++) {
      const v128 o = v128_load_aligned(src);
      const v128 a = v128_load_aligned(src - (y != -y0) * sstride);
      const v128 b = v128_load_unaligned(src - 2);
      const v128 c = v128_load_unaligned(src - 1);
      const v128 d = v128_shuffle_8(o, d_shuff);
      const v128 e = v128_shuffle_8(o, e_shuff);
      const v128 f = v128_load_aligned(src + (y - 1 != bottom) * sstride);
      calc_delta_hbd(o, a, b, c, d, e, f, dst, sp, sm);
      src += sstride;
      dst += dstride;
    }
  } else {  // No left/right clipping
    for (y = 0; y < sizey; y++) {
      const v128 o = v128_load_aligned(src);
      const v128 a = v128_load_aligned(src - (y != -y0) * sstride);
      const v128 b = v128_load_unaligned(src - 2);
      const v128 c = v128_load_unaligned(src - 1);
      const v128 d = v128_load_unaligned(src + 1);
      const v128 e = v128_load_unaligned(src + 2);
      const v128 f = v128_load_aligned(src + (y - 1 != bottom) * sstride);
      calc_delta_hbd(o, a, b, c, d, e, f, dst, sp, sm);
      src += sstride;
      dst += dstride;
    }
  }
}

void SIMD_FUNC(aom_clpf_block_hbd)(const uint16_t *src, uint16_t *dst,
                                   int sstride, int dstride, int x0, int y0,
                                   int sizex, int sizey, int width, int height,
                                   unsigned int strength) {
  if (sizex != 8 || width < 16 || y0 + 8 > height || x0 + 8 > width) {
    // Fallback to C for odd sizes
    aom_clpf_block_hbd_c(src, dst, sstride, dstride, x0, y0, sizex, sizey,
                         width, height, strength);
  } else {
    clpf_block_hbd(src, dst, sstride, dstride, x0, y0, sizey, width, height,
                   strength);
  }
}
#endif
