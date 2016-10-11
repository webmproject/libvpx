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
#include "aom_dsp/aom_simd.h"
#include "aom_ports/mem.h"

SIMD_INLINE void calc_diff(v128 o, v128 *a, v128 *b, v128 *c, v128 *d, v128 *e,
                           v128 *f) {
  // The difference will be 9 bit, offset by 128 so we can use saturated
  // sub to avoid going to 16 bit temporarily before "strength" clipping.
  const v128 c128 = v128_dup_8(128);
  v128 x = v128_add_8(c128, o);
  *a = v128_ssub_s8(v128_add_8(c128, *a), x);
  *b = v128_ssub_s8(v128_add_8(c128, *b), x);
  *c = v128_ssub_s8(v128_add_8(c128, *c), x);
  *d = v128_ssub_s8(v128_add_8(c128, *d), x);
  *e = v128_ssub_s8(v128_add_8(c128, *e), x);
  *f = v128_ssub_s8(v128_add_8(c128, *f), x);
}

SIMD_INLINE v128 delta_kernel(v128 o, v128 a, v128 b, v128 c, v128 d, v128 e,
                              v128 f, v128 sp, v128 sm) {
  const v128 tmp = v128_add_8(v128_max_s8(v128_min_s8(c, sp), sm),
                              v128_max_s8(v128_min_s8(d, sp), sm));
  const v128 delta = v128_add_8(
      v128_add_8(v128_shl_8(v128_add_8(v128_max_s8(v128_min_s8(a, sp), sm),
                                       v128_max_s8(v128_min_s8(f, sp), sm)),
                            2),
                 v128_add_8(v128_max_s8(v128_min_s8(b, sp), sm),
                            v128_max_s8(v128_min_s8(e, sp), sm))),
      v128_add_8(v128_add_8(tmp, tmp), tmp));

  return v128_add_8(
      o, v128_shr_s8(
             v128_add_8(v128_dup_8(8),
                        v128_add_8(delta, v128_cmplt_s8(delta, v128_zero()))),
             4));
}

SIMD_INLINE v128 calc_delta(v128 o, v128 a, v128 b, v128 c, v128 d, v128 e,
                            v128 f, v128 sp, v128 sm) {
  calc_diff(o, &a, &b, &c, &d, &e, &f);
  return delta_kernel(o, a, b, c, d, e, f, sp, sm);
}

SIMD_INLINE void clip_sides(v128 *b, v128 *c, v128 *d, v128 *e, int left,
                            int right) {
  DECLARE_ALIGNED(16, static const uint64_t,
                  b_shuff[]) = { 0x0504030201000000LL, 0x0d0c0b0a09080808LL };
  DECLARE_ALIGNED(16, static const uint64_t,
                  c_shuff[]) = { 0x0605040302010000LL, 0x0e0d0c0b0a090808LL };
  DECLARE_ALIGNED(16, static const uint64_t,
                  d_shuff[]) = { 0x0707060504030201LL, 0x0f0f0e0d0c0b0a09LL };
  DECLARE_ALIGNED(16, static const uint64_t,
                  e_shuff[]) = { 0x0707070605040302LL, 0x0f0f0f0e0d0c0b0aLL };

  if (!left) {  // Left clipping
    *b = v128_shuffle_8(*b, v128_load_aligned(b_shuff));
    *c = v128_shuffle_8(*c, v128_load_aligned(c_shuff));
  }
  if (!right) {  // Right clipping
    *d = v128_shuffle_8(*d, v128_load_aligned(d_shuff));
    *e = v128_shuffle_8(*e, v128_load_aligned(e_shuff));
  }
}

SIMD_INLINE void read_two_lines(const uint8_t *rec, const uint8_t *org,
                                int rstride, int ostride, int x0, int y0,
                                int bottom, int right, int y, v128 *o, v128 *r,
                                v128 *a, v128 *b, v128 *c, v128 *d, v128 *e,
                                v128 *f) {
  const v64 k1 = v64_load_aligned(org);
  const v64 k2 = v64_load_aligned(org + ostride);
  const v64 l1 = v64_load_aligned(rec);
  const v64 l2 = v64_load_aligned(rec + rstride);
  *o = v128_from_v64(k1, k2);
  *r = v128_from_v64(l1, l2);
  *a = v128_from_v64(v64_load_aligned(rec - (y != -y0) * rstride), l1);
  *f = v128_from_v64(l2, v64_load_aligned(rec + ((y != bottom) + 1) * rstride));
  *b = v128_from_v64(v64_load_unaligned(rec - 2 * !!x0),
                     v64_load_unaligned(rec - 2 * !!x0 + rstride));
  *c = v128_from_v64(v64_load_unaligned(rec - !!x0),
                     v64_load_unaligned(rec - !!x0 + rstride));
  *d = v128_from_v64(v64_load_unaligned(rec + !!right),
                     v64_load_unaligned(rec + !!right + rstride));
  *e = v128_from_v64(v64_load_unaligned(rec + 2 * !!right),
                     v64_load_unaligned(rec + 2 * !!right + rstride));
  clip_sides(b, c, d, e, x0, right);
}

void SIMD_FUNC(aom_clpf_detect)(const uint8_t *rec, const uint8_t *org,
                                int rstride, int ostride, int x0, int y0,
                                int width, int height, int *sum0, int *sum1,
                                unsigned int strength, int size) {
  const v128 sp = v128_dup_8(strength);
  const v128 sm = v128_dup_8(-(int)strength);
  const int right = width - 8 - x0;
  const int bottom = height - 2 - y0;
  ssd128_internal ssd0 = v128_ssd_u8_init();
  ssd128_internal ssd1 = v128_ssd_u8_init();
  int y;

  if (size != 8) {  // Fallback to plain C
    aom_clpf_detect_c(rec, org, rstride, ostride, x0, y0, width, height, sum0,
                      sum1, strength, size);
    return;
  }

  rec += x0 + y0 * rstride;
  org += x0 + y0 * ostride;

  for (y = 0; y < 8; y += 2) {
    v128 a, b, c, d, e, f, o, r;
    read_two_lines(rec, org, rstride, ostride, x0, y0, bottom, right, y, &o, &r,
                   &a, &b, &c, &d, &e, &f);
    ssd0 = v128_ssd_u8(ssd0, o, r);
    ssd1 = v128_ssd_u8(ssd1, o, calc_delta(r, a, b, c, d, e, f, sp, sm));
    rec += rstride * 2;
    org += ostride * 2;
  }
  *sum0 += v128_ssd_u8_sum(ssd0);
  *sum1 += v128_ssd_u8_sum(ssd1);
}

SIMD_INLINE void calc_delta_multi(v128 r, v128 o, v128 a, v128 b, v128 c,
                                  v128 d, v128 e, v128 f, ssd128_internal *ssd1,
                                  ssd128_internal *ssd2,
                                  ssd128_internal *ssd3) {
  calc_diff(r, &a, &b, &c, &d, &e, &f);
  *ssd1 = v128_ssd_u8(*ssd1, o, delta_kernel(r, a, b, c, d, e, f, v128_dup_8(1),
                                             v128_dup_8(-1)));
  *ssd2 = v128_ssd_u8(*ssd2, o, delta_kernel(r, a, b, c, d, e, f, v128_dup_8(2),
                                             v128_dup_8(-2)));
  *ssd3 = v128_ssd_u8(*ssd3, o, delta_kernel(r, a, b, c, d, e, f, v128_dup_8(4),
                                             v128_dup_8(-4)));
}

// Test multiple filter strengths at once.
void SIMD_FUNC(aom_clpf_detect_multi)(const uint8_t *rec, const uint8_t *org,
                                      int rstride, int ostride, int x0, int y0,
                                      int width, int height, int *sum,
                                      int size) {
  const int bottom = height - 2 - y0;
  const int right = width - 8 - x0;
  ssd128_internal ssd0 = v128_ssd_u8_init();
  ssd128_internal ssd1 = v128_ssd_u8_init();
  ssd128_internal ssd2 = v128_ssd_u8_init();
  ssd128_internal ssd3 = v128_ssd_u8_init();
  int y;

  if (size != 8) {  // Fallback to plain C
    aom_clpf_detect_multi_c(rec, org, rstride, ostride, x0, y0, width, height,
                            sum, size);
    return;
  }

  rec += x0 + y0 * rstride;
  org += x0 + y0 * ostride;

  for (y = 0; y < 8; y += 2) {
    v128 a, b, c, d, e, f, o, r;
    read_two_lines(rec, org, rstride, ostride, x0, y0, bottom, right, y, &o, &r,
                   &a, &b, &c, &d, &e, &f);
    ssd0 = v128_ssd_u8(ssd0, o, r);
    calc_delta_multi(r, o, a, b, c, d, e, f, &ssd1, &ssd2, &ssd3);
    rec += 2 * rstride;
    org += 2 * ostride;
  }
  sum[0] += v128_ssd_u8_sum(ssd0);
  sum[1] += v128_ssd_u8_sum(ssd1);
  sum[2] += v128_ssd_u8_sum(ssd2);
  sum[3] += v128_ssd_u8_sum(ssd3);
}

#if CONFIG_AOM_HIGHBITDEPTH
SIMD_INLINE void read_two_lines_hbd(const uint16_t *rec, const uint16_t *org,
                                    int rstride, int ostride, int x0, int y0,
                                    int bottom, int right, int y, v128 *o,
                                    v128 *r, v128 *a, v128 *b, v128 *c, v128 *d,
                                    v128 *e, v128 *f, int shift) {
  const v128 n1 = v128_shr_u16(v128_load_aligned(rec), shift);
  const v128 n2 = v128_shr_u16(v128_load_aligned(rec + rstride), shift);
  *o = v128_unziplo_8(v128_shr_u16(v128_load_aligned(org), shift),
                      v128_shr_u16(v128_load_aligned(org + ostride), shift));
  *r = v128_unziplo_8(n1, n2);
  *a = v128_unziplo_8(
      v128_shr_u16(v128_load_aligned(rec - (y != -y0) * rstride), shift), n1);
  *f = v128_unziplo_8(
      n2, v128_shr_u16(v128_load_unaligned(rec + ((y != bottom) + 1) * rstride),
                       shift));
  *b = v128_unziplo_8(
      v128_shr_u16(v128_load_unaligned(rec - 2 * !!x0), shift),
      v128_shr_u16(v128_load_unaligned(rec - 2 * !!x0 + rstride), shift));
  *c = v128_unziplo_8(
      v128_shr_u16(v128_load_unaligned(rec - !!x0), shift),
      v128_shr_u16(v128_load_unaligned(rec - !!x0 + rstride), shift));
  *d = v128_unziplo_8(
      v128_shr_u16(v128_load_unaligned(rec + !!right), shift),
      v128_shr_u16(v128_load_unaligned(rec + !!right + rstride), shift));
  *e = v128_unziplo_8(
      v128_shr_u16(v128_load_unaligned(rec + 2 * !!right), shift),
      v128_shr_u16(v128_load_unaligned(rec + 2 * !!right + rstride), shift));
  clip_sides(b, c, d, e, x0, right);
}

void SIMD_FUNC(aom_clpf_detect_hbd)(const uint16_t *rec, const uint16_t *org,
                                    int rstride, int ostride, int x0, int y0,
                                    int width, int height, int *sum0, int *sum1,
                                    unsigned int strength, int shift,
                                    int size) {
  const v128 sp = v128_dup_8(strength >> shift);
  const v128 sm = v128_dup_8(-(int)(strength >> shift));
  const int bottom = height - 2 - y0;
  const int right = width - 8 - x0;
  ssd128_internal ssd0 = v128_ssd_u8_init();
  ssd128_internal ssd1 = v128_ssd_u8_init();
  int y;

  if (size != 8) {  // Fallback to plain C
    aom_clpf_detect_hbd_c(rec, org, rstride, ostride, x0, y0, width, height,
                          sum0, sum1, strength, shift, size);
    return;
  }

  rec += x0 + y0 * rstride;
  org += x0 + y0 * ostride;

  for (y = 0; y < 8; y += 2) {
    v128 a, b, c, d, e, f, o, r;
    read_two_lines_hbd(rec, org, rstride, ostride, x0, y0, bottom, right, y, &o,
                       &r, &a, &b, &c, &d, &e, &f, shift);
    ssd0 = v128_ssd_u8(ssd0, o, r);
    ssd1 = v128_ssd_u8(ssd1, o, calc_delta(r, a, b, c, d, e, f, sp, sm));
    rec += rstride * 2;
    org += ostride * 2;
  }
  *sum0 += v128_ssd_u8_sum(ssd0);
  *sum1 += v128_ssd_u8_sum(ssd1);
}

void SIMD_FUNC(aom_clpf_detect_multi_hbd)(const uint16_t *rec,
                                          const uint16_t *org, int rstride,
                                          int ostride, int x0, int y0,
                                          int width, int height, int *sum,
                                          int shift, int size) {
  const int bottom = height - 2 - y0;
  const int right = width - 8 - x0;
  ssd128_internal ssd0 = v128_ssd_u8_init();
  ssd128_internal ssd1 = v128_ssd_u8_init();
  ssd128_internal ssd2 = v128_ssd_u8_init();
  ssd128_internal ssd3 = v128_ssd_u8_init();
  int y;

  if (size != 8) {  // Fallback to plain C
    aom_clpf_detect_multi_hbd_c(rec, org, rstride, ostride, x0, y0, width,
                                height, sum, shift, size);
    return;
  }

  rec += x0 + y0 * rstride;
  org += x0 + y0 * ostride;

  for (y = 0; y < 8; y += 2) {
    v128 a, b, c, d, e, f, o, r;
    read_two_lines_hbd(rec, org, rstride, ostride, x0, y0, bottom, right, y, &o,
                       &r, &a, &b, &c, &d, &e, &f, shift);
    ssd0 = v128_ssd_u8(ssd0, o, r);
    calc_delta_multi(r, o, a, b, c, d, e, f, &ssd1, &ssd2, &ssd3);
    rec += 2 * rstride;
    org += 2 * ostride;
  }
  sum[0] += v128_ssd_u8_sum(ssd0);
  sum[1] += v128_ssd_u8_sum(ssd1);
  sum[2] += v128_ssd_u8_sum(ssd2);
  sum[3] += v128_ssd_u8_sum(ssd3);
}
#endif
