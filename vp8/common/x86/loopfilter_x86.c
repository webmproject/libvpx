/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <emmintrin.h>  // SSE2
#include "vpx_config.h"
#include "vp8/common/loopfilter.h"

prototype_loopfilter(vp8_mbloop_filter_vertical_edge_mmx);
prototype_loopfilter(vp8_loop_filter_vertical_edge_mmx);
prototype_loopfilter(vp8_loop_filter_horizontal_edge_mmx);

prototype_loopfilter(vp8_loop_filter_vertical_edge_sse2);
prototype_loopfilter(vp8_loop_filter_horizontal_edge_sse2);
prototype_loopfilter(vp8_mbloop_filter_vertical_edge_sse2);

extern loop_filter_uvfunction vp8_loop_filter_horizontal_edge_uv_sse2;
extern loop_filter_uvfunction vp8_loop_filter_vertical_edge_uv_sse2;
extern loop_filter_uvfunction vp8_mbloop_filter_vertical_edge_uv_sse2;

#if HAVE_MMX
/* Horizontal MB filtering */
void vp8_loop_filter_mbh_mmx(unsigned char *y_ptr, unsigned char *u_ptr, unsigned char *v_ptr,
                             int y_stride, int uv_stride, struct loop_filter_info *lfi) {
}


/* Vertical MB Filtering */
void vp8_loop_filter_mbv_mmx(unsigned char *y_ptr, unsigned char *u_ptr, unsigned char *v_ptr,
                             int y_stride, int uv_stride, struct loop_filter_info *lfi) {
  vp8_mbloop_filter_vertical_edge_mmx(y_ptr, y_stride, lfi->mblim, lfi->lim, lfi->hev_thr, 2);

  if (u_ptr)
    vp8_mbloop_filter_vertical_edge_mmx(u_ptr, uv_stride, lfi->mblim, lfi->lim, lfi->hev_thr, 1);

  if (v_ptr)
    vp8_mbloop_filter_vertical_edge_mmx(v_ptr, uv_stride, lfi->mblim, lfi->lim, lfi->hev_thr, 1);
}


/* Horizontal B Filtering */
void vp8_loop_filter_bh_mmx(unsigned char *y_ptr, unsigned char *u_ptr, unsigned char *v_ptr,
                            int y_stride, int uv_stride, struct loop_filter_info *lfi) {

}


void vp8_loop_filter_bhs_mmx(unsigned char *y_ptr, int y_stride, const unsigned char *blimit) {
  vp8_loop_filter_simple_horizontal_edge_mmx(y_ptr + 4 * y_stride, y_stride, blimit);
  vp8_loop_filter_simple_horizontal_edge_mmx(y_ptr + 8 * y_stride, y_stride, blimit);
  vp8_loop_filter_simple_horizontal_edge_mmx(y_ptr + 12 * y_stride, y_stride, blimit);
}


/* Vertical B Filtering */
void vp8_loop_filter_bv_mmx(unsigned char *y_ptr, unsigned char *u_ptr, unsigned char *v_ptr,
                            int y_stride, int uv_stride, struct loop_filter_info *lfi) {
  vp8_loop_filter_vertical_edge_mmx(y_ptr + 4, y_stride, lfi->blim, lfi->lim, lfi->hev_thr, 2);
  vp8_loop_filter_vertical_edge_mmx(y_ptr + 8, y_stride, lfi->blim, lfi->lim, lfi->hev_thr, 2);
  vp8_loop_filter_vertical_edge_mmx(y_ptr + 12, y_stride, lfi->blim, lfi->lim, lfi->hev_thr, 2);

  if (u_ptr)
    vp8_loop_filter_vertical_edge_mmx(u_ptr + 4, uv_stride, lfi->blim, lfi->lim, lfi->hev_thr, 1);

  if (v_ptr)
    vp8_loop_filter_vertical_edge_mmx(v_ptr + 4, uv_stride, lfi->blim, lfi->lim, lfi->hev_thr, 1);
}


void vp8_loop_filter_bvs_mmx(unsigned char *y_ptr, int y_stride, const unsigned char *blimit) {
  vp8_loop_filter_simple_vertical_edge_mmx(y_ptr + 4, y_stride, blimit);
  vp8_loop_filter_simple_vertical_edge_mmx(y_ptr + 8, y_stride, blimit);
  vp8_loop_filter_simple_vertical_edge_mmx(y_ptr + 12, y_stride, blimit);
}
#endif


#if HAVE_SSE2
void vp8_mbloop_filter_horizontal_edge_c_sse2
(
  unsigned char *s,
  int p,
  const unsigned char *_blimit,
  const unsigned char *_limit,
  const unsigned char *_thresh,
  int count
) {
  DECLARE_ALIGNED(16, unsigned char, flat_op2[16]);
  DECLARE_ALIGNED(16, unsigned char, flat_op1[16]);
  DECLARE_ALIGNED(16, unsigned char, flat_op0[16]);
  DECLARE_ALIGNED(16, unsigned char, flat_oq2[16]);
  DECLARE_ALIGNED(16, unsigned char, flat_oq1[16]);
  DECLARE_ALIGNED(16, unsigned char, flat_oq0[16]);
  __m128i mask, hev, flat;
  __m128i thresh, limit, blimit;
  const __m128i zero = _mm_set1_epi16(0);
  __m128i p4, p3, p2, p1, p0, q0, q1, q2, q3, q4;

  thresh = _mm_shuffle_epi32(_mm_cvtsi32_si128(_thresh[0] * 0x01010101), 0);
  limit = _mm_shuffle_epi32(_mm_cvtsi32_si128(_limit[0] * 0x01010101), 0);
  blimit = _mm_shuffle_epi32(_mm_cvtsi32_si128(_blimit[0] * 0x01010101), 0);

  p4 = _mm_loadu_si128((__m128i *)(s - 5 * p));
  p3 = _mm_loadu_si128((__m128i *)(s - 4 * p));
  p2 = _mm_loadu_si128((__m128i *)(s - 3 * p));
  p1 = _mm_loadu_si128((__m128i *)(s - 2 * p));
  p0 = _mm_loadu_si128((__m128i *)(s - 1 * p));
  q0 = _mm_loadu_si128((__m128i *)(s - 0 * p));
  q1 = _mm_loadu_si128((__m128i *)(s + 1 * p));
  q2 = _mm_loadu_si128((__m128i *)(s + 2 * p));
  q3 = _mm_loadu_si128((__m128i *)(s + 3 * p));
  q4 = _mm_loadu_si128((__m128i *)(s + 4 * p));
  {
    const __m128i abs_p1p0 = _mm_or_si128(_mm_subs_epu8(p1, p0),
                                          _mm_subs_epu8(p0, p1));
    const __m128i abs_q1q0 = _mm_or_si128(_mm_subs_epu8(q1, q0),
                                          _mm_subs_epu8(q0, q1));
    const __m128i one = _mm_set1_epi8(1);
    const __m128i fe = _mm_set1_epi8(0xfe);
    const __m128i ff = _mm_cmpeq_epi8(abs_p1p0, abs_p1p0);
    __m128i abs_p0q0 = _mm_or_si128(_mm_subs_epu8(p0, q0),
                                    _mm_subs_epu8(q0, p0));
    __m128i abs_p1q1 = _mm_or_si128(_mm_subs_epu8(p1, q1),
                                    _mm_subs_epu8(q1, p1));
    __m128i work;
    flat = _mm_max_epu8(abs_p1p0, abs_q1q0);
    hev = _mm_subs_epu8(flat, thresh);
    hev = _mm_xor_si128(_mm_cmpeq_epi8(hev, zero), ff);

    abs_p0q0 =_mm_adds_epu8(abs_p0q0, abs_p0q0);
    abs_p1q1 = _mm_srli_epi16(_mm_and_si128(abs_p1q1, fe), 1);
    mask = _mm_subs_epu8(_mm_adds_epu8(abs_p0q0, abs_p1q1), blimit);
    mask = _mm_xor_si128(_mm_cmpeq_epi8(mask, zero), ff);
    // mask |= (abs(p0 - q0) * 2 + abs(p1 - q1) / 2  > blimit) * -1;
    mask = _mm_max_epu8(flat, mask);
    // mask |= (abs(p1 - p0) > limit) * -1;
    // mask |= (abs(q1 - q0) > limit) * -1;
    work = _mm_max_epu8(_mm_or_si128(_mm_subs_epu8(p2, p1),
                                     _mm_subs_epu8(p1, p2)),
                         _mm_or_si128(_mm_subs_epu8(p3, p2),
                                      _mm_subs_epu8(p2, p3)));
    mask = _mm_max_epu8(work, mask);
    work = _mm_max_epu8(_mm_or_si128(_mm_subs_epu8(q2, q1),
                                     _mm_subs_epu8(q1, q2)),
                         _mm_or_si128(_mm_subs_epu8(q3, q2),
                                      _mm_subs_epu8(q2, q3)));
    mask = _mm_max_epu8(work, mask);
    mask = _mm_subs_epu8(mask, limit);
    mask = _mm_cmpeq_epi8(mask, zero);

    work = _mm_max_epu8(_mm_or_si128(_mm_subs_epu8(p2, p0),
                                     _mm_subs_epu8(p0, p2)),
                         _mm_or_si128(_mm_subs_epu8(q2, q0),
                                      _mm_subs_epu8(q0, q2)));
    flat = _mm_max_epu8(work, flat);
    work = _mm_max_epu8(_mm_or_si128(_mm_subs_epu8(p3, p0),
                                     _mm_subs_epu8(p0, p3)),
                         _mm_or_si128(_mm_subs_epu8(q3, q0),
                                      _mm_subs_epu8(q0, q3)));
    flat = _mm_max_epu8(work, flat);
    work = _mm_max_epu8(_mm_or_si128(_mm_subs_epu8(p4, p0),
                                     _mm_subs_epu8(p0, p4)),
                         _mm_or_si128(_mm_subs_epu8(q4, q0),
                                      _mm_subs_epu8(q0, q4)));
    flat = _mm_max_epu8(work, flat);
    flat = _mm_subs_epu8(flat, one);
    flat = _mm_cmpeq_epi8(flat, zero);
    flat = _mm_and_si128(flat, mask);
  }
  {
    const __m128i four = _mm_set1_epi16(4);
    unsigned char *src = s;
    int i = 0;
    do {
      __m128i workp_a, workp_b, workp_shft;
      p4 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src - 5 * p)), zero);
      p3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src - 4 * p)), zero);
      p2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src - 3 * p)), zero);
      p1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src - 2 * p)), zero);
      p0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src - 1 * p)), zero);
      q0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src - 0 * p)), zero);
      q1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src + 1 * p)), zero);
      q2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src + 2 * p)), zero);
      q3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src + 3 * p)), zero);
      q4 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src + 4 * p)), zero);

      workp_a = _mm_add_epi16(_mm_add_epi16(p4, p3), _mm_add_epi16(p2, p1));
      workp_a = _mm_add_epi16(_mm_add_epi16(workp_a, four), p0);
      workp_b = _mm_add_epi16(_mm_add_epi16(q0, p2), p4);
      workp_shft = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
      _mm_storel_epi64((__m128i *)&flat_op2[i*8],
                       _mm_packus_epi16(workp_shft, workp_shft));

      workp_b = _mm_add_epi16(_mm_add_epi16(q0, q1), p1);
      workp_shft = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
      _mm_storel_epi64((__m128i *)&flat_op1[i*8],
                       _mm_packus_epi16(workp_shft, workp_shft));

      workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, p4), q2);
      workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, p1), p0);
      workp_shft = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
      _mm_storel_epi64((__m128i *)&flat_op0[i*8],
                       _mm_packus_epi16(workp_shft, workp_shft));

      workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, p3), q3);
      workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, p0), q0);
      workp_shft = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
      _mm_storel_epi64((__m128i *)&flat_oq0[i*8],
                       _mm_packus_epi16(workp_shft, workp_shft));

      workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, p2), q4);
      workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, q0), q1);
      workp_shft = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
      _mm_storel_epi64((__m128i *)&flat_oq1[i*8],
                       _mm_packus_epi16(workp_shft, workp_shft));

      workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, p1), q4);
      workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, q1), q2);
      workp_shft = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
      _mm_storel_epi64((__m128i *)&flat_oq2[i*8],
                       _mm_packus_epi16(workp_shft, workp_shft));

      src += 8;
    } while (++i < count);
  }
  // lp filter
  {
    const __m128i t4 = _mm_set1_epi8(4);
    const __m128i t3 = _mm_set1_epi8(3);
    const __m128i t80 = _mm_set1_epi8(0x80);
    const __m128i te0 = _mm_set1_epi8(0xe0);
    const __m128i t1f = _mm_set1_epi8(0x1f);
    const __m128i t1 = _mm_set1_epi8(0x1);
    const __m128i t7f = _mm_set1_epi8(0x7f);

    const __m128i ps1 = _mm_xor_si128(_mm_loadu_si128((__m128i *)(s - 2 * p)),
                                      t80);
    const __m128i ps0 = _mm_xor_si128(_mm_loadu_si128((__m128i *)(s - 1 * p)),
                                      t80);
    const __m128i qs0 = _mm_xor_si128(_mm_loadu_si128((__m128i *)(s + 0 * p)),
                                      t80);
    const __m128i qs1 = _mm_xor_si128(_mm_loadu_si128((__m128i *)(s + 1 * p)),
                                      t80);
    __m128i vp8_filt;
    __m128i work_a;
    __m128i filter1, filter2;

    vp8_filt = _mm_and_si128(_mm_subs_epi8(ps1, qs1), hev);
    work_a = _mm_subs_epi8(qs0, ps0);
    vp8_filt = _mm_adds_epi8(vp8_filt, work_a);
    vp8_filt = _mm_adds_epi8(vp8_filt, work_a);
    vp8_filt = _mm_adds_epi8(vp8_filt, work_a);
    /* (vp8_filter + 3 * (qs0 - ps0)) & mask */
    vp8_filt = _mm_and_si128(vp8_filt, mask);

    filter1 = _mm_adds_epi8(vp8_filt, t4);
    filter2 = _mm_adds_epi8(vp8_filt, t3);

    /* Filter1 >> 3 */
    work_a = _mm_cmpgt_epi8(zero, filter1);
    filter1 = _mm_srli_epi16(filter1, 3);
    work_a = _mm_and_si128(work_a, te0);
    filter1 = _mm_and_si128(filter1, t1f);
    filter1 = _mm_or_si128(filter1, work_a);

    /* Filter2 >> 3 */
    work_a = _mm_cmpgt_epi8(zero, filter2);
    filter2 = _mm_srli_epi16(filter2, 3);
    work_a = _mm_and_si128(work_a, te0);
    filter2 = _mm_and_si128(filter2, t1f);
    filter2 = _mm_or_si128(filter2, work_a);

    /* vp8_filt >> 1 */
    vp8_filt = _mm_adds_epi8(filter1, t1);
    work_a = _mm_cmpgt_epi8(zero, vp8_filt);
    vp8_filt = _mm_srli_epi16(vp8_filt, 1);
    work_a = _mm_and_si128(work_a, t80);
    vp8_filt = _mm_and_si128(vp8_filt, t7f);
    vp8_filt = _mm_or_si128(vp8_filt, work_a);

    vp8_filt = _mm_andnot_si128(hev, vp8_filt);

    work_a = _mm_xor_si128(_mm_subs_epi8(qs0, filter1), t80);
    q0 = _mm_load_si128((__m128i *)flat_oq0);
    work_a = _mm_andnot_si128(flat, work_a);
    q0 = _mm_and_si128(flat, q0);
    q0 = _mm_or_si128(work_a, q0);

    work_a = _mm_xor_si128(_mm_subs_epi8(qs1, vp8_filt), t80);
    q1 = _mm_load_si128((__m128i *)flat_oq1);
    work_a = _mm_andnot_si128(flat, work_a);
    q1 = _mm_and_si128(flat, q1);
    q1 = _mm_or_si128(work_a, q1);

    work_a = _mm_loadu_si128((__m128i *)(s + 2 * p));
    q2 = _mm_load_si128((__m128i *)flat_oq2);
    work_a = _mm_andnot_si128(flat, work_a);
    q2 = _mm_and_si128(flat, q2);
    q2 = _mm_or_si128(work_a, q2);

    work_a = _mm_xor_si128(_mm_adds_epi8(ps0, filter2), t80);
    p0 = _mm_load_si128((__m128i *)flat_op0);
    work_a = _mm_andnot_si128(flat, work_a);
    p0 = _mm_and_si128(flat, p0);
    p0 = _mm_or_si128(work_a, p0);

    work_a = _mm_xor_si128(_mm_adds_epi8(ps1, vp8_filt), t80);
    p1 = _mm_load_si128((__m128i *)flat_op1);
    work_a = _mm_andnot_si128(flat, work_a);
    p1 = _mm_and_si128(flat, p1);
    p1 = _mm_or_si128(work_a, p1);

    work_a = _mm_loadu_si128((__m128i *)(s - 3 * p));
    p2 = _mm_load_si128((__m128i *)flat_op2);
    work_a = _mm_andnot_si128(flat, work_a);
    p2 = _mm_and_si128(flat, p2);
    p2 = _mm_or_si128(work_a, p2);

    if (count == 1) {
      _mm_storel_epi64((__m128i *)(s - 3 * p), p2);
      _mm_storel_epi64((__m128i *)(s - 2 * p), p1);
      _mm_storel_epi64((__m128i *)(s - 1 * p), p0);
      _mm_storel_epi64((__m128i *)(s + 0 * p), q0);
      _mm_storel_epi64((__m128i *)(s + 1 * p), q1);
      _mm_storel_epi64((__m128i *)(s + 2 * p), q2);
    } else {
      _mm_storeu_si128((__m128i *)(s - 3 * p), p2);
      _mm_storeu_si128((__m128i *)(s - 2 * p), p1);
      _mm_storeu_si128((__m128i *)(s - 1 * p), p0);
      _mm_storeu_si128((__m128i *)(s + 0 * p), q0);
      _mm_storeu_si128((__m128i *)(s + 1 * p), q1);
      _mm_storeu_si128((__m128i *)(s + 2 * p), q2);
    }
  }
}

/* Horizontal MB filtering */
void vp8_loop_filter_mbh_sse2(unsigned char *y_ptr, unsigned char *u_ptr, unsigned char *v_ptr,
                              int y_stride, int uv_stride, struct loop_filter_info *lfi) {

  vp8_mbloop_filter_horizontal_edge_c_sse2(y_ptr, y_stride, lfi->mblim,
                                           lfi->lim, lfi->hev_thr, 2);

  /* TODO: write sse2 version with u,v interleaved */
  if (u_ptr)
    vp8_mbloop_filter_horizontal_edge_c_sse2(u_ptr, uv_stride, lfi->mblim,
                                             lfi->lim, lfi->hev_thr, 1);

  if (v_ptr)
    vp8_mbloop_filter_horizontal_edge_c_sse2(v_ptr, uv_stride, lfi->mblim,
                                             lfi->lim, lfi->hev_thr, 1);
}

void vp8_loop_filter_bh8x8_sse2(unsigned char *y_ptr, unsigned char *u_ptr,
                             unsigned char *v_ptr, int y_stride, int uv_stride,
                             struct loop_filter_info *lfi) {
  vp8_mbloop_filter_horizontal_edge_c_sse2(
    y_ptr + 8 * y_stride, y_stride, lfi->blim, lfi->lim, lfi->hev_thr, 2);
}

/* Vertical MB Filtering */
void vp8_loop_filter_mbv_sse2(unsigned char *y_ptr, unsigned char *u_ptr, unsigned char *v_ptr,
                              int y_stride, int uv_stride, struct loop_filter_info *lfi) {
  vp8_mbloop_filter_vertical_edge_sse2(y_ptr, y_stride, lfi->mblim, lfi->lim, lfi->hev_thr, 2);

  if (u_ptr)
    vp8_mbloop_filter_vertical_edge_uv_sse2(u_ptr, uv_stride, lfi->mblim, lfi->lim, lfi->hev_thr, v_ptr);
}


/* Horizontal B Filtering */
void vp8_loop_filter_bh_sse2(unsigned char *y_ptr, unsigned char *u_ptr, unsigned char *v_ptr,
                             int y_stride, int uv_stride, struct loop_filter_info *lfi) {
  vp8_loop_filter_horizontal_edge_sse2(y_ptr + 4 * y_stride, y_stride, lfi->blim, lfi->lim, lfi->hev_thr, 2);
  vp8_loop_filter_horizontal_edge_sse2(y_ptr + 8 * y_stride, y_stride, lfi->blim, lfi->lim, lfi->hev_thr, 2);
  vp8_loop_filter_horizontal_edge_sse2(y_ptr + 12 * y_stride, y_stride, lfi->blim, lfi->lim, lfi->hev_thr, 2);

  if (u_ptr)
    vp8_loop_filter_horizontal_edge_uv_sse2(u_ptr + 4 * uv_stride, uv_stride, lfi->blim, lfi->lim, lfi->hev_thr, v_ptr + 4 * uv_stride);
}


void vp8_loop_filter_bhs_sse2(unsigned char *y_ptr, int y_stride, const unsigned char *blimit) {
  vp8_loop_filter_simple_horizontal_edge_sse2(y_ptr + 4 * y_stride, y_stride, blimit);
  vp8_loop_filter_simple_horizontal_edge_sse2(y_ptr + 8 * y_stride, y_stride, blimit);
  vp8_loop_filter_simple_horizontal_edge_sse2(y_ptr + 12 * y_stride, y_stride, blimit);
}


/* Vertical B Filtering */
void vp8_loop_filter_bv_sse2(unsigned char *y_ptr, unsigned char *u_ptr, unsigned char *v_ptr,
                             int y_stride, int uv_stride, struct loop_filter_info *lfi) {
  vp8_loop_filter_vertical_edge_sse2(y_ptr + 4, y_stride, lfi->blim, lfi->lim, lfi->hev_thr, 2);
  vp8_loop_filter_vertical_edge_sse2(y_ptr + 8, y_stride, lfi->blim, lfi->lim, lfi->hev_thr, 2);
  vp8_loop_filter_vertical_edge_sse2(y_ptr + 12, y_stride, lfi->blim, lfi->lim, lfi->hev_thr, 2);

  if (u_ptr)
    vp8_loop_filter_vertical_edge_uv_sse2(u_ptr + 4, uv_stride, lfi->blim, lfi->lim, lfi->hev_thr, v_ptr + 4);
}


void vp8_loop_filter_bvs_sse2(unsigned char *y_ptr, int y_stride, const unsigned char *blimit) {
  vp8_loop_filter_simple_vertical_edge_sse2(y_ptr + 4, y_stride, blimit);
  vp8_loop_filter_simple_vertical_edge_sse2(y_ptr + 8, y_stride, blimit);
  vp8_loop_filter_simple_vertical_edge_sse2(y_ptr + 12, y_stride, blimit);
}

#endif
