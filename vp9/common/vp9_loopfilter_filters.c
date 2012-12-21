/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <stdlib.h>
#include "vpx_config.h"
#include "vp9/common/vp9_loopfilter.h"
#include "vp9/common/vp9_onyxc_int.h"

typedef unsigned char uc;

static __inline signed char signed_char_clamp(int t) {
  t = (t < -128 ? -128 : t);
  t = (t > 127 ? 127 : t);
  return (signed char) t;
}


/* should we apply any filter at all ( 11111111 yes, 00000000 no) */
static __inline signed char filter_mask(uc limit, uc blimit,
                                        uc p3, uc p2, uc p1, uc p0,
                                        uc q0, uc q1, uc q2, uc q3) {
  signed char mask = 0;
  mask |= (abs(p3 - p2) > limit) * -1;
  mask |= (abs(p2 - p1) > limit) * -1;
  mask |= (abs(p1 - p0) > limit) * -1;
  mask |= (abs(q1 - q0) > limit) * -1;
  mask |= (abs(q2 - q1) > limit) * -1;
  mask |= (abs(q3 - q2) > limit) * -1;
  mask |= (abs(p0 - q0) * 2 + abs(p1 - q1) / 2  > blimit) * -1;
  mask = ~mask;
  return mask;
}

/* is there high variance internal edge ( 11111111 yes, 00000000 no) */
static __inline signed char hevmask(uc thresh, uc p1, uc p0, uc q0, uc q1) {
  signed char hev = 0;
  hev  |= (abs(p1 - p0) > thresh) * -1;
  hev  |= (abs(q1 - q0) > thresh) * -1;
  return hev;
}

static __inline void filter(signed char mask, uc hev, uc *op1,
                            uc *op0, uc *oq0, uc *oq1)

{
  signed char ps0, qs0;
  signed char ps1, qs1;
  signed char filter, Filter1, Filter2;
  signed char u;

  ps1 = (signed char) * op1 ^ 0x80;
  ps0 = (signed char) * op0 ^ 0x80;
  qs0 = (signed char) * oq0 ^ 0x80;
  qs1 = (signed char) * oq1 ^ 0x80;

  /* add outer taps if we have high edge variance */
  filter = signed_char_clamp(ps1 - qs1);
  filter &= hev;

  /* inner taps */
  filter = signed_char_clamp(filter + 3 * (qs0 - ps0));
  filter &= mask;

  /* save bottom 3 bits so that we round one side +4 and the other +3
   * if it equals 4 we'll set to adjust by -1 to account for the fact
   * we'd round 3 the other way
   */
  Filter1 = signed_char_clamp(filter + 4);
  Filter2 = signed_char_clamp(filter + 3);
  Filter1 >>= 3;
  Filter2 >>= 3;
  u = signed_char_clamp(qs0 - Filter1);
  *oq0 = u ^ 0x80;
  u = signed_char_clamp(ps0 + Filter2);
  *op0 = u ^ 0x80;
  filter = Filter1;

  /* outer tap adjustments */
  filter += 1;
  filter >>= 1;
  filter &= ~hev;

  u = signed_char_clamp(qs1 - filter);
  *oq1 = u ^ 0x80;
  u = signed_char_clamp(ps1 + filter);
  *op1 = u ^ 0x80;

}

void vp9_loop_filter_horizontal_edge_c
(
  unsigned char *s,
  int p, /* pitch */
  const unsigned char *blimit,
  const unsigned char *limit,
  const unsigned char *thresh,
  int count
) {
  int  hev = 0; /* high edge variance */
  signed char mask = 0;
  int i = 0;

  /* loop filter designed to work using chars so that we can make maximum use
   * of 8 bit simd instructions.
   */
  do {
    mask = filter_mask(limit[0], blimit[0],
                       s[-4 * p], s[-3 * p], s[-2 * p], s[-1 * p],
                       s[0 * p], s[1 * p], s[2 * p], s[3 * p]);

    hev = hevmask(thresh[0], s[-2 * p], s[-1 * p], s[0 * p], s[1 * p]);

    filter(mask, hev, s - 2 * p, s - 1 * p, s, s + 1 * p);

    ++s;
  } while (++i < count * 8);
}

void vp9_loop_filter_vertical_edge_c(unsigned char *s,
                                     int p,
                                     const unsigned char *blimit,
                                     const unsigned char *limit,
                                     const unsigned char *thresh,
                                     int count) {
  int  hev = 0; /* high edge variance */
  signed char mask = 0;
  int i = 0;

  /* loop filter designed to work using chars so that we can make maximum use
   * of 8 bit simd instructions.
   */
  do {
    mask = filter_mask(limit[0], blimit[0],
                       s[-4], s[-3], s[-2], s[-1],
                       s[0], s[1], s[2], s[3]);

    hev = hevmask(thresh[0], s[-2], s[-1], s[0], s[1]);

    filter(mask, hev, s - 2, s - 1, s, s + 1);

    s += p;
  } while (++i < count * 8);
}
static __inline signed char flatmask(uc thresh,
                                     uc p4, uc p3, uc p2, uc p1, uc p0,
                                     uc q0, uc q1, uc q2, uc q3, uc q4) {
  signed char flat = 0;
  flat |= (abs(p1 - p0) > 1) * -1;
  flat |= (abs(q1 - q0) > 1) * -1;
  flat |= (abs(p0 - p2) > 1) * -1;
  flat |= (abs(q0 - q2) > 1) * -1;
  flat |= (abs(p3 - p0) > 1) * -1;
  flat |= (abs(q3 - q0) > 1) * -1;
  flat |= (abs(p4 - p0) > 1) * -1;
  flat |= (abs(q4 - q0) > 1) * -1;
  flat = ~flat;
  return flat;
}

static __inline void mbfilter(signed char mask, uc hev, uc flat,
                              uc *op4, uc *op3, uc *op2, uc *op1, uc *op0,
                              uc *oq0, uc *oq1, uc *oq2, uc *oq3, uc *oq4) {
  /* use a 7 tap filter [1, 1, 1, 2, 1, 1, 1] for flat line */
  if (flat && mask) {
    unsigned char p0, q0;
    unsigned char p1, q1;
    unsigned char p2, q2;
    unsigned char p3, q3;
    unsigned char p4, q4;

    p4 = *op4;
    p3 = *op3;
    p2 = *op2;
    p1 = *op1;
    p0 = *op0;
    q0 = *oq0;
    q1 = *oq1;
    q2 = *oq2;
    q3 = *oq3;
    q4 = *oq4;

    *op2 = (p4 + p4 + p3 + p2 + p2 + p1 + p0 + q0 + 4) >> 3;
    *op1 = (p4 + p3 + p2 + p1 + p1 + p0 + q0 + q1 + 4) >> 3;
    *op0 = (p3 + p2 + p1 + p0 + p0 + q0 + q1 + q2 + 4) >> 3;
    *oq0 = (p2 + p1 + p0 + q0 + q0 + q1 + q2 + q3 + 4) >> 3;
    *oq1 = (p1 + p0 + q0 + q1 + q1 + q2 + q3 + q4 + 4) >> 3;
    *oq2 = (p0 + q0 + q1 + q2 + q2 + q3 + q4 + q4 + 4) >> 3;
  } else {
    signed char ps0, qs0;
    signed char ps1, qs1;
    signed char filter, Filter1, Filter2;
    signed char u;

    ps1 = (signed char) * op1 ^ 0x80;
    ps0 = (signed char) * op0 ^ 0x80;
    qs0 = (signed char) * oq0 ^ 0x80;
    qs1 = (signed char) * oq1 ^ 0x80;

    /* add outer taps if we have high edge variance */
    filter = signed_char_clamp(ps1 - qs1);
    filter &= hev;

    /* inner taps */
    filter = signed_char_clamp(filter + 3 * (qs0 - ps0));
    filter &= mask;

    Filter1 = signed_char_clamp(filter + 4);
    Filter2 = signed_char_clamp(filter + 3);
    Filter1 >>= 3;
    Filter2 >>= 3;

    u = signed_char_clamp(qs0 - Filter1);
    *oq0 = u ^ 0x80;
    u = signed_char_clamp(ps0 + Filter2);
    *op0 = u ^ 0x80;
    filter = Filter1;

    /* outer tap adjustments */
    filter += 1;
    filter >>= 1;
    filter &= ~hev;

    u = signed_char_clamp(qs1 - filter);
    *oq1 = u ^ 0x80;
    u = signed_char_clamp(ps1 + filter);
    *op1 = u ^ 0x80;
  }
}
void vp9_mbloop_filter_horizontal_edge_c
(
  unsigned char *s,
  int p,
  const unsigned char *blimit,
  const unsigned char *limit,
  const unsigned char *thresh,
  int count
) {
  signed char hev = 0; /* high edge variance */
  signed char mask = 0;
  signed char flat = 0;
  int i = 0;

  /* loop filter designed to work using chars so that we can make maximum use
   * of 8 bit simd instructions.
   */
  do {

    mask = filter_mask(limit[0], blimit[0],
                       s[-4 * p], s[-3 * p], s[-2 * p], s[-1 * p],
                       s[ 0 * p], s[ 1 * p], s[ 2 * p], s[ 3 * p]);

    hev = hevmask(thresh[0], s[-2 * p], s[-1 * p], s[0 * p], s[1 * p]);

    flat = flatmask(thresh[0],
                    s[-5 * p], s[-4 * p], s[-3 * p], s[-2 * p], s[-1 * p],
                    s[ 0 * p], s[ 1 * p], s[ 2 * p], s[ 3 * p], s[ 4 * p]);
    mbfilter(mask, hev, flat,
             s - 5 * p, s - 4 * p, s - 3 * p, s - 2 * p, s - 1 * p,
             s,       s + 1 * p, s + 2 * p, s + 3 * p, s + 4 * p);

    ++s;
  } while (++i < count * 8);

}
void vp9_mbloop_filter_vertical_edge_c
(
  unsigned char *s,
  int p,
  const unsigned char *blimit,
  const unsigned char *limit,
  const unsigned char *thresh,
  int count
) {
  signed char hev = 0; /* high edge variance */
  signed char mask = 0;
  signed char flat = 0;
  int i = 0;

  do {

    mask = filter_mask(limit[0], blimit[0],
                       s[-4], s[-3], s[-2], s[-1],
                       s[0], s[1], s[2], s[3]);

    hev = hevmask(thresh[0], s[-2], s[-1], s[0], s[1]);
    flat = flatmask(thresh[0],
                    s[-5], s[-4], s[-3], s[-2], s[-1],
                    s[ 0], s[ 1], s[ 2], s[ 3], s[ 4]);
    mbfilter(mask, hev, flat,
             s - 5, s - 4, s - 3, s - 2, s - 1,
             s,     s + 1, s + 2, s + 3, s + 4);
    s += p;
  } while (++i < count * 8);

}

/* should we apply any filter at all ( 11111111 yes, 00000000 no) */
static __inline signed char simple_filter_mask(uc blimit,
                                               uc p1, uc p0,
                                               uc q0, uc q1) {
  /* Why does this cause problems for win32?
   * error C2143: syntax error : missing ';' before 'type'
   *  (void) limit;
   */
  signed char mask = (abs(p0 - q0) * 2 + abs(p1 - q1) / 2  <= blimit) * -1;
  return mask;
}

static __inline void simple_filter(signed char mask,
                                   uc *op1, uc *op0,
                                   uc *oq0, uc *oq1) {
  signed char filter, Filter1, Filter2;
  signed char p1 = (signed char) * op1 ^ 0x80;
  signed char p0 = (signed char) * op0 ^ 0x80;
  signed char q0 = (signed char) * oq0 ^ 0x80;
  signed char q1 = (signed char) * oq1 ^ 0x80;
  signed char u;

  filter = signed_char_clamp(p1 - q1);
  filter = signed_char_clamp(filter + 3 * (q0 - p0));
  filter &= mask;

  /* save bottom 3 bits so that we round one side +4 and the other +3 */
  Filter1 = signed_char_clamp(filter + 4);
  Filter1 >>= 3;
  u = signed_char_clamp(q0 - Filter1);
  *oq0  = u ^ 0x80;

  Filter2 = signed_char_clamp(filter + 3);
  Filter2 >>= 3;
  u = signed_char_clamp(p0 + Filter2);
  *op0 = u ^ 0x80;
}

void vp9_loop_filter_simple_horizontal_edge_c
(
  unsigned char *s,
  int p,
  const unsigned char *blimit
) {
  signed char mask = 0;
  int i = 0;

  do {
    mask = simple_filter_mask(blimit[0],
                              s[-2 * p], s[-1 * p],
                              s[0 * p], s[1 * p]);
    simple_filter(mask,
                  s - 2 * p, s - 1 * p,
                  s, s + 1 * p);
    ++s;
  } while (++i < 16);
}

void vp9_loop_filter_simple_vertical_edge_c
(
  unsigned char *s,
  int p,
  const unsigned char *blimit
) {
  signed char mask = 0;
  int i = 0;

  do {
    mask = simple_filter_mask(blimit[0], s[-2], s[-1], s[0], s[1]);
    simple_filter(mask, s - 2, s - 1, s, s + 1);
    s += p;
  } while (++i < 16);

}

/* Vertical MB Filtering */
void vp9_loop_filter_mbv_c(unsigned char *y_ptr, unsigned char *u_ptr,
                           unsigned char *v_ptr, int y_stride, int uv_stride,
                           struct loop_filter_info *lfi) {
  vp9_mbloop_filter_vertical_edge_c(y_ptr, y_stride,
                                    lfi->mblim, lfi->lim, lfi->hev_thr, 2);

  if (u_ptr)
    vp9_mbloop_filter_vertical_edge_c(u_ptr, uv_stride,
                                      lfi->mblim, lfi->lim, lfi->hev_thr, 1);

  if (v_ptr)
    vp9_mbloop_filter_vertical_edge_c(v_ptr, uv_stride,
                                      lfi->mblim, lfi->lim, lfi->hev_thr, 1);
}

/* Vertical B Filtering */
void vp9_loop_filter_bv_c(unsigned char *y_ptr, unsigned char *u_ptr,
                          unsigned char *v_ptr, int y_stride, int uv_stride,
                          struct loop_filter_info *lfi) {
  vp9_loop_filter_vertical_edge_c(y_ptr + 4, y_stride,
                                  lfi->blim, lfi->lim, lfi->hev_thr, 2);
  vp9_loop_filter_vertical_edge_c(y_ptr + 8, y_stride,
                                  lfi->blim, lfi->lim, lfi->hev_thr, 2);
  vp9_loop_filter_vertical_edge_c(y_ptr + 12, y_stride,
                                  lfi->blim, lfi->lim, lfi->hev_thr, 2);

  if (u_ptr)
    vp9_loop_filter_vertical_edge_c(u_ptr + 4, uv_stride,
                                    lfi->blim, lfi->lim, lfi->hev_thr, 1);

  if (v_ptr)
    vp9_loop_filter_vertical_edge_c(v_ptr + 4, uv_stride,
                                    lfi->blim, lfi->lim, lfi->hev_thr, 1);
}

/* Horizontal MB filtering */
void vp9_loop_filter_mbh_c(unsigned char *y_ptr, unsigned char *u_ptr,
                           unsigned char *v_ptr, int y_stride, int uv_stride,
                           struct loop_filter_info *lfi) {
  vp9_mbloop_filter_horizontal_edge_c(y_ptr, y_stride,
                                      lfi->mblim, lfi->lim, lfi->hev_thr, 2);

  if (u_ptr)
    vp9_mbloop_filter_horizontal_edge_c(u_ptr, uv_stride,
                                        lfi->mblim, lfi->lim, lfi->hev_thr, 1);

  if (v_ptr)
    vp9_mbloop_filter_horizontal_edge_c(v_ptr, uv_stride,
                                        lfi->mblim, lfi->lim, lfi->hev_thr, 1);
}

/* Horizontal B Filtering */
void vp9_loop_filter_bh_c(unsigned char *y_ptr, unsigned char *u_ptr,
                          unsigned char *v_ptr, int y_stride, int uv_stride,
                          struct loop_filter_info *lfi) {
  vp9_loop_filter_horizontal_edge_c(y_ptr + 4 * y_stride, y_stride,
                                    lfi->blim, lfi->lim, lfi->hev_thr, 2);
  vp9_loop_filter_horizontal_edge_c(y_ptr + 8 * y_stride, y_stride,
                                    lfi->blim, lfi->lim, lfi->hev_thr, 2);
  vp9_loop_filter_horizontal_edge_c(y_ptr + 12 * y_stride, y_stride,
                                    lfi->blim, lfi->lim, lfi->hev_thr, 2);

  if (u_ptr)
    vp9_loop_filter_horizontal_edge_c(u_ptr + 4 * uv_stride, uv_stride,
                                      lfi->blim, lfi->lim, lfi->hev_thr, 1);

  if (v_ptr)
    vp9_loop_filter_horizontal_edge_c(v_ptr + 4 * uv_stride, uv_stride,
                                      lfi->blim, lfi->lim, lfi->hev_thr, 1);
}

void vp9_loop_filter_bh8x8_c(unsigned char *y_ptr, unsigned char *u_ptr,
                             unsigned char *v_ptr, int y_stride, int uv_stride,
                             struct loop_filter_info *lfi) {
  vp9_mbloop_filter_horizontal_edge_c(
    y_ptr + 8 * y_stride, y_stride, lfi->blim, lfi->lim, lfi->hev_thr, 2);
}

void vp9_loop_filter_bhs_c(unsigned char *y_ptr, int y_stride,
                           const unsigned char *blimit) {
  vp9_loop_filter_simple_horizontal_edge_c(y_ptr + 4 * y_stride,
                                           y_stride, blimit);
  vp9_loop_filter_simple_horizontal_edge_c(y_ptr + 8 * y_stride,
                                           y_stride, blimit);
  vp9_loop_filter_simple_horizontal_edge_c(y_ptr + 12 * y_stride,
                                           y_stride, blimit);
}

void vp9_loop_filter_bv8x8_c(unsigned char *y_ptr, unsigned char *u_ptr,
                             unsigned char *v_ptr, int y_stride, int uv_stride,
                             struct loop_filter_info *lfi) {
  vp9_mbloop_filter_vertical_edge_c(
    y_ptr + 8, y_stride, lfi->blim, lfi->lim, lfi->hev_thr, 2);
}

void vp9_loop_filter_bvs_c(unsigned char *y_ptr, int y_stride,
                           const unsigned char *blimit) {
  vp9_loop_filter_simple_vertical_edge_c(y_ptr + 4, y_stride, blimit);
  vp9_loop_filter_simple_vertical_edge_c(y_ptr + 8, y_stride, blimit);
  vp9_loop_filter_simple_vertical_edge_c(y_ptr + 12, y_stride, blimit);
}
