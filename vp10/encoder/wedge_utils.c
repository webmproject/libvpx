/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>

#include "vpx/vpx_integer.h"

#include "vpx_ports/mem.h"

#include "vpx_dsp/vpx_dsp_common.h"

#include "vp10/common/reconinter.h"

#define MAX_MASK_VALUE  (1 << WEDGE_WEIGHT_BITS)

/**
 * Computes SSE of a compound predictor constructed from 2 fundamental
 * predictors p0 and p1 using blending with mask.
 *
 * r1:  Residuals of p1.
 *      (source - p1)
 * d:   Difference of p1 and p0.
 *      (p1 - p0)
 * m:   The blending mask
 * N:   Number of pixels
 *
 * 'r1', 'd', and 'm' are contiguous.
 *
 * Computes:
 *  Sum((MAX_MASK_VALUE*r1 + mask*d)**2), which is equivalent to:
 *  Sum((mask*r0 + (MAX_MASK_VALUE-mask)*r1)**2),
 *    where r0 is (source - p0), and r1 is (source - p1), which is in turn
 *    is equivalent to:
 *  Sum((source*MAX_MASK_VALUE - (mask*p0 + (MAX_MASK_VALUE-mask)*p1))**2),
 *    which is the SSE of the residuals of the compound predictor scaled up by
 *    MAX_MASK_VALUE**2.
 *
 * Note that we clamp the partial term in the loop to 16 bits signed. This is
 * to facilitate equivalent SIMD implementation. It should have no effect if
 * residuals are within 16 - WEDGE_WEIGHT_BITS (=10) signed, which always
 * holds for 8 bit input, and on real input, it should hold practically always,
 * as residuals are expected to be small.
 */
uint64_t vp10_wedge_sse_from_residuals_c(const int16_t *r1,
                                         const int16_t *d,
                                         const uint8_t *m,
                                         int N) {
  uint64_t csse = 0;
  int i;
  assert(N % 64 == 0);
  for (i = 0 ; i < N ; i++) {
    int32_t t = MAX_MASK_VALUE*r1[i] + m[i]*d[i];
    t = clamp(t, INT16_MIN, INT16_MAX);
    csse += t*t;
  }
  return ROUND_POWER_OF_TWO(csse, 2 * WEDGE_WEIGHT_BITS);
}

/**
 * Choose the mask sign for a compound predictor.
 *
 * ds:    Difference of the squares of the residuals.
 *        r0**2 - r1**2
 * m:     The blending mask
 * N:     Number of pixels
 * limit: Pre-computed threshold value.
 *        MAX_MASK_VALUE/2 * (sum(r0**2) - sum(r1**2))
 *
 * 'ds' and 'm' are contiguous.
 *
 * Returns true if the negated mask has lower SSE compared to the positive
 * mask. Computation is based on:
 *  Sum((mask*r0 + (MAX_MASK_VALUE-mask)*r1)**2)
 *                                     >
 *                                Sum(((MAX_MASK_VALUE-mask)*r0 + mask*r1)**2)
 *
 *  which can be simplified to:
 *
 *  Sum(mask*(r0**2 - r1**2)) > MAX_MASK_VALUE/2 * (sum(r0**2) - sum(r1**2))
 *
 *  The right hand side does not depend on the mask, and needs to be passed as
 *  the 'limit' parameter.
 *
 *  After pre-computing (r0**2 - r1**2), which is passed in as 'ds', the left
 *  hand side is simply a scalar product between an int16_t and uint8_t vector.
 *
 *  Note that for efficiency, ds is stored on 16 bits. Real input residuals
 *  being small, this should not cause a noticeable issue.
 */
int vp10_wedge_sign_from_residuals_c(const int16_t *ds,
                                     const uint8_t *m,
                                     int N,
                                     int64_t limit) {
  int64_t acc = 0;

  assert(N % 64 == 0);

  do {
    acc += *ds++ * *m++;
  } while (--N);

  return acc > limit;
}

/**
 * Compute the element-wise difference of the squares of 2 arrays.
 *
 * d: Difference of the squares of the inputs: a**2 - b**2
 * a: First input array
 * b: Second input array
 * N: Number of elements
 *
 * 'd', 'a', and 'b' are contiguous.
 *
 * The result is saturated to signed 16 bits.
 */
void vp10_wedge_compute_delta_squares_c(int16_t *d,
                                        const int16_t *a,
                                        const int16_t *b,
                                        int N) {
  int i;

  assert(N % 64 == 0);

  for (i = 0 ; i < N ; i++)
    d[i] = clamp(a[i]*a[i] - b[i]*b[i], INT16_MIN, INT16_MAX);
}

