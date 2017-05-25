/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp8/common/skin_detection.h"

#define MODEL_MODE 1

// Fixed-point skin color model parameters.
static const int skin_mean[5][2] = { { 7463, 9614 },
                                     { 6400, 10240 },
                                     { 7040, 10240 },
                                     { 8320, 9280 },
                                     { 6800, 9614 } };
static const int skin_inv_cov[4] = { 4107, 1663, 1663, 2157 };  // q16
static const int skin_threshold[6] = { 1570636, 1400000, 800000,
                                       800000,  800000,  800000 };  // q18

// Thresholds on luminance.
static const int y_low = 40;
static const int y_high = 220;

// Evaluates the Mahalanobis distance measure for the input CbCr values.
static int evaluate_skin_color_difference(int cb, int cr, int idx) {
  const int cb_q6 = cb << 6;
  const int cr_q6 = cr << 6;
  const int cb_diff_q12 =
      (cb_q6 - skin_mean[idx][0]) * (cb_q6 - skin_mean[idx][0]);
  const int cbcr_diff_q12 =
      (cb_q6 - skin_mean[idx][0]) * (cr_q6 - skin_mean[idx][1]);
  const int cr_diff_q12 =
      (cr_q6 - skin_mean[idx][1]) * (cr_q6 - skin_mean[idx][1]);
  const int cb_diff_q2 = (cb_diff_q12 + (1 << 9)) >> 10;
  const int cbcr_diff_q2 = (cbcr_diff_q12 + (1 << 9)) >> 10;
  const int cr_diff_q2 = (cr_diff_q12 + (1 << 9)) >> 10;
  const int skin_diff =
      skin_inv_cov[0] * cb_diff_q2 + skin_inv_cov[1] * cbcr_diff_q2 +
      skin_inv_cov[2] * cbcr_diff_q2 + skin_inv_cov[3] * cr_diff_q2;
  return skin_diff;
}

// Checks if the input yCbCr values corresponds to skin color.
int skin_pixel(int y, int cb, int cr, int motion) {
  if (y < y_low || y > y_high) {
    return 0;
  } else {
    if (MODEL_MODE == 0) {
      return (evaluate_skin_color_difference(cb, cr, 0) < skin_threshold[0]);
    } else {
      int i = 0;
      // Exit on grey.
      if (cb == 128 && cr == 128) return 0;
      // Exit on very strong cb.
      if (cb > 150 && cr < 110) return 0;
      for (; i < 5; ++i) {
        int skin_color_diff = evaluate_skin_color_difference(cb, cr, i);
        if (skin_color_diff < skin_threshold[i + 1]) {
          if (y < 60 && skin_color_diff > 3 * (skin_threshold[i + 1] >> 2)) {
            return 0;
          } else if (motion == 0 &&
                     skin_color_diff > (skin_threshold[i + 1] >> 1)) {
            return 0;
          } else {
            return 1;
          }
        }
        // Exit if difference is much large than the threshold.
        if (skin_color_diff > (skin_threshold[i + 1] << 3)) {
          return 0;
        }
      }
      return 0;
    }
  }
}

int compute_skin_block(const uint8_t *y, const uint8_t *u, const uint8_t *v,
                       int stride, int strideuv, int consec_zeromv,
                       int curr_motion_magn) {
  // No skin if block has been zero/small motion for long consecutive time.
  if (consec_zeromv > 60 && curr_motion_magn == 0) {
    return 0;
  } else {
    int motion = 1;
    // Take center pixel in block to determine is_skin.
    const int ysource = (y[7 * stride + 7] + y[7 * stride + 8] +
                         y[8 * stride + 7] + y[8 * stride + 8]) >>
                        2;
    const int usource = (u[3 * strideuv + 3] + u[3 * strideuv + 4] +
                         u[4 * strideuv + 3] + u[4 * strideuv + 4]) >>
                        2;
    const int vsource = (v[3 * strideuv + 3] + v[3 * strideuv + 4] +
                         v[4 * strideuv + 3] + v[4 * strideuv + 4]) >>
                        2;
    if (consec_zeromv > 25 && curr_motion_magn == 0) motion = 0;
    return skin_pixel(ysource, usource, vsource, motion);
  }
}
