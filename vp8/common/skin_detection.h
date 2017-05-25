/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP8_ENCODER_SKIN_DETECTION_H_
#define VP8_ENCODER_SKIN_DETECTION_H_

#include "vpx/vpx_integer.h"

#ifdef __cplusplus
extern "C" {
#endif

int skin_pixel(int y, int cb, int cr, int motion);

int compute_skin_block(const uint8_t *y, const uint8_t *u, const uint8_t *v,
                       int stride, int strideuv, int consec_zeromv,
                       int curr_motion_magn);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP8_ENCODER_SKIN_DETECTION_H_
