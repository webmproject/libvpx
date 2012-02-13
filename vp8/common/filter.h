/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef FILTER_H
#define FILTER_H

#include "vpx_config.h"
#include "vpx_scale/yv12config.h"

#define BLOCK_HEIGHT_WIDTH 4
#define VP8_FILTER_WEIGHT 128
#define VP8_FILTER_SHIFT  7

/* whether to use a special filter for edge pixels */
#define EDGE_PIXEL_FILTER 0

extern const short vp8_bilinear_filters[8][2];
extern const short vp8_sub_pel_filters[8][INTERP_EXTEND*2];

#if EDGE_PIXEL_FILTER > 0
#define EDGE_PIXEL_FILTER_EXTEND 2
extern const short vp8_sub_pel_filters_ns[64][4*EDGE_PIXEL_FILTER_EXTEND*EDGE_PIXEL_FILTER_EXTEND];
#endif

#endif //FILTER_H
