/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_ENCODER_CORNER_MATCH_H_
#define AV1_ENCODER_CORNER_MATCH_H_

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

typedef struct {
  double x, y;
  double rx, ry;
} Correspondence;

int determine_correspondence(unsigned char *frm, int *frm_corners,
                             int num_frm_corners, unsigned char *ref,
                             int *ref_corners, int num_ref_corners, int width,
                             int height, int frm_stride, int ref_stride,
                             double *correspondence_pts);

#endif  // AV1_ENCODER_CORNER_MATCH_H_
