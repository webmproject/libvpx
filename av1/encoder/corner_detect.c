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
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <assert.h>

#include "av1/encoder/corner_detect.h"
#include "third_party/fastfeat/fast.h"

// Fast_9 wrapper
#define FAST_BARRIER 40
int FastCornerDetect(unsigned char *buf, int width, int height, int stride,
                     int *points, int max_points) {
  int num_points;
  xy *frm_corners_xy = fast9_detect_nonmax(buf, width, height, stride,
                                           FAST_BARRIER, &num_points);
  num_points = (num_points <= max_points ? num_points : max_points);
  if (num_points > 0 && frm_corners_xy) {
    memcpy(points, frm_corners_xy, sizeof(xy) * num_points);
    free(frm_corners_xy);
    return num_points;
  } else {
    return 0;
  }
}
