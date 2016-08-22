/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_CORNER_DETECT_H_
#define VP10_ENCODER_CORNER_DETECT_H_

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

int FastCornerDetect(unsigned char *buf, int width, int height, int stride,
                     int *points, int max_points);

#endif  // VP10_ENCODER_CORNER_DETECT_H
