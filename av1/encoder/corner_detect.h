/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_ENCODER_CORNER_DETECT_H_
#define AV1_ENCODER_CORNER_DETECT_H_

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

int fast_corner_detect(unsigned char *buf, int width, int height, int stride,
                       int *points, int max_points);

#endif  // AV1_ENCODER_CORNER_DETECT_H_
