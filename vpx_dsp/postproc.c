/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <stdlib.h>

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"

#include "vpx/vpx_integer.h"
#include "vpx_ports/mem.h"

void vpx_plane_add_noise_c(uint8_t *start, char *noise,
                           char blackclamp[16],
                           char whiteclamp[16],
                           char bothclamp[16],
                           unsigned int width, unsigned int height, int pitch) {
  unsigned int i, j;

  // TODO(jbb): why does simd code use both but c doesn't,  normalize and
  // fix..
  (void) bothclamp;
  for (i = 0; i < height; i++) {
    uint8_t *pos = start + i * pitch;
    char  *ref = (char *)(noise + (rand() & 0xff));  // NOLINT

    for (j = 0; j < width; j++) {
      if (pos[j] < blackclamp[0])
        pos[j] = blackclamp[0];

      if (pos[j] > 255 - whiteclamp[0])
        pos[j] = 255 - whiteclamp[0];

      pos[j] += ref[j];
    }
  }
}
