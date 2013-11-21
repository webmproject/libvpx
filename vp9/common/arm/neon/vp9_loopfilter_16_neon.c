/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vp9_rtcd.h"

void vp9_mbloop_filter_horizontal_edge_16_neon(uint8_t *s, int p /* pitch */,
                                               const uint8_t *blimit0,
                                               const uint8_t *limit0,
                                               const uint8_t *thresh0,
                                               const uint8_t *blimit1,
                                               const uint8_t *limit1,
                                               const uint8_t *thresh1) {
  vp9_mbloop_filter_horizontal_edge(s, p, blimit0, limit0, thresh0, 1);
  vp9_mbloop_filter_horizontal_edge(s + 8, p, blimit1, limit1, thresh1, 1);
}
