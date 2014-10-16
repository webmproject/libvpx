/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include "vp9/common/vp9_common.h"
#include "vpx_ports/mem.h"

unsigned int vp9_avg_8x8_c(const uint8_t *s, int p) {
  int i, j;
  int sum = 0;
  for (i = 0; i < 8; ++i, s+=p)
    for (j = 0; j < 8; sum += s[j], ++j) {}

  return (sum + 32) >> 6;
}

#if CONFIG_VP9_HIGHBITDEPTH
unsigned int vp9_highbd_avg_8x8_c(const uint8_t *s8, int p) {
  int i, j;
  int sum = 0;
  const uint16_t* s = CONVERT_TO_SHORTPTR(s8);
  for (i = 0; i < 8; ++i, s+=p)
    for (j = 0; j < 8; sum += s[j], ++j) {}

  return (sum + 32) >> 6;
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

