/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "aom_dsp/mips/common_dspr2.h"

#if HAVE_DSPR2
uint8_t aom_ff_cropTbl_a[256 + 2 * CROP_WIDTH];
uint8_t *aom_ff_cropTbl;

void aom_dsputil_static_init(void) {
  int i;

  for (i = 0; i < 256; i++) aom_ff_cropTbl_a[i + CROP_WIDTH] = i;

  for (i = 0; i < CROP_WIDTH; i++) {
    aom_ff_cropTbl_a[i] = 0;
    aom_ff_cropTbl_a[i + CROP_WIDTH + 256] = 255;
  }

  aom_ff_cropTbl = &aom_ff_cropTbl_a[CROP_WIDTH];
}

#endif
