/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/common/mips/dspr2/vp9_common_dspr2.h"
#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_ports/mem.h"

#if HAVE_DSPR2
uint8_t vp9_ff_cropTbl_a[256 + 2 * CROP_WIDTH];
uint8_t *vp9_ff_cropTbl;

void vp9_dsputil_static_init(void) {
  int i;

  for (i = 0; i < 256; i++) vp9_ff_cropTbl_a[i + CROP_WIDTH] = i;

  for (i = 0; i < CROP_WIDTH; i++) {
    vp9_ff_cropTbl_a[i] = 0;
    vp9_ff_cropTbl_a[i + CROP_WIDTH + 256] = 255;
  }

  vp9_ff_cropTbl = &vp9_ff_cropTbl_a[CROP_WIDTH];
}

#endif
