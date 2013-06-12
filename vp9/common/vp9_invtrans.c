/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/common/vp9_invtrans.h"
#include "./vp9_rtcd.h"

void vp9_inverse_transform_b_4x4_add(MACROBLOCKD *xd, int eob, int16_t *dqcoeff,
                                     uint8_t *dest, int stride) {
  if (eob <= 1)
    xd->inv_txm4x4_1_add(dqcoeff, dest, stride);
  else
    xd->inv_txm4x4_add(dqcoeff, dest, stride);
}
