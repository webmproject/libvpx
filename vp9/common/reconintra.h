/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef __INC_RECONINTRA_H
#define __INC_RECONINTRA_H

#include "blockd.h"

extern void vp9_recon_intra_mbuv(MACROBLOCKD *xd);
extern B_PREDICTION_MODE vp9_find_dominant_direction(unsigned char *ptr,
                                                     int stride, int n);
extern B_PREDICTION_MODE vp9_find_bpred_context(BLOCKD *x);

#endif  // __INC_RECONINTRA_H
