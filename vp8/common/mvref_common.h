/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "onyxc_int.h"
#include "blockd.h"

// MR reference entropy header file.
#if CONFIG_NEWBESTREFMV

#ifndef __INC_MVREF_COMMON_H
#define __INC_MVREF_COMMON_H

unsigned int mv_distance(int_mv *mv1, int_mv *mv2);

void find_mv_refs(
  MACROBLOCKD *xd,
  MODE_INFO *here,
  MODE_INFO *lf_here,
  MV_REFERENCE_FRAME ref_frame,
  int_mv * mv_ref_list,
  int *ref_sign_bias
);

#endif

#endif
