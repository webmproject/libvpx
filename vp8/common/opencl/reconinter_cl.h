/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_RECONINTER_CL_H
#define __INC_RECONINTER_CL_H

#include "blockd_cl.h"
#include "subpixel_cl.h"
#include "filter_cl.h"

extern void vp8_build_inter_predictors_mb_cl(MACROBLOCKD *x);
extern void vp8_build_inter_predictors_mbuv_cl(MACROBLOCKD *x);

extern void vp8_build_inter_predictors_mb_s_cl(MACROBLOCKD *x);
//extern void vp8_build_inter_predictors_b_cl(BLOCKD *d, int pitch);

#endif
