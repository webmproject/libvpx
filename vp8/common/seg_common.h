/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "type_aliases.h"
#include "vp8/common/blockd.h"

#ifndef __INC_SEG_COMMON_H__
#define __INC_SEG_COMMON_H__ 1

int segfeature_active( MACROBLOCKD *xd,
                       int segment_id,
                       SEG_LVL_FEATURES feature_id );

void enable_segfeature( MACROBLOCKD *xd,
                        int segment_id,
                        SEG_LVL_FEATURES feature_id );

void disable_segfeature( MACROBLOCKD *xd,
                         int segment_id,
                         SEG_LVL_FEATURES feature_id );

#endif /* __INC_SEG_COMMON_H__ */

