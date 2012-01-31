/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "type_aliases.h"
#include "onyxc_int.h"
#include "vp8/common/blockd.h"

#ifndef __INC_SEG_COMMON_H__
#define __INC_SEG_COMMON_H__ 1

int segfeature_active( MACROBLOCKD *xd,
                       int segment_id,
                       SEG_LVL_FEATURES feature_id );

void clearall_segfeatures( MACROBLOCKD *xd );

void enable_segfeature( MACROBLOCKD *xd,
                        int segment_id,
                        SEG_LVL_FEATURES feature_id );

void disable_segfeature( MACROBLOCKD *xd,
                         int segment_id,
                         SEG_LVL_FEATURES feature_id );

int seg_feature_data_bits( SEG_LVL_FEATURES feature_id );

int is_segfeature_signed( SEG_LVL_FEATURES feature_id );

void clear_segdata( MACROBLOCKD *xd,
                    int segment_id,
                    SEG_LVL_FEATURES feature_id);

void set_segdata( MACROBLOCKD *xd,
                  int segment_id,
                  SEG_LVL_FEATURES feature_id,
                  int seg_data );

int get_segdata( MACROBLOCKD *xd,
                 int segment_id,
                 SEG_LVL_FEATURES feature_id );

void clear_segref( MACROBLOCKD *xd, int segment_id );

void set_segref( MACROBLOCKD *xd,
                 int segment_id,
                 MV_REFERENCE_FRAME ref_frame );

int check_segref( MACROBLOCKD *xd,
                  int segment_id,
                  MV_REFERENCE_FRAME ref_frame );

int check_segref_inter(MACROBLOCKD *xd, int segment_id);

int get_seg_tx_type(MACROBLOCKD *xd, int segment_id);

#endif /* __INC_SEG_COMMON_H__ */

