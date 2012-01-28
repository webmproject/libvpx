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
#include "onyxc_int.h"
#include "vp8/common/blockd.h"

#ifndef __INC_PRED_COMMON_H__
#define __INC_PRED_COMMON_H__ 1


// Predicted items
typedef enum
{
    PRED_SEG_ID = 0,               // Segment identifier

#if CONFIG_COMPRED
    PRED_REF = 1
#endif

} PRED_ID;


extern unsigned char get_pred_context( VP8_COMMON *const cm,
                                       MACROBLOCKD *const xd,
                                       PRED_ID pred_id );

extern vp8_prob get_pred_prob( VP8_COMMON *const cm,
                               MACROBLOCKD *const xd,
                               PRED_ID pred_id );

extern unsigned char get_pred_flag( MACROBLOCKD *const xd,
                                    PRED_ID pred_id );

extern void set_pred_flag( MACROBLOCKD *const xd,
                           PRED_ID pred_id,
                           unsigned char pred_flag);


extern unsigned char get_pred_mb_segid( VP8_COMMON *const cm, int MbIndex );

#if CONFIG_COMPRED
extern MV_REFERENCE_FRAME get_pred_ref( VP8_COMMON *const cm,
                                        MACROBLOCKD *const xd );
extern void compute_mod_refprobs( VP8_COMMON *const cm );

#endif

#endif /* __INC_PRED_COMMON_H__ */
