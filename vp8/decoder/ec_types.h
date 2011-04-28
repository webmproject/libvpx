/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP8_DEC_EC_TYPES_H
#define VP8_DEC_EC_TYPES_H

#define MAX_OVERLAPS 16

typedef struct
{
    int overlap;
    B_MODE_INFO *bmi;
    MV_REFERENCE_FRAME ref_frame;
} OVERLAP_NODE;

typedef struct
{
    /* TODO(holmer): This array should be exchanged for a linked list */
    OVERLAP_NODE overlaps[MAX_OVERLAPS];
} B_OVERLAP;

typedef struct
{
    B_OVERLAP overlaps[16];
} MB_OVERLAP;

typedef struct
{
    MV mv;
    MV_REFERENCE_FRAME ref_frame;
} EC_BLOCK;

#endif /* VP8_DEC_EC_TYPES_H */
