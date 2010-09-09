/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


/****************************************************************************
*
*   Module Title :     boolhuff.h
*
*   Description  :     Bool Coder header file.
*
****************************************************************************/
#ifndef __INC_BOOLHUFF_H
#define __INC_BOOLHUFF_H


typedef struct
{
    unsigned int lowvalue;
    unsigned int range;
    unsigned int value;
    int count;
    unsigned int pos;
    unsigned char *buffer;

    // Variables used to track bit costs without outputing to the bitstream
    unsigned int  measure_cost;
    unsigned long bit_counter;
} BOOL_CODER;

extern void vp8_start_encode(BOOL_CODER *bc, unsigned char *buffer);
extern void vp8_encode_bool(BOOL_CODER *bc, int x, int context);
extern void vp8_encode_value(BOOL_CODER *br, int data, int bits);
extern void vp8_stop_encode(BOOL_CODER *bc);
extern const unsigned int vp8_prob_cost[256];

#endif
