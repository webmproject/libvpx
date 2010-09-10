/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef _bigend_h
#define _bigend_h

#if defined(__cplusplus)
extern "C" {
#endif

#define invert2(x) ( (((x)>>8)&0x00ff) | (((x)<<8)&0xff00) )
#define invert4(x) ( ((invert2(x)&0x0000ffff)<<16) | (invert2((x>>16))&0x0000ffff) )

#define high_byte(x) (unsigned char)x
#define mid2Byte(x) (unsigned char)(x >> 8)
#define mid1Byte(x) (unsigned char)(x >> 16)
#define low_byte(x) (unsigned char)(x >> 24)

#define SWAPENDS 1

#if defined(__cplusplus)
}
#endif
#endif
