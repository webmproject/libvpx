/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef _littlend_h
#define _littlend_h

#if defined(__cplusplus)
extern "C" {
#endif

#define invert2(x) (x)
#define invert4(x) (x)

#define low_byte(x) (unsigned char)x
#define mid1Byte(x) (unsigned char)(x >> 8)
#define mid2Byte(x) (unsigned char)(x >> 16)
#define high_byte(x) (unsigned char)(x >> 24)

#define SWAPENDS 0

#if defined(__cplusplus)
}
#endif

#endif
