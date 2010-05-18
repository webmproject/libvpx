/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


/****************************************************************************
*
*   Module Title :     xprintf.h
*
*   Description  :     Debug print interface header file.
*
****************************************************************************/
#ifndef __INC_XPRINTF_H
#define __INC_XPRINTF_H

/****************************************************************************
*  Header Files
****************************************************************************/

/****************************************************************************
*  Functions
****************************************************************************/

// Display a printf style message on the current video frame
extern int onyx_xprintf(unsigned char *ppbuffer, long n_pixel, long n_size, long n_stride, const char *format, ...);

#endif
