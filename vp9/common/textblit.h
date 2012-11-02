/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef __INC_TEXTBLIT_H
#define __INC_TEXTBLIT_H

extern void vp9_blit_text(const char *msg, unsigned char *address,
                          const int pitch);
extern void vp9_blit_line(int x0, int x1, int y0, int y1,
                          unsigned char *image, const int pitch);

#endif  // __INC_TEXTBLIT_H
