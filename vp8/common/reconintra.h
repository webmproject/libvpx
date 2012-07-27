/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_RECONINTRA_H
#define __INC_RECONINTRA_H

#include "blockd.h"
void d45_predictor(unsigned char *ypred_ptr, int y_stride, int n,
                   unsigned char *yabove_row, unsigned char *yleft_col);
void d135_predictor(unsigned char *ypred_ptr, int y_stride, int n,
                    unsigned char *yabove_row, unsigned char *yleft_col);
void d116_predictor(unsigned char *ypred_ptr, int y_stride, int n,
                    unsigned char *yabove_row, unsigned char *yleft_col);
void d153_predictor(unsigned char *ypred_ptr, int y_stride, int n,
                    unsigned char *yabove_row, unsigned char *yleft_col);
void d27_predictor(unsigned char *ypred_ptr, int y_stride, int n,
                   unsigned char *yabove_row, unsigned char *yleft_col);
void d64_predictor(unsigned char *ypred_ptr, int y_stride, int n,
                   unsigned char *yabove_row, unsigned char *yleft_col);

extern void init_intra_left_above_pixels(MACROBLOCKD *x);

#endif
