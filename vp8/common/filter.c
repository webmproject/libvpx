/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <stdlib.h>
#include "filter.h"
#include "vpx_ports/mem.h"

DECLARE_ALIGNED(16, const short, vp8_bilinear_filters[8][2]) =
{
    { 128,   0 },
    { 112,  16 },
    {  96,  32 },
    {  80,  48 },
    {  64,  64 },
    {  48,  80 },
    {  32,  96 },
    {  16, 112 }
};

#if CONFIG_ENHANCED_INTERP
#define FILTER_ALPHA 875
DECLARE_ALIGNED(16, const short, vp8_sub_pel_filters[8][2*INTERP_EXTEND]) =
{
    /* Generated using MATLAB:
     * alpha = 0.875;
     * b=intfilt(8,4,alpha);
     * bi=round(128*b);
     * ba=flipud(reshape([bi 0], 8, 8));
     * % Now normalize the powers of the polyphase components
     * disp(num2str(ba, '%d,'))
     */
#if FILTER_ALPHA == 90
    /* alpha = 0.9 */
    { 0,   0,   0, 128,   0,   0,   0,  0},
    {-3,   6, -13, 126,  18,  -8,   5, -3},
    {-6,  11, -22, 118,  39, -16,   9, -5},
    {-8,  14, -27, 104,  62, -23,  13, -7},
    {-8,  14, -27,  85,  85, -27,  14, -8},
    {-7,  13, -23,  62, 104, -27,  14, -8},
    {-5,   9, -16,  39, 118, -22,  11, -6},
    {-3,   5,  -8,  18, 126, -13,   6, -3}
#elif FILTER_ALPHA == 875
    /* alpha = 0.875 */
    { 0,   0,   0, 128,   0,   0,   0,  0},
    {-3,   6, -13, 126,  18,  -8,   4, -2},
    {-5,  10, -22, 118,  39, -16,   9, -5},
    {-7,  13, -26, 104,  62, -23,  12, -7},
    {-7,  13, -27,  85,  85, -27,  13, -7},
    {-7,  12, -23,  62, 104, -26,  13, -7},
    {-5,   9, -16,  39, 118, -22,  10, -5},
    {-2,   4,  -8,  18, 126, -13,   6, -3}
#elif FILTER_ALPHA == 85
    /* alpha = 0.85 */
    { 0,   0,   0, 128,   0,   0,   0,  0},
    {-2,   5, -12, 124,  18,  -7,   4, -2},
    {-4,  10, -20, 114,  39, -15,   8, -4},
    {-5,  12, -24, 100,  60, -21,  11, -5},
    {-5,  12, -24,  81,  81, -24,  12, -5},
    {-5,  11, -21,  60, 100, -24,  12, -5},
    {-4,   8, -15,  39, 114, -20,  10, -4},
    {-2,   4,  -7,  18, 124, -12,   5, -2}
#endif
};

#if EDGE_PIXEL_FILTER > 0

#define EDGE_SIMPLE_THRESH 128
#define EDGE_GRAD_THRESH 128
#define EDGE_GRADS2X2_THRESH 4
/* TODO: Refine these filters */
DECLARE_ALIGNED(16, const short, vp8_sub_pel_filters_ns[64][4*EDGE_PIXEL_FILTER_EXTEND*EDGE_PIXEL_FILTER_EXTEND]) =
{
#if EDGE_PIXEL_FILTER_EXTEND == 4
    {0,   0,   0,   0,   0, 128,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   0,   0,  -8, 124,  13,  -1,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   0,   0, -11, 112,  30,  -3,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   0,   0, -11,  94,  51,  -6,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,  0,  0,  0, -9, 73, 73, -9,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,   0,   0,  -6,  51,  94, -11,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   0,   0,  -3,  30, 112, -11,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   0,   0,  -1,  13, 124,  -8,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,  -8,   0,   0,   0, 124,   0,   0,   0,  13,   0,   0,   0,  -1,   0, 0},
    {0,  -8,  -1,   0,  -7, 120,  12,   0,  -1,  12,   1,   0,   0,   0,   0, 0},
    {1,  -7,  -2,   0, -11, 108,  29,  -3,  -1,  11,   3,   0,   0,   0,   0, 0},
    {0,  -6,  -3,   0, -10,  91,  49,  -6,  -1,   9,   5,   0,   0,   0,   0, 0},
    {0, -4, -4,  0, -8, 70, 70, -8, -1,  7,  7, -1,  0,  0,  0, 0},
    {0,  -3,  -6,   0,  -6,  49,  91, -10,   0,   5,   9,  -1,   0,   0,   0, 0},
    {0,  -2,  -7,   1,  -3,  29, 108, -11,   0,   3,  11,  -1,   0,   0,   0, 0},
    {0,  -1,  -8,   0,   0,  12, 120,  -7,   0,   1,  12,  -1,   0,   0,   0, 0},
    {0, -11,   0,   0,   0, 112,   0,   0,   0,  30,   0,   0,   0,  -3,   0, 0},
    {1, -11,  -1,   0,  -7, 108,  11,   0,  -2,  29,   3,   0,   0,  -3,   0, 0},
    {1, -10,  -2,   0, -10,  98,  26,  -3,  -2,  26,   7,   0,   0,  -3,   0, 0},
    {1,  -8,  -4,   0, -10,  82,  44,  -5,  -2,  22,  12,  -1,   0,  -2,  -1, 0},
    {0, -6, -6,  0, -7, 64, 64, -7, -2, 17, 17, -2,  0, -2, -2, 0},
    {0,  -4,  -8,   1,  -5,  44,  82, -10,  -1,  12,  22,  -2,   0,  -1,  -2, 0},
    {0,  -2, -10,   1,  -3,  26,  98, -10,   0,   7,  26,  -2,   0,   0,  -3, 0},
    {0,  -1, -11,   1,   0,  11, 108,  -7,   0,   3,  29,  -2,   0,   0,  -3, 0},
    {0, -11,   0,   0,   0,  94,   0,   0,   0,  51,   0,   0,   0,  -6,   0, 0},
    {0, -10,  -1,   0,  -6,  91,   9,   0,  -3,  49,   5,   0,   0,  -6,   0, 0},
    {1, -10,  -2,   0,  -8,  82,  22,  -2,  -4,  44,  12,  -1,   0,  -5,  -1, 0},
    {0, -8, -4,  0, -8, 70, 37, -4, -4, 37, 20, -2,  0, -4, -2, 0},
    {0, -6, -6,  0, -6, 54, 54, -6, -3, 28, 28, -3,  0, -3, -3, 0},
    {0, -4, -8,  0, -4, 37, 70, -8, -2, 20, 37, -4,  0, -2, -4, 0},
    {0,  -2, -10,   1,  -2,  22,  82,  -8,  -1,  12,  44,  -4,   0,  -1,  -5, 0},
    {0,  -1, -10,   0,   0,   9,  91,  -6,   0,   5,  49,  -3,   0,   0,  -6, 0},
    {0, -9,  0,  0,  0, 73,  0,  0,  0, 73,  0,  0,  0, -9,  0, 0},
    {0, -8, -1,  0, -4, 70,  7,  0, -4, 70,  7,  0,  0, -8, -1, 0},
    {0, -7, -2,  0, -6, 64, 17, -2, -6, 64, 17, -2,  0, -7, -2, 0},
    {0, -6, -3,  0, -6, 54, 28, -3, -6, 54, 28, -3,  0, -6, -3, 0},
    {0, -5, -5,  0, -5, 42, 42, -5, -5, 42, 42, -5,  0, -5, -5, 0},
    {0, -3, -6,  0, -3, 28, 54, -6, -3, 28, 54, -6,  0, -3, -6, 0},
    {0, -2, -7,  0, -2, 17, 64, -6, -2, 17, 64, -6,  0, -2, -7, 0},
    {0, -1, -8,  0,  0,  7, 70, -4,  0,  7, 70, -4,  0, -1, -8, 0},
    {0,  -6,   0,   0,   0,  51,   0,   0,   0,  94,   0,   0,   0, -11,   0, 0},
    {0,  -6,   0,   0,  -3,  49,   5,   0,  -6,  91,   9,   0,   0, -10,  -1, 0},
    {0,  -5,  -1,   0,  -4,  44,  12,  -1,  -8,  82,  22,  -2,   1, -10,  -2, 0},
    {0, -4, -2,  0, -4, 37, 20, -2, -8, 70, 37, -4,  0, -8, -4, 0},
    {0, -3, -3,  0, -3, 28, 28, -3, -6, 54, 54, -6,  0, -6, -6, 0},
    {0, -2, -4,  0, -2, 20, 37, -4, -4, 37, 70, -8,  0, -4, -8, 0},
    {0,  -1,  -5,   0,  -1,  12,  44,  -4,  -2,  22,  82,  -8,   0,  -2, -10, 1},
    {0,   0,  -6,   0,   0,   5,  49,  -3,   0,   9,  91,  -6,   0,  -1, -10, 0},
    {0,  -3,   0,   0,   0,  30,   0,   0,   0, 112,   0,   0,   0, -11,   0, 0},
    {0,  -3,   0,   0,  -2,  29,   3,   0,  -7, 108,  11,   0,   1, -11,  -1, 0},
    {0,  -3,   0,   0,  -2,  26,   7,   0, -10,  98,  26,  -3,   1, -10,  -2, 0},
    {0,  -2,  -1,   0,  -2,  22,  12,  -1, -10,  82,  44,  -5,   1,  -8,  -4, 0},
    {0, -2, -2,  0, -2, 17, 17, -2, -7, 64, 64, -7,  0, -6, -6, 0},
    {0,  -1,  -2,   0,  -1,  12,  22,  -2,  -5,  44,  82, -10,   0,  -4,  -8, 1},
    {0,   0,  -3,   0,   0,   7,  26,  -2,  -3,  26,  98, -10,   0,  -2, -10, 1},
    {0,   0,  -3,   0,   0,   3,  29,  -2,   0,  11, 108,  -7,   0,  -1, -11, 1},
    {0,  -1,   0,   0,   0,  13,   0,   0,   0, 124,   0,   0,   0,  -8,   0, 0},
    {0,   0,   0,   0,  -1,  12,   1,   0,  -8, 120,  12,   0,   0,  -7,  -1, 0},
    {0,   0,   0,   0,  -1,  11,   3,   0, -11, 108,  29,  -3,   1,  -7,  -2, 0},
    {0,   0,   0,   0,  -1,   9,   5,   0, -10,  91,  49,  -6,   0,  -6,  -3, 0},
    {0,  0,  0,  0, -1,  7,  7, -1, -8, 70, 70, -8,  0, -4, -4, 0},
    {0,   0,   0,   0,   0,   5,   9,  -1,  -6,  49,  91, -10,   0,  -3,  -6, 0},
    {0,   0,   0,   0,   0,   3,  11,  -1,  -3,  29, 108, -11,   0,  -2,  -7, 1},
    {0,   0,   0,   0,   0,   1,  12,  -1,   0,  12, 120,  -8,   0,  -1,  -7, 0}
#elif EDGE_PIXEL_FILTER_EXTEND == 6
    {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 128,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3, -11, 124,  15,  -4,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   4, -17, 114,  34,  -9,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   4, -19,  98,  56, -14,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   4, -18,  78,  78, -18,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3, -14,  56,  98, -19,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,  -9,  34, 114, -17,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  -4,  15, 124, -11,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   3,   0,   0,   0,   0,   0, -11,   0,   0,   0,   0,   0, 124,   0,   0,   0,   0,   0,  15,   0,   0,   0,   0,   0,  -4,   0,   0,   0,   0,   0,   1,   0,   0, 0},
    {0,   0,   3,   0,   0,   0,   0,   1, -11,  -1,   0,   0,   3, -11, 121,  15,  -4,   0,   0,  -1,  15,   2,   0,   0,   0,   0,  -4,   0,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   3,   1,   0,   0,   0,   1, -10,  -3,   1,   0,   4, -17, 111,  34,  -9,   2,   0,  -2,  14,   4,  -1,   0,   0,   0,  -4,  -1,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   2,   1,   0,   0,   0,   2,  -8,  -5,   1,   0,   4, -18,  95,  54, -13,   3,   0,  -2,  12,   7,  -2,   0,   0,   0,  -3,  -2,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   2,   2,   0,   0,   0,   1,  -7,  -7,   1,   0,   4, -17,  76,  76, -17,   4,   0,  -2,   9,   9,  -2,   0,   0,   0,  -2,  -2,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   1,   2,   0,   0,   0,   1,  -5,  -8,   2,   0,   3, -13,  54,  95, -18,   4,   0,  -2,   7,  12,  -2,   0,   0,   0,  -2,  -3,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   1,   3,   0,   0,   0,   1,  -3, -10,   1,   0,   2,  -9,  34, 111, -17,   4,   0,  -1,   4,  14,  -2,   0,   0,   0,  -1,  -4,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   0,   3,   0,   0,   0,   0,  -1, -11,   1,   0,   0,  -4,  15, 121, -11,   3,   0,   0,   2,  15,  -1,   0,   0,   0,   0,  -4,   0,   0,   0,   0,   0,   0,   0, 0},
    {0,   0,   4,   0,   0,   0,   0,   0, -17,   0,   0,   0,   0,   0, 114,   0,   0,   0,   0,   0,  34,   0,   0,   0,   0,   0,  -9,   0,   0,   0,   0,   0,   2,   0,   0, 0},
    {0,   0,   4,   0,   0,   0,   0,   1, -17,  -2,   0,   0,   3, -10, 111,  14,  -4,   0,   1,  -3,  34,   4,  -1,   0,   0,   1,  -9,  -1,   0,   0,   0,   0,   2,   0,   0, 0},
    {0,   0,   4,   1,   0,   0,  -1,   2, -15,  -5,   1,   0,   4, -15, 101,  31,  -8,   1,   1,  -5,  31,   9,  -2,   0,   0,   1,  -8,  -2,   1,   0,   0,   0,   1,   0,   0, 0},
    {0,  -1,   3,   2,   0,   0,  -1,   2, -13,  -7,   2,   0,   4, -17,  87,  50, -12,   2,   1,  -5,  26,  15,  -4,   1,   0,   1,  -7,  -4,   1,   0,   0,   0,   1,   1,   0, 0},
    {0,  -1,   3,   3,  -1,   0,   0,   2, -10, -10,   2,   0,   3, -16,  69,  69, -16,   3,   1,  -5,  21,  21,  -5,   1,   0,   1,  -5,  -5,   1,   0,   0,   0,   1,   1,   0, 0},
    {0,   0,   2,   3,  -1,   0,   0,   2,  -7, -13,   2,  -1,   2, -12,  50,  87, -17,   4,   1,  -4,  15,  26,  -5,   1,   0,   1,  -4,  -7,   1,   0,   0,   0,   1,   1,   0, 0},
    {0,   0,   1,   4,   0,   0,   0,   1,  -5, -15,   2,  -1,   1,  -8,  31, 101, -15,   4,   0,  -2,   9,  31,  -5,   1,   0,   1,  -2,  -8,   1,   0,   0,   0,   0,   1,   0, 0},
    {0,   0,   0,   4,   0,   0,   0,   0,  -2, -17,   1,   0,   0,  -4,  14, 111, -10,   3,   0,  -1,   4,  34,  -3,   1,   0,   0,  -1,  -9,   1,   0,   0,   0,   0,   2,   0, 0},
    {0,   0,   4,   0,   0,   0,   0,   0, -19,   0,   0,   0,   0,   0,  98,   0,   0,   0,   0,   0,  56,   0,   0,   0,   0,   0, -14,   0,   0,   0,   0,   0,   3,   0,   0, 0},
    {0,   0,   4,   0,   0,   0,   0,   2, -18,  -2,   0,   0,   2,  -8,  95,  12,  -3,   0,   1,  -5,  54,   7,  -2,   0,   0,   1, -13,  -2,   0,   0,   0,   0,   3,   0,   0, 0},
    {0,  -1,   4,   1,   0,   0,  -1,   2, -17,  -5,   1,   0,   3, -13,  87,  26,  -7,   1,   2,  -7,  50,  15,  -4,   1,   0,   2, -12,  -4,   1,   0,   0,   0,   2,   1,   0, 0},
    {0,  -1,   3,   2,   0,   0,  -1,   3, -15,  -8,   2,   0,   3, -15,  75,  43, -11,   2,   2,  -8,  43,  24,  -6,   1,   0,   2, -10,  -6,   1,   0,   0,   0,   2,   1,   0, 0},
    {0,  -1,   3,   3,  -1,   0,  -1,   2, -12, -12,   2,  -1,   3, -13,  59,  59, -13,   3,   2,  -8,  34,  34,  -8,   2,   0,   2,  -8,  -8,   2,   0,   0,   0,   2,   2,   0, 0},
    {0,   0,   2,   3,  -1,   0,   0,   2,  -8, -15,   3,  -1,   2, -11,  43,  75, -15,   3,   1,  -6,  24,  43,  -8,   2,   0,   1,  -6, -10,   2,   0,   0,   0,   1,   2,   0, 0},
    {0,   0,   1,   4,  -1,   0,   0,   1,  -5, -17,   2,  -1,   1,  -7,  26,  87, -13,   3,   1,  -4,  15,  50,  -7,   2,   0,   1,  -4, -12,   2,   0,   0,   0,   1,   2,   0, 0},
    {0,   0,   0,   4,   0,   0,   0,   0,  -2, -18,   2,   0,   0,  -3,  12,  95,  -8,   2,   0,  -2,   7,  54,  -5,   1,   0,   0,  -2, -13,   1,   0,   0,   0,   0,   3,   0, 0},
    {0,   0,   4,   0,   0,   0,   0,   0, -18,   0,   0,   0,   0,   0,  78,   0,   0,   0,   0,   0,  78,   0,   0,   0,   0,   0, -18,   0,   0,   0,   0,   0,   4,   0,   0, 0},
    {0,   0,   4,   0,   0,   0,   0,   1, -17,  -2,   0,   0,   2,  -7,  76,   9,  -2,   0,   2,  -7,  76,   9,  -2,   0,   0,   1, -17,  -2,   0,   0,   0,   0,   4,   0,   0, 0},
    {0,   0,   3,   1,   0,   0,  -1,   2, -16,  -5,   1,   0,   3, -10,  69,  21,  -5,   1,   3, -10,  69,  21,  -5,   1,  -1,   2, -16,  -5,   1,   0,   0,   0,   3,   1,   0, 0},
    {0,  -1,   3,   2,   0,   0,  -1,   2, -13,  -8,   2,   0,   3, -12,  59,  34,  -8,   2,   3, -12,  59,  34,  -8,   2,  -1,   2, -13,  -8,   2,   0,   0,  -1,   3,   2,   0, 0},
    {0,   0,   2,   2,   0,   0,   0,   2, -10, -10,   2,   0,   2, -10,  47,  47, -10,   2,   2, -11,  47,  47, -11,   2,   0,   2, -11, -11,   2,   0,   0,   0,   2,   2,   0, 0},
    {0,   0,   2,   3,  -1,   0,   0,   2,  -8, -13,   2,  -1,   2,  -8,  34,  59, -12,   3,   2,  -8,  34,  59, -12,   3,   0,   2,  -8, -13,   2,  -1,   0,   0,   2,   3,  -1, 0},
    {0,   0,   1,   3,   0,   0,   0,   1,  -5, -16,   2,  -1,   1,  -5,  21,  69, -10,   3,   1,  -5,  21,  69, -10,   3,   0,   1,  -5, -16,   2,  -1,   0,   0,   1,   3,   0, 0},
    {0,   0,   0,   4,   0,   0,   0,   0,  -2, -17,   1,   0,   0,  -2,   9,  76,  -7,   2,   0,  -2,   9,  76,  -7,   2,   0,   0,  -2, -17,   1,   0,   0,   0,   0,   4,   0, 0},
    {0,   0,   3,   0,   0,   0,   0,   0, -14,   0,   0,   0,   0,   0,  56,   0,   0,   0,   0,   0,  98,   0,   0,   0,   0,   0, -19,   0,   0,   0,   0,   0,   4,   0,   0, 0},
    {0,   0,   3,   0,   0,   0,   0,   1, -13,  -2,   0,   0,   1,  -5,  54,   7,  -2,   0,   2,  -8,  95,  12,  -3,   0,   0,   2, -18,  -2,   0,   0,   0,   0,   4,   0,   0, 0},
    {0,   0,   2,   1,   0,   0,   0,   2, -12,  -4,   1,   0,   2,  -7,  50,  15,  -4,   1,   3, -13,  87,  26,  -7,   1,  -1,   2, -17,  -5,   1,   0,   0,  -1,   4,   1,   0, 0},
    {0,   0,   2,   1,   0,   0,   0,   2, -11,  -6,   1,   0,   2,  -8,  43,  24,  -6,   1,   3, -15,  75,  43, -10,   2,  -1,   3, -15,  -8,   2,   0,   0,  -1,   3,   2,   0, 0},
    {0,   0,   2,   2,   0,   0,   0,   2,  -8,  -8,   2,   0,   2,  -8,  34,  34,  -8,   2,   3, -13,  59,  59, -13,   3,  -1,   2, -12, -12,   2,  -1,   0,  -1,   3,   3,  -1, 0},
    {0,   0,   1,   2,   0,   0,   0,   1,  -6, -11,   2,   0,   1,  -6,  24,  43,  -8,   2,   2, -10,  43,  75, -15,   3,   0,   2,  -8, -15,   3,  -1,   0,   0,   2,   3,  -1, 0},
    {0,   0,   1,   2,   0,   0,   0,   1,  -4, -12,   2,   0,   1,  -4,  15,  50,  -7,   2,   1,  -7,  26,  87, -13,   3,   0,   1,  -5, -17,   2,  -1,   0,   0,   1,   4,  -1, 0},
    {0,   0,   0,   3,   0,   0,   0,   0,  -2, -13,   1,   0,   0,  -2,   7,  54,  -5,   1,   0,  -3,  12,  95,  -8,   2,   0,   0,  -2, -18,   2,   0,   0,   0,   0,   4,   0, 0},
    {0,   0,   2,   0,   0,   0,   0,   0,  -9,   0,   0,   0,   0,   0,  34,   0,   0,   0,   0,   0, 114,   0,   0,   0,   0,   0, -17,   0,   0,   0,   0,   0,   4,   0,   0, 0},
    {0,   0,   2,   0,   0,   0,   0,   1,  -9,  -1,   0,   0,   1,  -3,  34,   4,  -1,   0,   3, -10, 111,  14,  -4,   0,   0,   1, -17,  -2,   0,   0,   0,   0,   4,   0,   0, 0},
    {0,   0,   1,   0,   0,   0,   0,   1,  -8,  -2,   1,   0,   1,  -5,  31,   9,  -2,   0,   4, -15, 101,  31,  -8,   1,   0,   2, -15,  -5,   1,   0,   0,  -1,   4,   1,   0, 0},
    {0,   0,   1,   1,   0,   0,   0,   1,  -7,  -4,   1,   0,   1,  -5,  26,  15,  -4,   1,   4, -17,  87,  50, -12,   2,  -1,   2, -13,  -7,   2,   0,   0,  -1,   3,   2,   0, 0},
    {0,   0,   1,   1,   0,   0,   0,   1,  -5,  -5,   1,   0,   1,  -5,  21,  21,  -5,   1,   3, -16,  69,  69, -16,   3,   0,   2, -10, -10,   2,   0,   0,  -1,   3,   3,  -1, 0},
    {0,   0,   1,   1,   0,   0,   0,   1,  -4,  -7,   1,   0,   1,  -4,  15,  26,  -5,   1,   2, -12,  50,  87, -17,   4,   0,   2,  -7, -13,   2,  -1,   0,   0,   2,   3,  -1, 0},
    {0,   0,   0,   1,   0,   0,   0,   1,  -2,  -8,   1,   0,   0,  -2,   9,  31,  -5,   1,   1,  -8,  31, 101, -15,   4,   0,   1,  -5, -15,   2,   0,   0,   0,   1,   4,  -1, 0},
    {0,   0,   0,   2,   0,   0,   0,   0,  -1,  -9,   1,   0,   0,  -1,   4,  34,  -3,   1,   0,  -4,  14, 111, -10,   3,   0,   0,  -2, -17,   1,   0,   0,   0,   0,   4,   0, 0},
    {0,   0,   1,   0,   0,   0,   0,   0,  -4,   0,   0,   0,   0,   0,  15,   0,   0,   0,   0,   0, 124,   0,   0,   0,   0,   0, -11,   0,   0,   0,   0,   0,   3,   0,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,  -4,   0,   0,   0,   0,  -1,  15,   2,   0,   0,   3, -11, 121,  15,  -4,   0,   0,   1, -11,  -1,   0,   0,   0,   0,   3,   0,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,  -4,  -1,   0,   0,   0,  -2,  14,   4,  -1,   0,   4, -17, 111,  34,  -9,   2,   0,   1, -10,  -3,   1,   0,   0,   0,   3,   1,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,  -3,  -2,   0,   0,   0,  -2,  12,   7,  -2,   0,   4, -18,  95,  54, -13,   3,   0,   2,  -8,  -5,   1,   0,   0,   0,   2,   1,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,  -2,  -2,   0,   0,   0,  -2,   9,   9,  -2,   0,   4, -17,  76,  76, -17,   4,   0,   1,  -7,  -7,   1,   0,   0,   0,   2,   2,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,  -2,  -3,   0,   0,   0,  -2,   7,  12,  -2,   0,   3, -13,  54,  95, -18,   4,   0,   1,  -5,  -8,   2,   0,   0,   0,   1,   2,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,  -1,  -4,   0,   0,   0,  -1,   4,  14,  -2,   0,   2,  -9,  34, 111, -17,   4,   0,   1,  -3, -10,   1,   0,   0,   0,   1,   3,   0, 0},
    {0,   0,   0,   0,   0,   0,   0,   0,   0,  -4,   0,   0,   0,   0,   2,  15,  -1,   0,   0,  -4,  15, 121, -11,   3,   0,   0,  -1, -11,   1,   0,   0,   0,   0,   3,   0, 0}
#endif
};

#endif  // EDGE_PIXEL_FILTER

#else  // CONFIG_ENHANCED_INTERP

DECLARE_ALIGNED(16, const short, vp8_sub_pel_filters[8][6]) =
{

    { 0,  0,  128,    0,   0,  0 },         /* note that 1/8 pel positions are just as per alpha -0.5 bicubic */
    { 0, -6,  123,   12,  -1,  0 },
    { 2, -11, 108,   36,  -8,  1 },         /* New 1/4 pel 6 tap filter */
    { 0, -9,   93,   50,  -6,  0 },
    { 3, -16,  77,   77, -16,  3 },         /* New 1/2 pel 6 tap filter */
    { 0, -6,   50,   93,  -9,  0 },
    { 1, -8,   36,  108, -11,  2 },         /* New 1/4 pel 6 tap filter */
    { 0, -1,   12,  123,  -6,  0 },
};
#endif  // CONFIG_ENHANCED_INTERP

static void filter_block2d_first_pass
(
    unsigned char *src_ptr,
    int *output_ptr,
    unsigned int src_pixels_per_line,
    unsigned int pixel_step,
    unsigned int output_height,
    unsigned int output_width,
    const short *vp8_filter
)
{
    unsigned int i, j;
    int  Temp;

    for (i = 0; i < output_height; i++)
    {
        for (j = 0; j < output_width; j++)
        {
#if INTERP_EXTEND == 3
            Temp = ((int)src_ptr[-2 * (int)pixel_step] * vp8_filter[0]) +
                   ((int)src_ptr[-1 * (int)pixel_step] * vp8_filter[1]) +
                   ((int)src_ptr[0]                    * vp8_filter[2]) +
                   ((int)src_ptr[pixel_step]           * vp8_filter[3]) +
                   ((int)src_ptr[2*pixel_step]         * vp8_filter[4]) +
                   ((int)src_ptr[3*pixel_step]         * vp8_filter[5]) +
                   (VP8_FILTER_WEIGHT >> 1);      /* Rounding */
#elif INTERP_EXTEND == 4
            Temp = ((int)src_ptr[-3 * (int)pixel_step] * vp8_filter[0]) +
                   ((int)src_ptr[-2 * (int)pixel_step] * vp8_filter[1]) +
                   ((int)src_ptr[-1 * (int)pixel_step] * vp8_filter[2]) +
                   ((int)src_ptr[0]                    * vp8_filter[3]) +
                   ((int)src_ptr[pixel_step]           * vp8_filter[4]) +
                   ((int)src_ptr[2 * pixel_step]       * vp8_filter[5]) +
                   ((int)src_ptr[3 * pixel_step]       * vp8_filter[6]) +
                   ((int)src_ptr[4 * pixel_step]       * vp8_filter[7]) +
                   (VP8_FILTER_WEIGHT >> 1);      /* Rounding */
#elif INTERP_EXTEND == 5
            Temp = ((int)src_ptr[-4 * (int)pixel_step] * vp8_filter[0]) +
                   ((int)src_ptr[-3 * (int)pixel_step] * vp8_filter[1]) +
                   ((int)src_ptr[-2 * (int)pixel_step] * vp8_filter[2]) +
                   ((int)src_ptr[-1 * (int)pixel_step] * vp8_filter[3]) +
                   ((int)src_ptr[0]                    * vp8_filter[4]) +
                   ((int)src_ptr[pixel_step]           * vp8_filter[5]) +
                   ((int)src_ptr[2 * pixel_step]       * vp8_filter[6]) +
                   ((int)src_ptr[3 * pixel_step]       * vp8_filter[7]) +
                   ((int)src_ptr[4 * pixel_step]       * vp8_filter[8]) +
                   ((int)src_ptr[5 * pixel_step]       * vp8_filter[9]) +
                   (VP8_FILTER_WEIGHT >> 1);      /* Rounding */
#endif

            /* Normalize back to 0-255 */
            Temp = Temp >> VP8_FILTER_SHIFT;

            if (Temp < 0)
                Temp = 0;
            else if (Temp > 255)
                Temp = 255;

            output_ptr[j] = Temp;
            src_ptr++;
        }

        /* Next row... */
        src_ptr    += src_pixels_per_line - output_width;
        output_ptr += output_width;
    }
}

static void filter_block2d_second_pass
(
    int *src_ptr,
    unsigned char *output_ptr,
    int output_pitch,
    unsigned int src_pixels_per_line,
    unsigned int pixel_step,
    unsigned int output_height,
    unsigned int output_width,
    const short *vp8_filter
)
{
    unsigned int i, j;
    int  Temp;

    for (i = 0; i < output_height; i++)
    {
        for (j = 0; j < output_width; j++)
        {
            /* Apply filter */
#if INTERP_EXTEND == 3
            Temp = ((int)src_ptr[-2 * (int)pixel_step] * vp8_filter[0]) +
                   ((int)src_ptr[-1 * (int)pixel_step] * vp8_filter[1]) +
                   ((int)src_ptr[0]                    * vp8_filter[2]) +
                   ((int)src_ptr[pixel_step]           * vp8_filter[3]) +
                   ((int)src_ptr[2*pixel_step]         * vp8_filter[4]) +
                   ((int)src_ptr[3*pixel_step]         * vp8_filter[5]) +
                   (VP8_FILTER_WEIGHT >> 1);   /* Rounding */
#elif INTERP_EXTEND == 4
            Temp = ((int)src_ptr[-3 * (int)pixel_step] * vp8_filter[0]) +
                   ((int)src_ptr[-2 * (int)pixel_step] * vp8_filter[1]) +
                   ((int)src_ptr[-1 * (int)pixel_step] * vp8_filter[2]) +
                   ((int)src_ptr[0]                    * vp8_filter[3]) +
                   ((int)src_ptr[pixel_step]           * vp8_filter[4]) +
                   ((int)src_ptr[2 * pixel_step]       * vp8_filter[5]) +
                   ((int)src_ptr[3 * pixel_step]       * vp8_filter[6]) +
                   ((int)src_ptr[4 * pixel_step]       * vp8_filter[7]) +
                   (VP8_FILTER_WEIGHT >> 1);      /* Rounding */
#elif INTERP_EXTEND == 5
            Temp = ((int)src_ptr[-4 * (int)pixel_step] * vp8_filter[0]) +
                   ((int)src_ptr[-3 * (int)pixel_step] * vp8_filter[1]) +
                   ((int)src_ptr[-2 * (int)pixel_step] * vp8_filter[2]) +
                   ((int)src_ptr[-1 * (int)pixel_step] * vp8_filter[3]) +
                   ((int)src_ptr[0]                    * vp8_filter[4]) +
                   ((int)src_ptr[pixel_step]           * vp8_filter[5]) +
                   ((int)src_ptr[2 * pixel_step]       * vp8_filter[6]) +
                   ((int)src_ptr[3 * pixel_step]       * vp8_filter[7]) +
                   ((int)src_ptr[4 * pixel_step]       * vp8_filter[8]) +
                   ((int)src_ptr[5 * pixel_step]       * vp8_filter[9]) +
                   (VP8_FILTER_WEIGHT >> 1);      /* Rounding */
#endif

            /* Normalize back to 0-255 */
            Temp = Temp >> VP8_FILTER_SHIFT;

            if (Temp < 0)
                Temp = 0;
            else if (Temp > 255)
                Temp = 255;

            output_ptr[j] = (unsigned char)Temp;
            src_ptr++;
        }

        /* Start next row */
        src_ptr    += src_pixels_per_line - output_width;
        output_ptr += output_pitch;
    }
}

#if EDGE_PIXEL_FILTER > 0
static void filter_non_separable
(

    unsigned char *src_ptr,
    unsigned char *output_ptr,
    unsigned int src_pixels_per_line,
    unsigned int pixel_step,
    const short *vp8_filter
)
{
    int Temp;
#if EDGE_PIXEL_FILTER_EXTEND == 2
    /* This code computes non-separable filtering of a pixel
     * using a 4x4 neighborhood as shown where F is the pixel
     * that src_ptr points to:
     *
     *    A B C D
     *    E F G H
     *    I J K L
     *    M N O P
     *
     * The 16 filter coefficients are in row by row order
     * */
    Temp = ((int)src_ptr[-1 * (int)src_pixels_per_line - 1 * (int)pixel_step] * vp8_filter[0]) +
           ((int)src_ptr[-1 * (int)src_pixels_per_line + 0 * (int)pixel_step] * vp8_filter[1]) +
           ((int)src_ptr[-1 * (int)src_pixels_per_line + 1 * (int)pixel_step] * vp8_filter[2]) +
           ((int)src_ptr[-1 * (int)src_pixels_per_line + 2 * (int)pixel_step] * vp8_filter[3]) +
           ((int)src_ptr[ 0 * (int)src_pixels_per_line - 1 * (int)pixel_step] * vp8_filter[4]) +
           ((int)src_ptr[ 0 * (int)src_pixels_per_line + 0 * (int)pixel_step] * vp8_filter[5]) +
           ((int)src_ptr[ 0 * (int)src_pixels_per_line + 1 * (int)pixel_step] * vp8_filter[6]) +
           ((int)src_ptr[ 0 * (int)src_pixels_per_line + 2 * (int)pixel_step] * vp8_filter[7]) +
           ((int)src_ptr[ 1 * (int)src_pixels_per_line - 1 * (int)pixel_step] * vp8_filter[8]) +
           ((int)src_ptr[ 1 * (int)src_pixels_per_line + 0 * (int)pixel_step] * vp8_filter[9]) +
           ((int)src_ptr[ 1 * (int)src_pixels_per_line + 1 * (int)pixel_step] * vp8_filter[10]) +
           ((int)src_ptr[ 1 * (int)src_pixels_per_line + 2 * (int)pixel_step] * vp8_filter[11]) +
           ((int)src_ptr[ 2 * (int)src_pixels_per_line - 1 * (int)pixel_step] * vp8_filter[12]) +
           ((int)src_ptr[ 2 * (int)src_pixels_per_line + 0 * (int)pixel_step] * vp8_filter[13]) +
           ((int)src_ptr[ 2 * (int)src_pixels_per_line + 1 * (int)pixel_step] * vp8_filter[14]) +
           ((int)src_ptr[ 2 * (int)src_pixels_per_line + 2 * (int)pixel_step] * vp8_filter[15]) +
           (VP8_FILTER_WEIGHT >> 1);      /* Rounding */
#elif EDGE_PIXEL_FILTER_EXTEND == 3
    /* This code computes non-separable filtering of a pixel
     * using a 6x6 neighborhood as shown where O is the pixel
     * that src_ptr points to:
     *
     *    A B C D E F
     *    G H I J K L
     *    M N O P Q R
     *    S T U V W X
     *    Y Z a b c d
     *    e f g h i j
     *
     * The 36 filter coefficients are in row by row order
     * */
    Temp = ((int)src_ptr[-2 * (int)src_pixels_per_line - 2 * (int)pixel_step] * vp8_filter[0]) +
           ((int)src_ptr[-2 * (int)src_pixels_per_line - 1 * (int)pixel_step] * vp8_filter[1]) +
           ((int)src_ptr[-2 * (int)src_pixels_per_line + 0 * (int)pixel_step] * vp8_filter[2]) +
           ((int)src_ptr[-2 * (int)src_pixels_per_line + 1 * (int)pixel_step] * vp8_filter[3]) +
           ((int)src_ptr[-2 * (int)src_pixels_per_line + 2 * (int)pixel_step] * vp8_filter[4]) +
           ((int)src_ptr[-2 * (int)src_pixels_per_line + 3 * (int)pixel_step] * vp8_filter[5]) +
           ((int)src_ptr[-1 * (int)src_pixels_per_line - 2 * (int)pixel_step] * vp8_filter[6]) +
           ((int)src_ptr[-1 * (int)src_pixels_per_line - 1 * (int)pixel_step] * vp8_filter[7]) +
           ((int)src_ptr[-1 * (int)src_pixels_per_line + 0 * (int)pixel_step] * vp8_filter[8]) +
           ((int)src_ptr[-1 * (int)src_pixels_per_line + 1 * (int)pixel_step] * vp8_filter[9]) +
           ((int)src_ptr[-1 * (int)src_pixels_per_line + 2 * (int)pixel_step] * vp8_filter[10]) +
           ((int)src_ptr[-1 * (int)src_pixels_per_line + 3 * (int)pixel_step] * vp8_filter[11]) +
           ((int)src_ptr[ 0 * (int)src_pixels_per_line - 2 * (int)pixel_step] * vp8_filter[12]) +
           ((int)src_ptr[ 0 * (int)src_pixels_per_line - 1 * (int)pixel_step] * vp8_filter[13]) +
           ((int)src_ptr[ 0 * (int)src_pixels_per_line + 0 * (int)pixel_step] * vp8_filter[14]) +
           ((int)src_ptr[ 0 * (int)src_pixels_per_line + 1 * (int)pixel_step] * vp8_filter[15]) +
           ((int)src_ptr[ 0 * (int)src_pixels_per_line + 2 * (int)pixel_step] * vp8_filter[16]) +
           ((int)src_ptr[ 0 * (int)src_pixels_per_line + 3 * (int)pixel_step] * vp8_filter[17]) +
           ((int)src_ptr[ 1 * (int)src_pixels_per_line - 2 * (int)pixel_step] * vp8_filter[18]) +
           ((int)src_ptr[ 1 * (int)src_pixels_per_line - 1 * (int)pixel_step] * vp8_filter[19]) +
           ((int)src_ptr[ 1 * (int)src_pixels_per_line + 0 * (int)pixel_step] * vp8_filter[20]) +
           ((int)src_ptr[ 1 * (int)src_pixels_per_line + 1 * (int)pixel_step] * vp8_filter[21]) +
           ((int)src_ptr[ 1 * (int)src_pixels_per_line + 2 * (int)pixel_step] * vp8_filter[22]) +
           ((int)src_ptr[ 1 * (int)src_pixels_per_line + 3 * (int)pixel_step] * vp8_filter[23]) +
           ((int)src_ptr[ 2 * (int)src_pixels_per_line - 2 * (int)pixel_step] * vp8_filter[24]) +
           ((int)src_ptr[ 2 * (int)src_pixels_per_line - 1 * (int)pixel_step] * vp8_filter[25]) +
           ((int)src_ptr[ 2 * (int)src_pixels_per_line + 0 * (int)pixel_step] * vp8_filter[26]) +
           ((int)src_ptr[ 2 * (int)src_pixels_per_line + 1 * (int)pixel_step] * vp8_filter[27]) +
           ((int)src_ptr[ 2 * (int)src_pixels_per_line + 2 * (int)pixel_step] * vp8_filter[28]) +
           ((int)src_ptr[ 2 * (int)src_pixels_per_line + 3 * (int)pixel_step] * vp8_filter[29]) +
           ((int)src_ptr[ 3 * (int)src_pixels_per_line - 2 * (int)pixel_step] * vp8_filter[30]) +
           ((int)src_ptr[ 3 * (int)src_pixels_per_line - 1 * (int)pixel_step] * vp8_filter[31]) +
           ((int)src_ptr[ 3 * (int)src_pixels_per_line + 0 * (int)pixel_step] * vp8_filter[32]) +
           ((int)src_ptr[ 3 * (int)src_pixels_per_line + 1 * (int)pixel_step] * vp8_filter[33]) +
           ((int)src_ptr[ 3 * (int)src_pixels_per_line + 2 * (int)pixel_step] * vp8_filter[34]) +
           ((int)src_ptr[ 3 * (int)src_pixels_per_line + 3 * (int)pixel_step] * vp8_filter[35]) +
           (VP8_FILTER_WEIGHT >> 1);      /* Rounding */
#endif
    Temp = Temp >> VP8_FILTER_SHIFT;

    if (Temp < 0)
        Temp = 0;
    else if (Temp > 255)
        Temp = 255;

    *output_ptr = Temp;
}

static void filter_edge_pixel
(
    unsigned char *src_ptr,
    unsigned char *output_ptr,
    unsigned int src_pixels_per_line,
    unsigned int pixel_step,
    int xoffset,
    int yoffset
)
{
    const short *vp8_filter=vp8_sub_pel_filters_ns[xoffset+8*yoffset];
    filter_non_separable(src_ptr, output_ptr, src_pixels_per_line, pixel_step, vp8_filter);
}

static void get_sobel_grads(unsigned char *src_ptr, int width, int height,
                            unsigned int src_pixels_per_line,
                            unsigned int *sum_g)
{
    /* Assume that the block always has extension of at least 1 */
    int i, j;
    int gx, gy, gd, ga;
    unsigned char *prev = src_ptr-src_pixels_per_line;
    unsigned char *prev2 = src_ptr-2*src_pixels_per_line;
    unsigned char *curr = src_ptr;
    unsigned char *next = src_ptr+src_pixels_per_line;
    unsigned char *next2 = src_ptr+2*src_pixels_per_line;
    sum_g[0] = sum_g[1] = sum_g[2] = sum_g[3] = 0;
    for (i=0; i<height; ++i)
    {
        for (j=0; j<width; ++j)
        {
            gx = abs((prev[1]-prev[-1])+((curr[1]-curr[-1])*2)+(next[1]-next[-1]));
            gy = abs((prev[-1]-next[-1])+((prev[0]-next[0])*2)+(prev[1]-next[1]));
            gd = abs((curr[2]-prev2[0])+((next[1]-prev[-1])*2)+(next2[0]-curr[-2]));
            ga = abs((curr[2]-next2[0])+((prev[1]-next[-1])*2)+(prev2[0]-curr[-2]));
            sum_g[0] += (gx>EDGE_GRAD_THRESH*4);
            sum_g[1] += (gy>EDGE_GRAD_THRESH*4);
            sum_g[2] += (gd>EDGE_GRAD_THRESH*4);
            sum_g[3] += (ga>EDGE_GRAD_THRESH*4);
            prev++;
            prev2++;
            curr++;
            next++;
            next2++;
        }
        prev  += src_pixels_per_line-width;
        curr  += src_pixels_per_line-width;
        next  += src_pixels_per_line-width;
        prev2 += src_pixels_per_line-width;
        next2 += src_pixels_per_line-width;
    }
}

static int edge_pixel_detected(unsigned char *src_ptr, int src_pitch)
{
    unsigned int ng[4];
    get_sobel_grads(src_ptr, 2, 2, src_pitch, ng);
    return (ng[0] + ng[1] + ng[2] + ng[3] > EDGE_GRADS2X2_THRESH);
}

static int edge_pixel_detected_simple(unsigned char *src_ptr, int src_pitch)
{
    int gmax, gmin, gmax2, gmin2;
    if (src_ptr[0]>src_ptr[1])
    {
        gmax=src_ptr[0];
        gmin=src_ptr[1];
    }
    else
    {
        gmax=src_ptr[1];
        gmin=src_ptr[0];
    }
    src_ptr += src_pitch;
    if (src_ptr[0]>src_ptr[1])
    {
        gmax2=src_ptr[0];
        gmin2=src_ptr[1];
    }
    else
    {
        gmax2=src_ptr[1];
        gmin2=src_ptr[0];
    }
    if (gmax2>gmax) gmax=gmax2;
    if (gmin2<gmin) gmin=gmin2;
    return (gmax - gmin > EDGE_SIMPLE_THRESH);
}

void vp8_edge_pixel_interpolation
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int width,
    int height,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int  dst_pitch
)
{
    unsigned char *sp = src_ptr;
    unsigned char *dp = dst_ptr;
    int i, j;
    for (i = 0; i < height; ++i, sp+=src_pixels_per_line-width, dp+=dst_pitch-width)
        for (j = 0; j < width; ++j, ++sp, ++dp)
        {
            if (edge_pixel_detected(sp, src_pixels_per_line))
            {
                filter_edge_pixel(sp, dp, src_pixels_per_line, 1, xoffset, yoffset);
            }
        }
}
#endif  // EDGE_PIXEL_FILTER

/*
 * The only functional difference between filter_block2d_second_pass()
 * and this function is that filter_block2d_second_pass() does a sixtap
 * filter on the input and stores it in the output. This function
 * (filter_block2d_second_pass_avg()) does a sixtap filter on the input,
 * and then averages that with the content already present in the output
 * ((filter_result + dest + 1) >> 1) and stores that in the output.
 */
static void filter_block2d_second_pass_avg
(
    int *src_ptr,
    unsigned char *output_ptr,
    int output_pitch,
    unsigned int src_pixels_per_line,
    unsigned int pixel_step,
    unsigned int output_height,
    unsigned int output_width,
    const short *vp8_filter
)
{
    unsigned int i, j;
    int  Temp;

    for (i = 0; i < output_height; i++)
    {
        for (j = 0; j < output_width; j++)
        {
            /* Apply filter */
#if INTERP_EXTEND == 3
            Temp = ((int)src_ptr[-2 * (int)pixel_step] * vp8_filter[0]) +
                   ((int)src_ptr[-1 * (int)pixel_step] * vp8_filter[1]) +
                   ((int)src_ptr[0]                    * vp8_filter[2]) +
                   ((int)src_ptr[pixel_step]           * vp8_filter[3]) +
                   ((int)src_ptr[2*pixel_step]         * vp8_filter[4]) +
                   ((int)src_ptr[3*pixel_step]         * vp8_filter[5]) +
                   (VP8_FILTER_WEIGHT >> 1);   /* Rounding */
#elif INTERP_EXTEND == 4
            Temp = ((int)src_ptr[-3 * (int)pixel_step] * vp8_filter[0]) +
                   ((int)src_ptr[-2 * (int)pixel_step] * vp8_filter[1]) +
                   ((int)src_ptr[-1 * (int)pixel_step] * vp8_filter[2]) +
                   ((int)src_ptr[0]                    * vp8_filter[3]) +
                   ((int)src_ptr[pixel_step]           * vp8_filter[4]) +
                   ((int)src_ptr[2 * pixel_step]       * vp8_filter[5]) +
                   ((int)src_ptr[3 * pixel_step]       * vp8_filter[6]) +
                   ((int)src_ptr[4 * pixel_step]       * vp8_filter[7]) +
                   (VP8_FILTER_WEIGHT >> 1);      /* Rounding */
#elif INTERP_EXTEND == 5
            Temp = ((int)src_ptr[-4 * (int)pixel_step] * vp8_filter[0]) +
                   ((int)src_ptr[-3 * (int)pixel_step] * vp8_filter[1]) +
                   ((int)src_ptr[-2 * (int)pixel_step] * vp8_filter[2]) +
                   ((int)src_ptr[-1 * (int)pixel_step] * vp8_filter[3]) +
                   ((int)src_ptr[0]                    * vp8_filter[4]) +
                   ((int)src_ptr[pixel_step]           * vp8_filter[5]) +
                   ((int)src_ptr[2 * pixel_step]       * vp8_filter[6]) +
                   ((int)src_ptr[3 * pixel_step]       * vp8_filter[7]) +
                   ((int)src_ptr[4 * pixel_step]       * vp8_filter[8]) +
                   ((int)src_ptr[5 * pixel_step]       * vp8_filter[9]) +
                   (VP8_FILTER_WEIGHT >> 1);      /* Rounding */
#endif

            /* Normalize back to 0-255 */
            Temp = Temp >> VP8_FILTER_SHIFT;

            if (Temp < 0)
                Temp = 0;
            else if (Temp > 255)
                Temp = 255;

            output_ptr[j] = (unsigned char) ((output_ptr[j] + Temp + 1) >> 1);
            src_ptr++;
        }

        /* Start next row */
        src_ptr    += src_pixels_per_line - output_width;
        output_ptr += output_pitch;
    }
}

static void filter_block2d
(
    unsigned char  *src_ptr,
    unsigned char  *output_ptr,
    unsigned int src_pixels_per_line,
    int output_pitch,
    const short  *HFilter,
    const short  *VFilter
)
{
    int FData[(3+INTERP_EXTEND*2)*4]; /* Temp data buffer used in filtering */

    /* First filter 1-D horizontally... */
    filter_block2d_first_pass(src_ptr - ((INTERP_EXTEND-1) * src_pixels_per_line), FData, src_pixels_per_line, 1,
                              3+INTERP_EXTEND*2, 4, HFilter);

    /* then filter verticaly... */
    filter_block2d_second_pass(FData + 4*(INTERP_EXTEND-1), output_ptr, output_pitch, 4, 4, 4, 4, VFilter);
}


void vp8_sixtap_predict_c
(
    unsigned char  *src_ptr,
    int   src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int dst_pitch
)
{
    const short  *HFilter;
    const short  *VFilter;

    HFilter = vp8_sub_pel_filters[xoffset];   /* 6 tap */
    VFilter = vp8_sub_pel_filters[yoffset];   /* 6 tap */

    filter_block2d(src_ptr, dst_ptr, src_pixels_per_line, dst_pitch, HFilter, VFilter);
#if CONFIG_ENHANCED_INTERP && EDGE_PIXEL_FILTER > 0
    vp8_edge_pixel_interpolation(src_ptr, src_pixels_per_line, 4, 4,
                                 xoffset, yoffset, dst_ptr, dst_pitch);
#endif
}
void vp8_sixtap_predict8x8_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int  dst_pitch
)
{
    const short  *HFilter;
    const short  *VFilter;
    // int FData[(7+INTERP_EXTEND*2)*16];   /* Temp data buffer used in filtering */
    int FData[(7+INTERP_EXTEND*2)*8];   /* Temp data buffer used in filtering */

    HFilter = vp8_sub_pel_filters[xoffset];   /* 6 tap */
    VFilter = vp8_sub_pel_filters[yoffset];   /* 6 tap */

    /* First filter 1-D horizontally... */
    filter_block2d_first_pass(src_ptr - ((INTERP_EXTEND-1) * src_pixels_per_line), FData, src_pixels_per_line, 1,
                              7+INTERP_EXTEND*2, 8, HFilter);


    /* then filter verticaly... */
    filter_block2d_second_pass(FData + 8*(INTERP_EXTEND-1), dst_ptr, dst_pitch, 8, 8, 8, 8, VFilter);

#if CONFIG_ENHANCED_INTERP && EDGE_PIXEL_FILTER > 0
    vp8_edge_pixel_interpolation(src_ptr, src_pixels_per_line, 8, 8,
                                 xoffset, yoffset, dst_ptr, dst_pitch);
#endif
}

void vp8_sixtap_predict_avg8x8_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int  dst_pitch
)
{
    const short  *HFilter;
    const short  *VFilter;
    // int FData[(7+INTERP_EXTEND*2)*16];   /* Temp data buffer used in filtering */
    int FData[(7+INTERP_EXTEND*2)*8];   /* Temp data buffer used in filtering */

    HFilter = vp8_sub_pel_filters[xoffset];   /* 6 tap */
    VFilter = vp8_sub_pel_filters[yoffset];   /* 6 tap */

    /* First filter 1-D horizontally... */
    filter_block2d_first_pass(src_ptr - ((INTERP_EXTEND-1) * src_pixels_per_line), FData, src_pixels_per_line, 1,
                              7+INTERP_EXTEND*2, 8, HFilter);

    /* then filter verticaly... */
    filter_block2d_second_pass_avg(FData + 8*(INTERP_EXTEND-1), dst_ptr, dst_pitch, 8, 8, 8, 8, VFilter);
#if CONFIG_ENHANCED_INTERP && EDGE_PIXEL_FILTER > 0
    vp8_edge_pixel_interpolation(src_ptr, src_pixels_per_line, 8, 8,
                                 xoffset, yoffset, dst_ptr, dst_pitch);
#endif
}

void vp8_sixtap_predict8x4_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int  dst_pitch
)
{
    const short  *HFilter;
    const short  *VFilter;
    // int FData[(7+INTERP_EXTEND*2)*16];   /* Temp data buffer used in filtering */
    int FData[(3+INTERP_EXTEND*2)*8];   /* Temp data buffer used in filtering */

    HFilter = vp8_sub_pel_filters[xoffset];   /* 6 tap */
    VFilter = vp8_sub_pel_filters[yoffset];   /* 6 tap */

    /* First filter 1-D horizontally... */
    filter_block2d_first_pass(src_ptr - ((INTERP_EXTEND-1) * src_pixels_per_line), FData, src_pixels_per_line, 1,
                              3+INTERP_EXTEND*2, 8, HFilter);


    /* then filter verticaly... */
    filter_block2d_second_pass(FData + 8*(INTERP_EXTEND-1), dst_ptr, dst_pitch, 8, 8, 4, 8, VFilter);

#if CONFIG_ENHANCED_INTERP && EDGE_PIXEL_FILTER > 0
    vp8_edge_pixel_interpolation(src_ptr, src_pixels_per_line, 8, 4,
                                 xoffset, yoffset, dst_ptr, dst_pitch);
#endif
}

void vp8_sixtap_predict16x16_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int  dst_pitch
)
{
    const short  *HFilter;
    const short  *VFilter;
    // int FData[(15+INTERP_EXTEND*2)*24];   /* Temp data buffer used in filtering */
    int FData[(15+INTERP_EXTEND*2)*16];  /* Temp data buffer used in filtering */


    HFilter = vp8_sub_pel_filters[xoffset];   /* 6 tap */
    VFilter = vp8_sub_pel_filters[yoffset];   /* 6 tap */

    /* First filter 1-D horizontally... */
    filter_block2d_first_pass(src_ptr - ((INTERP_EXTEND-1) * src_pixels_per_line), FData, src_pixels_per_line, 1,
                              15+INTERP_EXTEND*2, 16, HFilter);

    /* then filter verticaly... */
    filter_block2d_second_pass(FData + 16*(INTERP_EXTEND-1), dst_ptr, dst_pitch, 16, 16, 16, 16, VFilter);

#if CONFIG_ENHANCED_INTERP && EDGE_PIXEL_FILTER > 0
    vp8_edge_pixel_interpolation(src_ptr, src_pixels_per_line, 16, 16,
                                 xoffset, yoffset, dst_ptr, dst_pitch);
#endif
}

void vp8_sixtap_predict_avg16x16_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int  dst_pitch
)
{
    const short  *HFilter;
    const short  *VFilter;
    // int FData[(15+INTERP_EXTEND*2)*24];   /* Temp data buffer used in filtering */
    int FData[(15+INTERP_EXTEND*2)*16];  /* Temp data buffer used in filtering */

    HFilter = vp8_sub_pel_filters[xoffset];   /* 6 tap */
    VFilter = vp8_sub_pel_filters[yoffset];   /* 6 tap */

    /* First filter 1-D horizontally... */
    filter_block2d_first_pass(src_ptr - ((INTERP_EXTEND-1) * src_pixels_per_line), FData,
                              src_pixels_per_line, 1, 15+INTERP_EXTEND*2, 16, HFilter);

    /* then filter verticaly... */
    filter_block2d_second_pass_avg(FData + 16*(INTERP_EXTEND-1), dst_ptr, dst_pitch,
                                   16, 16, 16, 16, VFilter);
#if CONFIG_ENHANCED_INTERP && EDGE_PIXEL_FILTER > 0
    vp8_edge_pixel_interpolation(src_ptr, src_pixels_per_line, 16, 16,
                                 xoffset, yoffset, dst_ptr, dst_pitch);
#endif
}

/****************************************************************************
 *
 *  ROUTINE       : filter_block2d_bil_first_pass
 *
 *  INPUTS        : UINT8  *src_ptr    : Pointer to source block.
 *                  UINT32  src_stride : Stride of source block.
 *                  UINT32  height     : Block height.
 *                  UINT32  width      : Block width.
 *                  INT32  *vp8_filter : Array of 2 bi-linear filter taps.
 *
 *  OUTPUTS       : INT32  *dst_ptr    : Pointer to filtered block.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Applies a 1-D 2-tap bi-linear filter to the source block
 *                  in the horizontal direction to produce the filtered output
 *                  block. Used to implement first-pass of 2-D separable filter.
 *
 *  SPECIAL NOTES : Produces INT32 output to retain precision for next pass.
 *                  Two filter taps should sum to VP8_FILTER_WEIGHT.
 *
 ****************************************************************************/
static void filter_block2d_bil_first_pass
(
    unsigned char  *src_ptr,
    unsigned short *dst_ptr,
    unsigned int    src_stride,
    unsigned int    height,
    unsigned int    width,
    const short    *vp8_filter
)
{
    unsigned int i, j;

    for (i = 0; i < height; i++)
    {
        for (j = 0; j < width; j++)
        {
            /* Apply bilinear filter */
            dst_ptr[j] = (((int)src_ptr[0] * vp8_filter[0]) +
                          ((int)src_ptr[1] * vp8_filter[1]) +
                          (VP8_FILTER_WEIGHT / 2)) >> VP8_FILTER_SHIFT;
            src_ptr++;
        }

        /* Next row... */
        src_ptr += src_stride - width;
        dst_ptr += width;
    }
}

/****************************************************************************
 *
 *  ROUTINE       : filter_block2d_bil_second_pass
 *
 *  INPUTS        : INT32  *src_ptr    : Pointer to source block.
 *                  UINT32  dst_pitch  : Destination block pitch.
 *                  UINT32  height     : Block height.
 *                  UINT32  width      : Block width.
 *                  INT32  *vp8_filter : Array of 2 bi-linear filter taps.
 *
 *  OUTPUTS       : UINT16 *dst_ptr    : Pointer to filtered block.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Applies a 1-D 2-tap bi-linear filter to the source block
 *                  in the vertical direction to produce the filtered output
 *                  block. Used to implement second-pass of 2-D separable filter.
 *
 *  SPECIAL NOTES : Requires 32-bit input as produced by filter_block2d_bil_first_pass.
 *                  Two filter taps should sum to VP8_FILTER_WEIGHT.
 *
 ****************************************************************************/
static void filter_block2d_bil_second_pass
(
    unsigned short *src_ptr,
    unsigned char  *dst_ptr,
    int             dst_pitch,
    unsigned int    height,
    unsigned int    width,
    const short    *vp8_filter
)
{
    unsigned int  i, j;
    int  Temp;

    for (i = 0; i < height; i++)
    {
        for (j = 0; j < width; j++)
        {
            /* Apply filter */
            Temp = ((int)src_ptr[0]     * vp8_filter[0]) +
                   ((int)src_ptr[width] * vp8_filter[1]) +
                   (VP8_FILTER_WEIGHT / 2);
            dst_ptr[j] = (unsigned int)(Temp >> VP8_FILTER_SHIFT);
            src_ptr++;
        }

        /* Next row... */
        dst_ptr += dst_pitch;
    }
}

/*
 * As before for filter_block2d_second_pass_avg(), the functional difference
 * between filter_block2d_bil_second_pass() and filter_block2d_bil_second_pass_avg()
 * is that filter_block2d_bil_second_pass() does a bilinear filter on input
 * and stores the result in output; filter_block2d_bil_second_pass_avg(),
 * instead, does a bilinear filter on input, averages the resulting value
 * with the values already present in the output and stores the result of
 * that back into the output ((filter_result + dest + 1) >> 1).
 */
static void filter_block2d_bil_second_pass_avg
(
    unsigned short *src_ptr,
    unsigned char  *dst_ptr,
    int             dst_pitch,
    unsigned int    height,
    unsigned int    width,
    const short    *vp8_filter
)
{
    unsigned int  i, j;
    int  Temp;

    for (i = 0; i < height; i++)
    {
        for (j = 0; j < width; j++)
        {
            /* Apply filter */
            Temp = ((int)src_ptr[0]     * vp8_filter[0]) +
                   ((int)src_ptr[width] * vp8_filter[1]) +
                   (VP8_FILTER_WEIGHT / 2);
            dst_ptr[j] = (unsigned int)(((Temp >> VP8_FILTER_SHIFT) + dst_ptr[j] + 1) >> 1);
            src_ptr++;
        }

        /* Next row... */
        dst_ptr += dst_pitch;
    }
}

/****************************************************************************
 *
 *  ROUTINE       : filter_block2d_bil
 *
 *  INPUTS        : UINT8  *src_ptr          : Pointer to source block.
 *                  UINT32  src_pitch        : Stride of source block.
 *                  UINT32  dst_pitch        : Stride of destination block.
 *                  INT32  *HFilter          : Array of 2 horizontal filter taps.
 *                  INT32  *VFilter          : Array of 2 vertical filter taps.
 *                  INT32  Width             : Block width
 *                  INT32  Height            : Block height
 *
 *  OUTPUTS       : UINT16 *dst_ptr       : Pointer to filtered block.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : 2-D filters an input block by applying a 2-tap
 *                  bi-linear filter horizontally followed by a 2-tap
 *                  bi-linear filter vertically on the result.
 *
 *  SPECIAL NOTES : The largest block size can be handled here is 16x16
 *
 ****************************************************************************/
static void filter_block2d_bil
(
    unsigned char *src_ptr,
    unsigned char *dst_ptr,
    unsigned int   src_pitch,
    unsigned int   dst_pitch,
    const short   *HFilter,
    const short   *VFilter,
    int            Width,
    int            Height
)
{

    unsigned short FData[17*16];    /* Temp data buffer used in filtering */

    /* First filter 1-D horizontally... */
    filter_block2d_bil_first_pass(src_ptr, FData, src_pitch, Height + 1, Width, HFilter);

    /* then 1-D vertically... */
    filter_block2d_bil_second_pass(FData, dst_ptr, dst_pitch, Height, Width, VFilter);
}

static void filter_block2d_bil_avg
(
    unsigned char *src_ptr,
    unsigned char *dst_ptr,
    unsigned int   src_pitch,
    unsigned int   dst_pitch,
    const short   *HFilter,
    const short   *VFilter,
    int            Width,
    int            Height
)
{
    unsigned short FData[17*16];    /* Temp data buffer used in filtering */

    /* First filter 1-D horizontally... */
    filter_block2d_bil_first_pass(src_ptr, FData, src_pitch, Height + 1, Width, HFilter);

    /* then 1-D vertically... */
    filter_block2d_bil_second_pass_avg(FData, dst_ptr, dst_pitch, Height, Width, VFilter);
}

void vp8_bilinear_predict4x4_c
(
    unsigned char  *src_ptr,
    int   src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int dst_pitch
)
{
    const short *HFilter;
    const short *VFilter;

    HFilter = vp8_bilinear_filters[xoffset];
    VFilter = vp8_bilinear_filters[yoffset];
#if 0
    {
        int i;
        unsigned char temp1[16];
        unsigned char temp2[16];

        bilinear_predict4x4_mmx(src_ptr, src_pixels_per_line, xoffset, yoffset, temp1, 4);
        filter_block2d_bil(src_ptr, temp2, src_pixels_per_line, 4, HFilter, VFilter, 4, 4);

        for (i = 0; i < 16; i++)
        {
            if (temp1[i] != temp2[i])
            {
                bilinear_predict4x4_mmx(src_ptr, src_pixels_per_line, xoffset, yoffset, temp1, 4);
                filter_block2d_bil(src_ptr, temp2, src_pixels_per_line, 4, HFilter, VFilter, 4, 4);
            }
        }
    }
#endif
    filter_block2d_bil(src_ptr, dst_ptr, src_pixels_per_line, dst_pitch, HFilter, VFilter, 4, 4);

}

void vp8_bilinear_predict8x8_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int  dst_pitch
)
{
    const short *HFilter;
    const short *VFilter;

    HFilter = vp8_bilinear_filters[xoffset];
    VFilter = vp8_bilinear_filters[yoffset];

    filter_block2d_bil(src_ptr, dst_ptr, src_pixels_per_line, dst_pitch, HFilter, VFilter, 8, 8);

}

void vp8_bilinear_predict_avg8x8_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int  dst_pitch
)
{
    const short *HFilter;
    const short *VFilter;

    HFilter = vp8_bilinear_filters[xoffset];
    VFilter = vp8_bilinear_filters[yoffset];

    filter_block2d_bil_avg(src_ptr, dst_ptr, src_pixels_per_line,
                           dst_pitch, HFilter, VFilter, 8, 8);
}

void vp8_bilinear_predict8x4_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int  dst_pitch
)
{
    const short *HFilter;
    const short *VFilter;

    HFilter = vp8_bilinear_filters[xoffset];
    VFilter = vp8_bilinear_filters[yoffset];

    filter_block2d_bil(src_ptr, dst_ptr, src_pixels_per_line, dst_pitch, HFilter, VFilter, 8, 4);

}

void vp8_bilinear_predict16x16_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int  dst_pitch
)
{
    const short *HFilter;
    const short *VFilter;

    HFilter = vp8_bilinear_filters[xoffset];
    VFilter = vp8_bilinear_filters[yoffset];

    filter_block2d_bil(src_ptr, dst_ptr, src_pixels_per_line, dst_pitch, HFilter, VFilter, 16, 16);
}

void vp8_bilinear_predict_avg16x16_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int  dst_pitch
)
{
    const short *HFilter;
    const short *VFilter;

    HFilter = vp8_bilinear_filters[xoffset];
    VFilter = vp8_bilinear_filters[yoffset];

    filter_block2d_bil_avg(src_ptr, dst_ptr, src_pixels_per_line,
                           dst_pitch, HFilter, VFilter, 16, 16);
}
