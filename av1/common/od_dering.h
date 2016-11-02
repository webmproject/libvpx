/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#if !defined(_dering_H)
#define _dering_H (1)

// clang-format off

#include "odintrin.h"

#if defined(DAALA_ODINTRIN)
#include "filter.h"
typedef int16_t od_dering_in;
#endif

#define OD_DERINGSIZES (2)

#define OD_DERING_NBLOCKS (OD_BSIZE_MAX / 8)

#define OD_FILT_BORDER (3)
#define OD_FILT_BSTRIDE (OD_BSIZE_MAX + 2 * OD_FILT_BORDER)

extern const int OD_DIRECTION_OFFSETS_TABLE[8][3];

typedef int (*od_filter_dering_direction_func)(int16_t *y, int ystride,
                                               const int16_t *in, int threshold,
                                               int dir);
typedef void (*od_filter_dering_orthogonal_func)(int16_t *y, int ystride,
                                                 const int16_t *in,
                                                 int threshold, int dir);
void copy_blocks_16bit(int16_t *dst, int dstride, int16_t *src,
    unsigned char (*bskip)[2], int dering_count, int bsize);

void od_dering(int16_t *y, const od_dering_in *x, int xstride,
               int nvb, int nhb, int sbx, int sby, int nhsb, int nvsb, int xdec,
               int dir[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS], int pli,
               unsigned char (*bskip)[2], int skip_stride, int threshold,
               int coeff_shift);
int od_filter_dering_direction_4x4_c(int16_t *y, int ystride, const int16_t *in,
                                     int threshold, int dir);
int od_filter_dering_direction_8x8_c(int16_t *y, int ystride, const int16_t *in,
                                     int threshold, int dir);
void od_filter_dering_orthogonal_4x4_c(int16_t *y, int ystride,
                                       const int16_t *in, int threshold,
                                       int dir);
void od_filter_dering_orthogonal_8x8_c(int16_t *y, int ystride,
                                       const int16_t *in, int threshold,
                                       int dir);
#endif
