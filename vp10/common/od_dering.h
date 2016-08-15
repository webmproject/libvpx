/*Daala video codec
Copyright (c) 2003-2010 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#if !defined(_dering_H)
# define _dering_H (1)

# include "odintrin.h"

# if defined(DAALA_ODINTRIN)
#  include "filter.h"
typedef int16_t od_dering_in;
# endif

#define OD_DERINGSIZES (2)

#define OD_DERING_NO_CHECK_OVERLAP (0)
#define OD_DERING_CHECK_OVERLAP (1)

#define OD_DERING_LEVELS (6)
extern const double OD_DERING_GAIN_TABLE[OD_DERING_LEVELS];

#define OD_DERING_NBLOCKS (OD_BSIZE_MAX/8)

#define OD_FILT_BORDER (3)
#define OD_FILT_BSTRIDE (OD_BSIZE_MAX + 2*OD_FILT_BORDER)

extern const int OD_DIRECTION_OFFSETS_TABLE[8][3];

typedef void (*od_filter_dering_direction_func)(int16_t *y, int ystride,
 const int16_t *in, int threshold, int dir);
typedef void (*od_filter_dering_orthogonal_func)(int16_t *y, int ystride,
 const int16_t *in, const od_dering_in *x, int xstride, int threshold,
 int dir);

struct od_dering_opt_vtbl {
  od_filter_dering_direction_func filter_dering_direction[OD_DERINGSIZES];
  od_filter_dering_orthogonal_func filter_dering_orthogonal[OD_DERINGSIZES];
};
typedef struct od_dering_opt_vtbl od_dering_opt_vtbl;


void od_dering(const od_dering_opt_vtbl *vtbl, int16_t *y, int ystride,
 const od_dering_in *x, int xstride, int nvb, int nhb, int sbx, int sby,
 int nhsb, int nvsb, int xdec, int dir[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS],
 int pli, unsigned char *bskip, int skip_stride, int threshold, int overlap,
 int coeff_shift);
void od_filter_dering_direction_c(int16_t *y, int ystride, const int16_t *in,
 int ln, int threshold, int dir);
void od_filter_dering_orthogonal_c(int16_t *y, int ystride, const int16_t *in,
 const od_dering_in *x, int xstride, int ln, int threshold, int dir);

extern const od_dering_opt_vtbl OD_DERING_VTBL_C;

void od_filter_dering_direction_4x4_c(int16_t *y, int ystride,
 const int16_t *in, int threshold, int dir);
void od_filter_dering_direction_8x8_c(int16_t *y, int ystride,
 const int16_t *in, int threshold, int dir);
void od_filter_dering_orthogonal_4x4_c(int16_t *y, int ystride,
 const int16_t *in, const od_dering_in *x, int xstride, int threshold,
 int dir);
void od_filter_dering_orthogonal_8x8_c(int16_t *y, int ystride,
 const int16_t *in, const od_dering_in *x, int xstride, int threshold,
 int dir);

#endif
