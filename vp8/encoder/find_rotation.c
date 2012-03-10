/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "block.h"
#include "variance.h"

#if CONFIG_ROTATION

int vp8_find_best_rotation(MACROBLOCK *x, BLOCK *b, BLOCKD *d, int_mv *bestmv,
                           int_mv *ref_mv, int *bri, int error_per_bit,
                           const vp8_variance_fn_ptr_t *vfp, int *mvcost[2],
                           int *distortion, unsigned int *sse1)
{
    unsigned char *z = (*(b->base_src) + b->src);

    int ri;

    int y_stride;

    unsigned int besterr;
    int br = bestmv->as_mv.row;
    int bc = bestmv->as_mv.col;
    unsigned char *y = *(d->base_pre) + d->pre + br * d->pre_stride + bc;
    y_stride = d->pre_stride;

    // calculate central point error
    besterr = vfp->vf(y, y_stride, z, b->src_stride, sse1);
    *distortion = besterr;

    // find the best matching rotation
    *bri = 5;
    for (ri = 0; ri < ROTATIONS; ri++)
    {
        unsigned int this_err;
        unsigned char pb[256];
        predict_rotated_16x16(ri, y, y_stride, pb, 16);
        this_err = vfp->vf(pb, 16, z, b->src_stride, sse1);

        if (this_err < besterr)
        {
            *bri = ri;
            besterr = this_err;
        }
    }
    *sse1 = besterr;
    *distortion = besterr;

    return 0;
}

#endif
