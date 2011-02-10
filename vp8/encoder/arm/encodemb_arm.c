/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp8/encoder/encodemb.h"
#include "vp8/common/reconinter.h"
#include "vp8/encoder/quantize.h"
#include "vp8/common/invtrans.h"
#include "vp8/common/recon.h"
#include "vp8/common/reconintra.h"
#include "vp8/encoder/dct.h"
#include "vpx_mem/vpx_mem.h"

extern void vp8_subtract_b_neon_func(short *diff, unsigned char *src, unsigned char *pred, int stride, int pitch);

void vp8_subtract_b_neon(BLOCK *be, BLOCKD *bd, int pitch)
{
    unsigned char *src_ptr = (*(be->base_src) + be->src);
    short *diff_ptr = be->src_diff;
    unsigned char *pred_ptr = bd->predictor;
    int src_stride = be->src_stride;

    vp8_subtract_b_neon_func(diff_ptr, src_ptr, pred_ptr, src_stride, pitch);
}
