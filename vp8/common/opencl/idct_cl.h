/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef IDCT_OPENCL_H
#define IDCT_OPENCL_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "vp8_opencl.h"
#include "vp8/common/blockd.h"

#define prototype_second_order_cl(sym) \
    void sym(BLOCKD *b)

#define prototype_idct_cl(sym) \
    void sym(BLOCKD *b, int pitch)

#define prototype_idct_scalar_add_cl(sym) \
    void sym(BLOCKD *b, cl_int use_diff, int diff_offset, int qcoeff_offset, \
             int pred_offset, unsigned char *output, cl_mem out_mem, int out_offset, size_t out_size, \
             int pitch, int stride)\


extern prototype_idct_cl(vp8_short_idct4x4llm_1_cl);
extern prototype_idct_cl(vp8_short_idct4x4llm_cl);
extern prototype_idct_scalar_add_cl(vp8_dc_only_idct_add_cl);

extern prototype_second_order_cl(vp8_short_inv_walsh4x4_1_cl);
extern prototype_second_order_cl(vp8_short_inv_walsh4x4_cl);

#ifdef	__cplusplus
}
#endif

#endif
