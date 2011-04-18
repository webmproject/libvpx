/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef DEQUANTIZE_CL_H
#define DEQUANTIZE_CL_H

#ifdef  __cplusplus
extern "C" {
#endif

#include "vp8/decoder/onyxd_int.h"
#include "vp8/decoder/dequantize.h"
#include "vp8/common/opencl/vp8_opencl.h"

#define prototype_dequant_block_cl(sym) \
    void sym(BLOCKD *x)

#define prototype_dequant_idct_add_cl(sym) \
    void sym(BLOCKD *b, unsigned char *dest_base,cl_mem dest_mem, int dest_offset, size_t dest_size, int q_offset, \
             int pred_offset, int pitch, int stride, \
             vp8_dequant_idct_add_fn_t idct_add)

#define prototype_dequant_dc_idct_add_cl(sym) \
    void sym(BLOCKD* b, int qcoeff_offset, \
             int pred_offset, unsigned char *dest_base, int dst_offset, \
             int pitch, int stride, \
             int dc)

#define prototype_dequant_dc_idct_add_y_block_cl(sym) \
    void sym(BLOCKD *b, \
             unsigned char *dst_base, cl_mem dst_mem, int dst_off,\
             int stride, char *eobs, int dc_offset)

#define prototype_dequant_idct_add_y_block_cl(sym) \
    void sym(VP8D_COMP *pbi, MACROBLOCKD *xd)

#define prototype_dequant_idct_add_uv_block_cl(sym) \
    void sym(VP8D_COMP *pbi, MACROBLOCKD *xd, \
        vp8_dequant_idct_add_uv_block_fn_t idct_add_uv_block)


    
extern prototype_dequant_block_cl(vp8_dequantize_b_cl);

//CL functions
extern prototype_dequant_idct_add_cl(vp8_dequant_idct_add_cl);

//C functions
extern prototype_dequant_dc_idct_add_cl(vp8_dequant_dc_idct_add_cl);


//Might be CL... check implementation.
extern prototype_dequant_dc_idct_add_y_block_cl(vp8_dequant_dc_idct_add_y_block_cl);
extern prototype_dequant_idct_add_y_block_cl(vp8_dequant_idct_add_y_block_cl);
extern prototype_dequant_idct_add_uv_block_cl(vp8_dequant_idct_add_uv_block_cl);



extern const char *dequantCompileOptions;
extern const char *dequant_cl_file_name;

#ifdef  __cplusplus
}
#endif

#endif
