/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp8/decoder/onyxd_int.h"
#include "vpx_ports/config.h"
#include "../../common/idct.h"
#include "vp8/common/opencl/blockd_cl.h"
#include "dequantize_cl.h"

//change q/dq/pre/eobs/dc to offsets
void vp8_dequant_dc_idct_add_y_block_cl(
    BLOCKD *b,
    unsigned char *dst_base, //xd->dst.buffer_alloc
    cl_mem dst_mem,
    int dst_off,
    int stride,         //xd->dst.y_stride
    char *eobs,         //xd->eobs
    int dc_offset       //xd->block[24].diff_offset
)
{
    int i, j;
    int q_offset = 0;
    int pre_offset = 0;
    int dst_offset = 0;
    unsigned char *dst = dst_base+dst_off;
    size_t dst_size = 16*(stride+1);

    vp8_cl_block_prep(b, QCOEFF|DEQUANT|DIFF|PREDICTOR);
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            if (*eobs++ > 1){
                vp8_dequant_dc_idct_add_cl (b, q_offset, pre_offset, dst, dst_offset, 16, stride, dc_offset);
            }
            else{
                vp8_dc_only_idct_add_cl(b, CL_TRUE, dc_offset, 0, pre_offset, dst, NULL, dst_offset, dst_size, 16, stride);
            }

            q_offset   += 16;
            pre_offset += 4;
            dst_offset += 4;
            dc_offset++;
        }

        pre_offset += 64 - 16;
        dst_offset += 4*stride - 16;
    }

    vp8_cl_block_finish(b, QCOEFF);

}

void vp8_dequant_idct_add_y_block_cl (VP8D_COMP *pbi, MACROBLOCKD *xd)
{
    int i, j;

    short *q = xd->qcoeff;
    int q_offset = 0;
    int pre_offset = 0;
    cl_mem dst_mem = xd->dst.buffer_mem;
    unsigned char *dst = xd->dst.buffer_alloc;
    int dst_offset = xd->dst.y_buffer - dst;
    int stride = xd->dst.y_stride;
    char *eobs = xd->eobs;
    int dst_size = 16 * (stride + 1);


    vp8_cl_mb_prep(xd,PREDICTOR|DIFF|QCOEFF);
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            if (*eobs++ > 1){
                vp8_cl_block_prep(&xd->block[0], DEQUANT);
                vp8_dequant_idct_add_cl(&xd->block[0], dst, dst_mem, dst_offset, dst_size+dst_offset, q_offset, pre_offset, 16, stride, pbi->dequant.idct_add);
                vp8_cl_block_finish(&xd->block[0], QCOEFF);
            }
            else
            {
                vp8_cl_block_prep(&xd->block[0], DEQUANT);
                vp8_dc_only_idct_add_cl(&xd->block[0], CL_FALSE, 0, q_offset, pre_offset, dst, dst_mem, dst_offset, dst_size+dst_offset, 16, stride);
                VP8_CL_FINISH(xd->cl_commands);
                ((int *)(q+q_offset))[0] = 0;
                vp8_cl_mb_prep(xd,QCOEFF);
            }

            q_offset   += 16;
            pre_offset += 4;
            dst_offset += 4;
        }

        pre_offset += 64 - 16;
        dst_offset += 4*stride - 16;
    }

}

void vp8_dequant_idct_add_uv_block_cl(VP8D_COMP *pbi, MACROBLOCKD *xd,
        vp8_dequant_idct_add_uv_block_fn_t idct_add_uv_block
)
{
    int i, j;

    int block_num = 16;
    BLOCKD b = xd->block[block_num];

    short *q = xd->qcoeff;

    cl_mem dst_mem = xd->dst.buffer_mem;
    unsigned char *dst = xd->dst.buffer_alloc;
    int u_off = xd->dst.u_buffer - dst;
    int v_off = xd->dst.v_buffer - dst;

    int stride = xd->dst.uv_stride;
    size_t dst_size = 8*(stride+1);
    char *eobs = xd->eobs+16;

    int pre_offset = block_num*16;
    int q_offset = block_num*16;
    int dst_offset = 0;

    vp8_cl_mb_prep(xd, DIFF|QCOEFF|PREDICTOR);
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            if (*eobs++ > 1){
                vp8_cl_block_prep(&xd->block[0], DEQUANT);
                vp8_dequant_idct_add_cl(&b, dst, dst_mem, u_off+dst_offset, u_off+dst_size, q_offset, pre_offset, 8, stride, DEQUANT_INVOKE (&pbi->dequant, idct_add));
            }
            else
            {
                vp8_cl_block_prep(&xd->block[block_num], DEQUANT);
                vp8_dc_only_idct_add_cl (&b, CL_FALSE, 0, q_offset, pre_offset, dst, dst_mem, u_off+dst_offset, u_off+dst_size, 8, stride);
                
                //Need round trip + finish until qcoeff set in CL
                vp8_cl_block_finish(&xd->block[0], QCOEFF);
                VP8_CL_FINISH(xd->cl_commands);
                ((int *)(q+q_offset))[0] = 0;
                vp8_cl_mb_prep(xd,QCOEFF);
            }

            q_offset    += 16;
            pre_offset  += 4;
            dst_offset += 4;
        }

        pre_offset  += 32 - 8;
        dst_offset += 4*stride - 8;
    }

    //Swap dstu out of cl_mem and dstv into it

    dst_offset = 0;
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            if (*eobs++ > 1){
                vp8_cl_block_prep(&b, DEQUANT);
                vp8_dequant_idct_add_cl (&b, dst, dst_mem, v_off+dst_offset, v_off+dst_size, q_offset,
                        pre_offset, 8, stride, DEQUANT_INVOKE (&pbi->dequant, idct_add));
            }
            else
            {
                vp8_cl_block_prep(&b, DEQUANT);
                vp8_dc_only_idct_add_cl (&b, CL_FALSE, 0, q_offset, pre_offset,
                        dst, dst_mem, v_off+dst_offset, v_off+dst_size, 8, stride);

                //Eventually replace with memset kernel call to prevent round trip
                vp8_cl_mb_finish(xd,QCOEFF);
                VP8_CL_FINISH(xd->cl_commands);
                ((int *)(q+q_offset))[0] = 0;
                vp8_cl_mb_prep(xd,QCOEFF);
            }

            q_offset    += 16;
            pre_offset  += 4;
            dst_offset += 4;
        }

        pre_offset  += 32 - 8;
        dst_offset += 4*stride - 8;
    }
    
    vp8_cl_mb_finish(xd,QCOEFF);

}
