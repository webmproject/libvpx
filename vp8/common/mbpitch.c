/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "blockd.h"

#include "stdio.h"
#include "vpx_config.h"
#if CONFIG_OPENCL
#include "opencl/vp8_opencl.h"
#endif

typedef enum
{
    PRED = 0,
    DEST = 1
} BLOCKSET;

static void setup_block
(
    BLOCKD *b,
    unsigned char **base,
    int Stride,
    int offset,
    BLOCKSET bs
)
{

    if (bs == DEST)
    {
        b->dst_stride = Stride;
        b->dst = offset;
        b->base_dst = base;
    }
    else
    {
        b->pre_stride = Stride;
        b->pre = offset;
        b->base_pre = base;
    }

}


static void setup_macroblock(MACROBLOCKD *x, BLOCKSET bs)
{
    int block;

    unsigned char **y, **u, **v;
    unsigned char **buf_base;
    int y_off, u_off, v_off;

    if (bs == DEST)
    {
        buf_base = &x->dst.buffer_alloc;
        y_off = x->dst.y_buffer - x->dst.buffer_alloc;
        u_off = x->dst.u_buffer - x->dst.buffer_alloc;
        v_off = x->dst.v_buffer - x->dst.buffer_alloc;
        y = &x->dst.y_buffer;
        u = &x->dst.u_buffer;
        v = &x->dst.v_buffer;
        y_off = 0;

        //y = buf_base;
        //y_off = x->dst.y_buffer - x->dst.buffer_alloc;
        
        u = buf_base;
        v = buf_base;

        u_off = x->dst.u_buffer - x->dst.buffer_alloc;
        v_off = x->dst.v_buffer - x->dst.buffer_alloc;

    }
    else
    {
        buf_base = &x->pre.buffer_alloc;
        y = &x->pre.y_buffer;
        u = &x->pre.u_buffer;
        v = &x->pre.v_buffer;
        y_off = u_off = v_off = 0;

        //y = buf_base;
        //y_off = x->pre.y_buffer - x->pre.buffer_alloc;
        //u = buf_base;
        //u_off = x->pre.u_buffer - x->pre.buffer_alloc;
        //v = buf_base;
        //v_off = x->pre.v_buffer - x->pre.buffer_alloc;
    }

    for (block = 0; block < 16; block++) /* y blocks */
    {
        setup_block(&x->block[block], y, x->dst.y_stride,
                        y_off + ((block >> 2) * 4 * x->dst.y_stride + (block & 3) * 4), bs);
    }

    for (block = 16; block < 20; block++) /* U and V blocks */
    {
        int block_off = ((block - 16) >> 1) * 4 * x->dst.uv_stride + (block & 1) * 4;

        setup_block(&x->block[block], u, x->dst.uv_stride,
                        u_off + block_off, bs);

        setup_block(&x->block[block+4], v, x->dst.uv_stride,
                        v_off + block_off, bs);
    }
}

void vp8_setup_block_dptrs(MACROBLOCKD *x)
{
    int r, c;
    unsigned int offset;

#if CONFIG_OPENCL && !ONE_CQ_PER_MB
    cl_command_queue y_cq, u_cq, v_cq;
    int err;
    if (cl_initialized == CL_SUCCESS){
        //Create command queue for Y/U/V Planes
        y_cq = clCreateCommandQueue(cl_data.context, cl_data.device_id, 0, &err);
        if (!y_cq || err != CL_SUCCESS) {
            printf("Error: Failed to create a command queue!\n");
            cl_destroy(NULL, VP8_CL_TRIED_BUT_FAILED);
        }
        u_cq = clCreateCommandQueue(cl_data.context, cl_data.device_id, 0, &err);
        if (!u_cq || err != CL_SUCCESS) {
            printf("Error: Failed to create a command queue!\n");
            cl_destroy(NULL, VP8_CL_TRIED_BUT_FAILED);
        }
        v_cq = clCreateCommandQueue(cl_data.context, cl_data.device_id, 0, &err);
        if (!v_cq || err != CL_SUCCESS) {
            printf("Error: Failed to create a command queue!\n");
            cl_destroy(NULL, VP8_CL_TRIED_BUT_FAILED);
        }
    }
#endif

    /* 16 Y blocks */
    for (r = 0; r < 4; r++)
    {
        for (c = 0; c < 4; c++)
        {
            offset = r * 4 * 16 + c * 4;
            x->block[r*4+c].diff_offset      = offset;
            x->block[r*4+c].predictor_offset = offset;
#if CONFIG_OPENCL && !ONE_CQ_PER_MB
            if (cl_initialized == CL_SUCCESS)
                x->block[r*4+c].cl_commands = y_cq;
#endif
        }
    }

    /* 4 U Blocks */
    for (r = 0; r < 2; r++)
    {
        for (c = 0; c < 2; c++)
        {
            offset = 256 + r * 4 * 8 + c * 4;
            x->block[16+r*2+c].diff_offset      = offset;
            x->block[16+r*2+c].predictor_offset = offset;

#if CONFIG_OPENCL && !ONE_CQ_PER_MB
            if (cl_initialized == CL_SUCCESS)
                x->block[16+r*2+c].cl_commands = u_cq;
#endif
        }
    }

    /* 4 V Blocks */
    for (r = 0; r < 2; r++)
    {
        for (c = 0; c < 2; c++)
        {
            offset = 320+ r * 4 * 8 + c * 4;
            x->block[20+r*2+c].diff_offset      = offset;
            x->block[20+r*2+c].predictor_offset = offset;

#if CONFIG_OPENCL && !ONE_CQ_PER_MB
            if (cl_initialized == CL_SUCCESS)
                x->block[20+r*2+c].cl_commands = v_cq;
#endif
        }
    }

    x->block[24].diff_offset = 384;

    for (r = 0; r < 25; r++)
    {
    	x->block[r].qcoeff_base = x->qcoeff;
    	x->block[r].qcoeff_offset = r * 16;
        x->block[r].dqcoeff_base = x->dqcoeff;
        x->block[r].dqcoeff_offset = r * 16;
        
        x->block[r].predictor_base = x->predictor;
        x->block[r].diff_base = x->diff;
        x->block[r].eobs_base = x->eobs;

#if CONFIG_OPENCL
        if (cl_initialized == CL_SUCCESS){
            /* Copy command queue reference from macroblock */
#if ONE_CQ_PER_MB
            x->block[r].cl_commands = x->cl_commands;
#endif

            /* Set up CL memory buffers as appropriate */
            x->block[r].cl_diff_mem = x->cl_diff_mem;
            x->block[r].cl_dqcoeff_mem = x->cl_dqcoeff_mem;
            x->block[r].cl_eobs_mem = x->cl_eobs_mem;
            x->block[r].cl_predictor_mem = x->cl_predictor_mem;
            x->block[r].cl_qcoeff_mem = x->cl_qcoeff_mem;
        }

        //Copy filter type to block.
        x->block[r].sixtap_filter = x->sixtap_filter;
#endif
    }

}

void vp8_build_block_doffsets(MACROBLOCKD *x)
{
    /* handle the destination pitch features */
    setup_macroblock(x, DEST);
    setup_macroblock(x, PRED);
}
