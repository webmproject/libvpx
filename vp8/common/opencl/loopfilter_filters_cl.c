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

#include <stdio.h>

#include "vpx_ports/config.h"
#include "vp8_opencl.h"
#include "blockd_cl.h"

//#include "loopfilter_cl.h"
//#include "../onyxc_int.h"

typedef unsigned char uc;

static void vp8_loop_filter_cl_run(
    cl_command_queue cq,
    cl_kernel kernel,
    cl_mem buf_mem,
    int s_off,
    int p,
    const signed char *flimit,
    const signed char *limit,
    const signed char *thresh,
    int count,
    int block_cnt
){
    size_t global[] = {count,block_cnt};
    int err;

    cl_mem flimit_mem;
    cl_mem limit_mem;
    cl_mem thresh_mem;

    VP8_CL_CREATE_BUF(cq, flimit_mem, , sizeof(uc)*16, flimit,, );
    VP8_CL_CREATE_BUF(cq, limit_mem, , sizeof(uc)*16, limit,, );
    VP8_CL_CREATE_BUF(cq, thresh_mem, , sizeof(uc)*16, thresh,, );

    err = 0;
    err = clSetKernelArg(kernel, 0, sizeof (cl_mem), &buf_mem);
    err |= clSetKernelArg(kernel, 1, sizeof (cl_int), &s_off);
    err |= clSetKernelArg(kernel, 2, sizeof (cl_int), &p);
    err |= clSetKernelArg(kernel, 3, sizeof (cl_mem), &flimit_mem);
    err |= clSetKernelArg(kernel, 4, sizeof (cl_mem), &limit_mem);
    err |= clSetKernelArg(kernel, 5, sizeof (cl_mem), &thresh_mem);
    err |= clSetKernelArg(kernel, 6, sizeof (cl_int), &block_cnt);
    VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",,
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel(cq, kernel, 2, NULL, global, NULL , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);,
    );

    clReleaseMemObject(flimit_mem);
    clReleaseMemObject(limit_mem);
    clReleaseMemObject(thresh_mem);

    VP8_CL_FINISH(cq);
}

void vp8_loop_filter_horizontal_edge_cl
(
    MACROBLOCKD *x,
    cl_mem s_base,
    int s_off,
    int p, /* pitch */
    const signed char *flimit,
    const signed char *limit,
    const signed char *thresh,
    int count,
    int block_cnt
)
{
    vp8_loop_filter_cl_run(x->cl_commands,
        cl_data.vp8_loop_filter_horizontal_edge_kernel, s_base, s_off,
        p, flimit, limit, thresh, count*8, block_cnt
    );
}

void vp8_loop_filter_vertical_edge_cl
(
    MACROBLOCKD *x,
    cl_mem s_base,
    int s_off,
    int p,
    const signed char *flimit,
    const signed char *limit,
    const signed char *thresh,
    int count,
    int block_cnt
)
{
    vp8_loop_filter_cl_run(x->cl_commands,
        cl_data.vp8_loop_filter_vertical_edge_kernel, s_base, s_off,
        p, flimit, limit, thresh, count*8, block_cnt
    );
}

void vp8_mbloop_filter_horizontal_edge_cl
(
    MACROBLOCKD *x,
    cl_mem s_base,
    int s_off,
    int p,
    const signed char *flimit,
    const signed char *limit,
    const signed char *thresh,
    int count,
    int block_cnt
)
{
    vp8_loop_filter_cl_run(x->cl_commands,
        cl_data.vp8_mbloop_filter_horizontal_edge_kernel, s_base, s_off,
        p, flimit, limit, thresh, count*8, block_cnt
    );
}


void vp8_mbloop_filter_vertical_edge_cl
(
    MACROBLOCKD *x,
    cl_mem s_base,
    int s_off,
    int p,
    const signed char *flimit,
    const signed char *limit,
    const signed char *thresh,
    int count,
    int block_cnt
)
{
    vp8_loop_filter_cl_run(x->cl_commands,
        cl_data.vp8_mbloop_filter_vertical_edge_kernel, s_base, s_off,
        p, flimit, limit, thresh, count*8, block_cnt
    );
}

void vp8_loop_filter_simple_horizontal_edge_cl
(
    MACROBLOCKD *x,
    cl_mem s_base,
    int s_off,
    int p,
    const signed char *flimit,
    const signed char *limit,
    const signed char *thresh,
    int count,
    int block_cnt
)
{
    vp8_loop_filter_cl_run(x->cl_commands,
        cl_data.vp8_loop_filter_simple_horizontal_edge_kernel, s_base, s_off,
        p, flimit, limit, thresh, count*8, block_cnt
    );
}

void vp8_loop_filter_simple_vertical_edge_cl
(
    MACROBLOCKD *x,
    cl_mem s_base,
    int s_off,
    int p,
    const signed char *flimit,
    const signed char *limit,
    const signed char *thresh,
    int count,
    int block_cnt
)
{
    vp8_loop_filter_cl_run(x->cl_commands,
        cl_data.vp8_loop_filter_simple_vertical_edge_kernel, s_base, s_off,
        p, flimit, limit, thresh, count*8, block_cnt
    );
}
