/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

//ACW: Remove me after debugging.
#include <stdio.h>
#include <string.h>

#include "vp8/common/opencl/blockd_cl.h"
#include "vp8/common/opencl/idct_cl.h"
#include "dequantize_cl.h"

const char *dequantCompileOptions = "";
const char *dequant_cl_file_name = "vp8/decoder/opencl/dequantize_cl.cl";

void cl_memset_short(short *s, int c, size_t n) {
    for (n /= sizeof(short); n > 0; --n)
        *s++ = c;
}

void vp8_memset_short_cl(cl_mem mem, int offset, short val){

}

int cl_destroy_dequant(){
    printf("Freeing dequant decoder resources\n");

    VP8_CL_RELEASE_KERNEL(cl_data.vp8_dequant_dc_idct_add_kernel);
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_dequant_idct_add_kernel);
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_dequantize_b_kernel);

    if (cl_data.dequant_program)
        clReleaseProgram(cl_data.dequant_program);
    cl_data.dequant_program = NULL;

    return CL_SUCCESS;
}

int cl_init_dequant() {
    int err;

    //printf("Initializing dequant program/kernels\n");

    // Create the compute program from the file-defined source code
    if (cl_load_program(&cl_data.dequant_program, dequant_cl_file_name,
            dequantCompileOptions) != CL_SUCCESS)
        return VP8_CL_TRIED_BUT_FAILED;

    // Create the compute kernels in the program we wish to run
    VP8_CL_CREATE_KERNEL(cl_data,dequant_program,vp8_dequant_dc_idct_add_kernel,"vp8_dequant_dc_idct_add_kernel");
    VP8_CL_CREATE_KERNEL(cl_data,dequant_program,vp8_dequant_idct_add_kernel,"vp8_dequant_idct_add_kernel");
    VP8_CL_CREATE_KERNEL(cl_data,dequant_program,vp8_dequantize_b_kernel,"vp8_dequantize_b_kernel");

    //printf("Created dequant kernels\n");

    return CL_SUCCESS;
}

void vp8_dequantize_b_cl(BLOCKD *d)
{
    int err;
    size_t global = 16;

    /* Set kernel arguments */
    err = 0;
    err = clSetKernelArg(cl_data.vp8_dequantize_b_kernel, 0, sizeof (cl_mem), &d->cl_dqcoeff_mem);
    err |= clSetKernelArg(cl_data.vp8_dequantize_b_kernel, 1, sizeof (cl_int), &d->dqcoeff_offset);
    err |= clSetKernelArg(cl_data.vp8_dequantize_b_kernel, 2, sizeof (cl_mem), &d->cl_qcoeff_mem);
    err |= clSetKernelArg(cl_data.vp8_dequantize_b_kernel, 3, sizeof (cl_int), &d->qcoeff_offset);
    err |= clSetKernelArg(cl_data.vp8_dequantize_b_kernel, 4, sizeof (cl_mem), &d->cl_dequant_mem);
    VP8_CL_CHECK_SUCCESS( d->cl_commands, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",
        vp8_dequantize_b_c(d),
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel( d->cl_commands, cl_data.vp8_dequantize_b_kernel, 1, NULL, &global, NULL , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( d->cl_commands, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);\
        vp8_dequantize_b_c(d),
    );

}

void vp8_dequant_idct_add_cl(BLOCKD *b, unsigned char *dest_base, cl_mem dest_mem, int dest_offset, size_t dst_size, int q_offset, int pred_offset, int pitch, int stride, vp8_dequant_idct_add_fn_t idct_add)
{
    int err;
    size_t global = 1;
    //cl_mem dest_mem = NULL;
    int free_mem = 0;

    if (dest_mem == NULL){
        //Initialize destination memory
        VP8_CL_CREATE_BUF(b->cl_commands, dest_mem, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
                dst_size, dest_base,,
        );
        free_mem = 1;
    }
    
    /* Set kernel arguments */
    err = 0;
    err = clSetKernelArg(cl_data.vp8_dequant_idct_add_kernel, 0, sizeof (cl_mem), &b->cl_qcoeff_mem);
    err |= clSetKernelArg(cl_data.vp8_dequant_idct_add_kernel, 1, sizeof (int), &q_offset);
    err |= clSetKernelArg(cl_data.vp8_dequant_idct_add_kernel, 2, sizeof (cl_mem), &b->cl_dequant_mem);
    err |= clSetKernelArg(cl_data.vp8_dequant_idct_add_kernel, 3, sizeof (cl_mem), &b->cl_predictor_mem);
    err |= clSetKernelArg(cl_data.vp8_dequant_idct_add_kernel, 4, sizeof (int), &pred_offset);
    err |= clSetKernelArg(cl_data.vp8_dequant_idct_add_kernel, 5, sizeof (cl_mem), &dest_mem);
    err |= clSetKernelArg(cl_data.vp8_dequant_idct_add_kernel, 6, sizeof (int), &dest_offset);
    err |= clSetKernelArg(cl_data.vp8_dequant_idct_add_kernel, 7, sizeof (int), &pitch);
    err |= clSetKernelArg(cl_data.vp8_dequant_idct_add_kernel, 8, sizeof (int), &stride);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",
        idct_add(b->qcoeff_base+q_offset, b->dequant,  b->predictor_base + pred_offset,
            dest_base + dest_offset, pitch, stride),
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel( b->cl_commands, cl_data.vp8_dequant_idct_add_kernel, 1, NULL, &global, NULL , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);\
        idct_add(b->qcoeff_base+q_offset, b->dequant,  b->predictor_base + pred_offset,
            dest_base + dest_offset, pitch, stride),
    );

    if (free_mem == 1){
        /* Read back the result data from the device */
        err = clEnqueueReadBuffer(b->cl_commands, dest_mem, CL_FALSE, 0, dst_size, dest_base, 0, NULL, NULL);
        VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
            "Error: Failed to read output array!\n",
            idct_add(b->qcoeff_base+q_offset, b->dequant,  b->predictor_base + pred_offset,
                dest_base + dest_offset, pitch, stride),
        );

        //CL Spec says this can be freed without clFinish first
        clReleaseMemObject(dest_mem);
    }

    return;
}

//Can modify arguments. Only called from vp8_dequant_dc_idct_add_y_block_cl.
void vp8_dequant_dc_idct_add_cl(
    BLOCKD *b,
    int qcoeff_offset,
    int pred_offset,
    unsigned char *dest_base,
    int dest_off,
    int pitch,
    int stride,
    int Dc_offset)
{
    int err;
    int dq_offset = 0;
    unsigned char *dest = dest_base + dest_off;
 
    cl_mem dest_mem = NULL;
    size_t dest_size;
    size_t global = 1;
    int dest_offset=0;

    //Initialize dest_mem
    dest_size = sizeof(cl_uchar)*(4*stride + dest_offset + 4);
    VP8_CL_CREATE_BUF(b->cl_commands, dest_mem, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
            dest_size, dest,,
    );

    //Assuming that all input cl_mem has been initialized outside of this Fn.

    /* Set kernel arguments */
    err = 0;
    err = clSetKernelArg(cl_data.vp8_dequant_dc_idct_add_kernel, 0, sizeof (cl_mem), &b->cl_qcoeff_mem);
    err |= clSetKernelArg(cl_data.vp8_dequant_dc_idct_add_kernel, 1, sizeof (int), &qcoeff_offset);
    err |= clSetKernelArg(cl_data.vp8_dequant_dc_idct_add_kernel, 2, sizeof (cl_mem), &b->cl_dequant_mem);
    err |= clSetKernelArg(cl_data.vp8_dequant_dc_idct_add_kernel, 3, sizeof(int), &dq_offset);
    err |= clSetKernelArg(cl_data.vp8_dequant_dc_idct_add_kernel, 4, sizeof (cl_mem), &b->cl_predictor_mem);
    err |= clSetKernelArg(cl_data.vp8_dequant_dc_idct_add_kernel, 5, sizeof (int), &pred_offset);
    err |= clSetKernelArg(cl_data.vp8_dequant_dc_idct_add_kernel, 6, sizeof (cl_mem), &b->cl_diff_mem);
    err |= clSetKernelArg(cl_data.vp8_dequant_dc_idct_add_kernel, 7, sizeof (int), &Dc_offset);
    err |= clSetKernelArg(cl_data.vp8_dequant_dc_idct_add_kernel, 8, sizeof (cl_mem), &dest_mem);
    err |= clSetKernelArg(cl_data.vp8_dequant_dc_idct_add_kernel, 9, sizeof (int), &pitch);
    err |= clSetKernelArg(cl_data.vp8_dequant_dc_idct_add_kernel, 10, sizeof (int), &stride);

    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",,
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel( b->cl_commands, cl_data.vp8_dequant_dc_idct_add_kernel, 1, NULL, &global, NULL , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);,
    );

    /* Read back the result data from the device */
    err = clEnqueueReadBuffer(b->cl_commands, dest_mem, CL_FALSE, 0, dest_size, dest, 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to read output array!\n",,
    );

    //CL Spec says this can be freed without clFinish first
    clReleaseMemObject(dest_mem);
    dest_mem = NULL;

    return;

}
