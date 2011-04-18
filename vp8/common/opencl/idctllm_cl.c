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

//ACW: Remove me after debugging.
#include <stdio.h>
#include <string.h>

#include "idct_cl.h"
#include "idctllm_cl.h"
#include "blockd_cl.h"

void cl_destroy_idct(){

    if (cl_data.idct_program)
        clReleaseProgram(cl_data.idct_program);

    cl_data.idct_program = NULL;
    
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_short_inv_walsh4x4_1_kernel);
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_short_inv_walsh4x4_1st_pass_kernel);
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_short_inv_walsh4x4_2nd_pass_kernel);
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_dc_only_idct_add_kernel);
    //VP8_CL_RELEASE_KERNEL(cl_data.vp8_short_idct4x4llm_1_kernel);
    //VP8_CL_RELEASE_KERNEL(cl_data.vp8_short_idct4x4llm_kernel);

}

int cl_init_idct() {
    int err;

    // Create the filter compute program from the file-defined source code
    if (cl_load_program(&cl_data.idct_program, idctllm_cl_file_name,
            idctCompileOptions) != CL_SUCCESS)
        return VP8_CL_TRIED_BUT_FAILED;

    // Create the compute kernel in the program we wish to run
    VP8_CL_CREATE_KERNEL(cl_data,idct_program,vp8_short_inv_walsh4x4_1_kernel,"vp8_short_inv_walsh4x4_1_kernel");
    VP8_CL_CREATE_KERNEL(cl_data,idct_program,vp8_short_inv_walsh4x4_1st_pass_kernel,"vp8_short_inv_walsh4x4_1st_pass_kernel");
    VP8_CL_CREATE_KERNEL(cl_data,idct_program,vp8_short_inv_walsh4x4_2nd_pass_kernel,"vp8_short_inv_walsh4x4_2nd_pass_kernel");
    VP8_CL_CREATE_KERNEL(cl_data,idct_program,vp8_dc_only_idct_add_kernel,"vp8_dc_only_idct_add_kernel");

    ////idct4x4llm kernels are only useful for the encoder
    //VP8_CL_CREATE_KERNEL(cl_data,idct_program,vp8_short_idct4x4llm_1_kernel,"vp8_short_idct4x4llm_1_kernel");
    //VP8_CL_CREATE_KERNEL(cl_data,idct_program,vp8_short_idct4x4llm_kernel,"vp8_short_idct4x4llm_kernel");

    return CL_SUCCESS;
}

#define max(x,y) (x > y ? x: y)
//#define NO_CL

/* Only useful for encoder... Untested... */
void vp8_short_idct4x4llm_cl(BLOCKD *b, int pitch)
{
    int err;

    short *input = b->dqcoeff_base + b->dqcoeff_offset;
    short *output = &b->diff_base[b->diff_offset];

    cl_mem src_mem, dst_mem;

    //1 instance for now. This should be split into 2-pass * 4 thread.
    size_t global = 1;

    if (cl_initialized != CL_SUCCESS){
        vp8_short_idct4x4llm_c(input,output,pitch);
        return;
    }

    VP8_CL_CREATE_BUF(b->cl_commands, src_mem,,
            sizeof(short)*16, input,
            vp8_short_idct4x4llm_c(input,output,pitch),
    );

    VP8_CL_CREATE_BUF(b->cl_commands, dst_mem,,
            sizeof(short)*(4+(pitch/2)*3), output,
            vp8_short_idct4x4llm_c(input,output,pitch),
    );

    //Set arguments and run kernel
    err = 0;
    err = clSetKernelArg(cl_data.vp8_short_idct4x4llm_kernel, 0, sizeof (cl_mem), &src_mem);
    err |= clSetKernelArg(cl_data.vp8_short_idct4x4llm_kernel, 1, sizeof (cl_mem), &dst_mem);
    err |= clSetKernelArg(cl_data.vp8_short_idct4x4llm_kernel, 2, sizeof (int), &pitch);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",
        vp8_short_idct4x4llm_c(input,output,pitch),
    );
    
    /* Execute the kernel */
    err = clEnqueueNDRangeKernel(b->cl_commands, cl_data.vp8_short_idct4x4llm_kernel, 1, NULL, &global, NULL , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);
        vp8_short_idct4x4llm_c(input,output,pitch),
    );

    /* Read back the result data from the device */
    err = clEnqueueReadBuffer(b->cl_commands, dst_mem, CL_FALSE, 0, sizeof(short)*(4+pitch/2*3), output, 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS(b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to read output array!\n",
        vp8_short_idct4x4llm_c(input,output,pitch),
    );

    clReleaseMemObject(src_mem);
    clReleaseMemObject(dst_mem);

    return;
}

/* Only useful for encoder... Untested... */
void vp8_short_idct4x4llm_1_cl(BLOCKD *b, int pitch)
{
    int err;
    size_t global = 4;

    short *input = b->dqcoeff_base + b->dqcoeff_offset;
    short *output = &b->diff_base[b->diff_offset];

    cl_mem src_mem, dst_mem;

    if (cl_initialized != CL_SUCCESS){
        vp8_short_idct4x4llm_1_c(input,output,pitch);
        return;
    }

    printf("vp8_short_idct4x4llm_1_cl\n");

    VP8_CL_CREATE_BUF(b->cl_commands, src_mem,,
            sizeof(short), input,
            vp8_short_idct4x4llm_1_c(input,output,pitch),
    );

    VP8_CL_CREATE_BUF(b->cl_commands, dst_mem,,
            sizeof(short)*(4+(pitch/2)*3), output,
            vp8_short_idct4x4llm_1_c(input,output,pitch),
    );

    //Set arguments and run kernel
    err = 0;
    err = clSetKernelArg(cl_data.vp8_short_idct4x4llm_1_kernel, 0, sizeof (cl_mem), &src_mem);
    err |= clSetKernelArg(cl_data.vp8_short_idct4x4llm_1_kernel, 1, sizeof (cl_mem), &dst_mem);
    err |= clSetKernelArg(cl_data.vp8_short_idct4x4llm_1_kernel, 2, sizeof (int), &pitch);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",
        vp8_short_idct4x4llm_1_c(input,output,pitch),
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel(b->cl_commands, cl_data.vp8_short_idct4x4llm_1_kernel, 1, NULL, &global, NULL , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);
        vp8_short_idct4x4llm_1_c(input,output,pitch),
    );

    /* Read back the result data from the device */
    err = clEnqueueReadBuffer(b->cl_commands, dst_mem, CL_FALSE, 0, sizeof(short)*(4+pitch/2*3), output, 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS(b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to read output array!\n",
        vp8_short_idct4x4llm_1_c(input,output,pitch),
    );

    clReleaseMemObject(src_mem);
    clReleaseMemObject(dst_mem);

    return;

}

void vp8_dc_only_idct_add_cl(BLOCKD *b, cl_int use_diff, int diff_offset, 
        int qcoeff_offset, int pred_offset,
        unsigned char *dst_base, cl_mem dst_mem, int dst_offset, size_t dest_size,
        int pitch, int stride
)
{
    
    int err;
    size_t global = 16;

    int free_mem = 0;
    //cl_mem dest_mem = NULL;

    if (dst_mem == NULL){
        VP8_CL_CREATE_BUF(b->cl_commands, dst_mem,,
                dest_size, dst_base,,
        );
        free_mem = 1;
    }

    //Set arguments and run kernel
    err =  clSetKernelArg(cl_data.vp8_dc_only_idct_add_kernel, 0, sizeof (cl_mem), &b->cl_predictor_mem);
    err |= clSetKernelArg(cl_data.vp8_dc_only_idct_add_kernel, 1, sizeof (int), &pred_offset);
    err |= clSetKernelArg(cl_data.vp8_dc_only_idct_add_kernel, 2, sizeof (cl_mem), &dst_mem);
    err |= clSetKernelArg(cl_data.vp8_dc_only_idct_add_kernel, 3, sizeof (int), &dst_offset);
    err |= clSetKernelArg(cl_data.vp8_dc_only_idct_add_kernel, 4, sizeof (int), &pitch);
    err |= clSetKernelArg(cl_data.vp8_dc_only_idct_add_kernel, 5, sizeof (int), &stride);
    err |= clSetKernelArg(cl_data.vp8_dc_only_idct_add_kernel, 6, sizeof (cl_int), &use_diff);
    err |= clSetKernelArg(cl_data.vp8_dc_only_idct_add_kernel, 7, sizeof (cl_mem), &b->cl_diff_mem);
    err |= clSetKernelArg(cl_data.vp8_dc_only_idct_add_kernel, 8, sizeof (int), &diff_offset);
    err |= clSetKernelArg(cl_data.vp8_dc_only_idct_add_kernel, 9, sizeof (cl_mem), &b->cl_qcoeff_mem);
    err |= clSetKernelArg(cl_data.vp8_dc_only_idct_add_kernel, 10, sizeof (int), &qcoeff_offset);
    err |= clSetKernelArg(cl_data.vp8_dc_only_idct_add_kernel, 11, sizeof (cl_mem), &b->cl_dequant_mem);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",,
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel(b->cl_commands, cl_data.vp8_dc_only_idct_add_kernel, 1, NULL, &global, NULL , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);,
    );


    if (free_mem == 1){
    /* Read back the result data from the device */
        err = clEnqueueReadBuffer(b->cl_commands, dst_mem, CL_FALSE, 0,
                dest_size, dst_base, 0, NULL, NULL);

        VP8_CL_CHECK_SUCCESS(b->cl_commands, err != CL_SUCCESS,
            "Error: Failed to read output array!\n",,
        );

        clReleaseMemObject(dst_mem);
    }

    return;
}

void vp8_short_inv_walsh4x4_cl(BLOCKD *b)
{
    int err;
    size_t global = 4;

    if (cl_initialized != CL_SUCCESS){
        vp8_short_inv_walsh4x4_c(b->dqcoeff_base+b->dqcoeff_offset,&b->diff_base[b->diff_offset]);
        return;
    }

    //Set arguments and run kernel
    err = 0;
    err = clSetKernelArg(cl_data.vp8_short_inv_walsh4x4_1st_pass_kernel, 0, sizeof (cl_mem), &b->cl_dqcoeff_mem);
    err |= clSetKernelArg(cl_data.vp8_short_inv_walsh4x4_1st_pass_kernel, 1, sizeof(int), &b->dqcoeff_offset);
    err |= clSetKernelArg(cl_data.vp8_short_inv_walsh4x4_1st_pass_kernel, 2, sizeof (cl_mem), &b->cl_diff_mem);
    err |= clSetKernelArg(cl_data.vp8_short_inv_walsh4x4_1st_pass_kernel, 3, sizeof(int), &b->diff_offset);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",
        vp8_short_inv_walsh4x4_c(b->dqcoeff_base+b->dqcoeff_offset, &b->diff_base[b->diff_offset]),
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel(b->cl_commands, cl_data.vp8_short_inv_walsh4x4_1st_pass_kernel, 1, NULL, &global, NULL , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);
        vp8_short_inv_walsh4x4_c(b->dqcoeff_base+b->dqcoeff_offset, &b->diff_base[b->diff_offset]),
    );

    //Second pass
    //Set arguments and run kernel
    err = 0;
    err = clSetKernelArg(cl_data.vp8_short_inv_walsh4x4_2nd_pass_kernel, 0, sizeof (cl_mem), &b->cl_diff_mem);
    err |= clSetKernelArg(cl_data.vp8_short_inv_walsh4x4_2nd_pass_kernel, 1, sizeof(int), &b->diff_offset);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",
        vp8_short_inv_walsh4x4_c(b->dqcoeff_base+b->dqcoeff_offset, &b->diff_base[b->diff_offset]),
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel(b->cl_commands, cl_data.vp8_short_inv_walsh4x4_2nd_pass_kernel, 1, NULL, &global, NULL , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);
        vp8_short_inv_walsh4x4_c(b->dqcoeff_base+b->dqcoeff_offset, &b->diff_base[b->diff_offset]),
    );

    return;
}

void vp8_short_inv_walsh4x4_1_cl(BLOCKD *b)
{
    
    int err;
    size_t global = 4;

    if (cl_initialized != CL_SUCCESS){
        vp8_short_inv_walsh4x4_1_c(b->dqcoeff_base + b->dqcoeff_offset,
            &b->diff_base[b->diff_offset]);
        return;
    }

    //Set arguments and run kernel
    err = 0;
    err = clSetKernelArg(cl_data.vp8_short_inv_walsh4x4_1_kernel, 0, sizeof (cl_mem), &b->cl_dqcoeff_mem);
    err |= clSetKernelArg(cl_data.vp8_short_inv_walsh4x4_1_kernel, 1, sizeof (int), &b->dqcoeff_offset);
    err |= clSetKernelArg(cl_data.vp8_short_inv_walsh4x4_1_kernel, 2, sizeof (cl_mem), &b->cl_diff_mem);
    err |= clSetKernelArg(cl_data.vp8_short_inv_walsh4x4_1_kernel, 3, sizeof (int), &b->diff_offset);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",
        vp8_short_inv_walsh4x4_1_c(b->dqcoeff_base + b->dqcoeff_offset,
            &b->diff_base[b->diff_offset]),
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel(b->cl_commands, cl_data.vp8_short_inv_walsh4x4_1_kernel, 1, NULL, &global, NULL , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( b->cl_commands, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);
        vp8_short_inv_walsh4x4_1_c(b->dqcoeff_base + b->dqcoeff_offset,
                &b->diff_base[b->diff_offset]),
    );

    return;
}
