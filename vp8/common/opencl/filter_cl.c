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

#include "vp8_opencl.h"
#include "filter_cl.h"
#include "../blockd.h"

#define SIXTAP_FILTER_LEN 6

const char *filterCompileOptions = "-Ivp8/common/opencl -DVP8_FILTER_WEIGHT=128 -DVP8_FILTER_SHIFT=7 -DFILTER_OFFSET";
const char *filter_cl_file_name = "vp8/common/opencl/filter_cl.cl";

#define STATIC_MEM 1
#if STATIC_MEM
static cl_mem int_mem = NULL;
#endif

void cl_destroy_filter(){

    if (cl_data.filter_program)
        clReleaseProgram(cl_data.filter_program);

    //VP8_CL_RELEASE_KERNEL(cl_data.vp8_block_variation_kernel);
#if !TWO_PASS_SIXTAP
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_sixtap_predict_kernel);
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_sixtap_predict8x8_kernel);
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_sixtap_predict8x4_kernel);
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_sixtap_predict16x16_kernel);
#else
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_filter_block2d_first_pass_kernel);
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_filter_block2d_second_pass_kernel);
#endif
    //VP8_CL_RELEASE_KERNEL(cl_data.vp8_bilinear_predict4x4_kernel);
    //VP8_CL_RELEASE_KERNEL(cl_data.vp8_bilinear_predict8x4_kernel);
    //VP8_CL_RELEASE_KERNEL(cl_data.vp8_bilinear_predict8x8_kernel);
    //VP8_CL_RELEASE_KERNEL(cl_data.vp8_bilinear_predict16x16_kernel);

#if MEM_COPY_KERNEL
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_memcpy_kernel);
#endif

    VP8_CL_RELEASE_KERNEL(cl_data.vp8_filter_block2d_bil_first_pass_kernel);
    VP8_CL_RELEASE_KERNEL(cl_data.vp8_filter_block2d_bil_second_pass_kernel);

#if STATIC_MEM
    if (int_mem != NULL)
        clReleaseMemObject(int_mem);
    int_mem = NULL;
#endif

    cl_data.filter_program = NULL;
}

int cl_init_filter() {
    int err;


    // Create the filter compute program from the file-defined source code
    if ( cl_load_program(&cl_data.filter_program, filter_cl_file_name,
            filterCompileOptions) != CL_SUCCESS )
        return VP8_CL_TRIED_BUT_FAILED;

    // Create the compute kernel in the program we wish to run
#if TWO_PASS_SIXTAP
    VP8_CL_CREATE_KERNEL(cl_data,filter_program,vp8_filter_block2d_first_pass_kernel,"vp8_filter_block2d_first_pass_kernel");
    VP8_CL_CREATE_KERNEL(cl_data,filter_program,vp8_filter_block2d_second_pass_kernel,"vp8_filter_block2d_second_pass_kernel");
    VP8_CL_CALC_LOCAL_SIZE(vp8_filter_block2d_first_pass_kernel,vp8_filter_block2d_first_pass_kernel_size);
    VP8_CL_CALC_LOCAL_SIZE(vp8_filter_block2d_second_pass_kernel,vp8_filter_block2d_second_pass_kernel_size);
#else
    VP8_CL_CREATE_KERNEL(cl_data,filter_program,vp8_sixtap_predict_kernel,"vp8_sixtap_predict_kernel");
    VP8_CL_CALC_LOCAL_SIZE(vp8_sixtap_predict_kernel,vp8_sixtap_predict_kernel_size);
    VP8_CL_CREATE_KERNEL(cl_data,filter_program,vp8_sixtap_predict8x8_kernel,"vp8_sixtap_predict8x8_kernel");
    VP8_CL_CALC_LOCAL_SIZE(vp8_sixtap_predict8x8_kernel,vp8_sixtap_predict8x8_kernel_size);
    VP8_CL_CREATE_KERNEL(cl_data,filter_program,vp8_sixtap_predict8x4_kernel,"vp8_sixtap_predict8x4_kernel");
    VP8_CL_CALC_LOCAL_SIZE(vp8_sixtap_predict8x4_kernel,vp8_sixtap_predict8x4_kernel_size);
    VP8_CL_CREATE_KERNEL(cl_data,filter_program,vp8_sixtap_predict16x16_kernel,"vp8_sixtap_predict16x16_kernel");
    VP8_CL_CALC_LOCAL_SIZE(vp8_sixtap_predict16x16_kernel,vp8_sixtap_predict16x16_kernel_size);
#endif
    
    //VP8_CL_CALC_LOCAL_SIZE(vp8_filter_block2d_bil_first_pass_kernel,vp8_filter_block2d_bil_first_pass_kernel_size);
    //VP8_CL_CALC_LOCAL_SIZE(vp8_filter_block2d_bil_second_pass_kernel,vp8_filter_block2d_bil_second_pass_kernel_size);
    VP8_CL_CREATE_KERNEL(cl_data,filter_program,vp8_filter_block2d_bil_first_pass_kernel,"vp8_filter_block2d_bil_first_pass_kernel");
    VP8_CL_CREATE_KERNEL(cl_data,filter_program,vp8_filter_block2d_bil_second_pass_kernel,"vp8_filter_block2d_bil_second_pass_kernel");


    //VP8_CL_CREATE_KERNEL(cl_data,filter_program,vp8_bilinear_predict4x4_kernel,"vp8_bilinear_predict4x4_kernel");
    //VP8_CL_CREATE_KERNEL(cl_data,filter_program,vp8_bilinear_predict8x4_kernel,"vp8_bilinear_predict8x4_kernel");
    //VP8_CL_CREATE_KERNEL(cl_data,filter_program,vp8_bilinear_predict8x8_kernel,"vp8_bilinear_predict8x8_kernel");
    //VP8_CL_CREATE_KERNEL(cl_data,filter_program,vp8_bilinear_predict16x16_kernel,"vp8_bilinear_predict16x16_kernel");

#if MEM_COPY_KERNEL
    VP8_CL_CREATE_KERNEL(cl_data,filter_program,vp8_memcpy_kernel,"vp8_memcpy_kernel");
    VP8_CL_CALC_LOCAL_SIZE(vp8_memcpy_kernel,vp8_memcpy_kernel_size);
#endif

#if STATIC_MEM
    VP8_CL_CREATE_BUF(NULL, int_mem, NULL, sizeof(cl_int)*21*16, NULL, ,err);
#endif

    return CL_SUCCESS;
}

void vp8_filter_block2d_first_pass_cl(
    cl_command_queue cq,
    cl_mem src_mem,
    int src_offset,
    cl_mem int_mem,
    unsigned int src_pixels_per_line,
    unsigned int int_height,
    unsigned int int_width,
    int xoffset
){
    int err;
    size_t global = int_width*int_height;
    size_t local = cl_data.vp8_filter_block2d_first_pass_kernel_size;
    if (local > global)
        local = global;

    err =  clSetKernelArg(cl_data.vp8_filter_block2d_first_pass_kernel, 0, sizeof (cl_mem), &src_mem);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_first_pass_kernel, 1, sizeof (int), &src_offset);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_first_pass_kernel, 2, sizeof (cl_mem), &int_mem);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_first_pass_kernel, 3, sizeof (cl_uint), &src_pixels_per_line);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_first_pass_kernel, 4, sizeof (cl_uint), &int_height);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_first_pass_kernel, 5, sizeof (cl_int), &int_width);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_first_pass_kernel, 6, sizeof (int), &xoffset);
    VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",
        ,
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel( cq, cl_data.vp8_filter_block2d_first_pass_kernel, 1, NULL, &global, &local , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);,
    );
}

void vp8_filter_block2d_second_pass_cl(
    cl_command_queue cq,
    cl_mem int_mem,
    int int_offset,
    cl_mem dst_mem,
    int dst_offset,
    int dst_pitch,
    unsigned int output_height,
    unsigned int output_width,
    int yoffset
){
    int err;
    size_t global = output_width*output_height;
    size_t local = cl_data.vp8_filter_block2d_second_pass_kernel_size;
    if (local > global){
        //printf("Local is now %ld\n",global);
        local = global;
    }

    /* Set kernel arguments */
    err =  clSetKernelArg(cl_data.vp8_filter_block2d_second_pass_kernel, 0, sizeof (cl_mem), &int_mem);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_second_pass_kernel, 1, sizeof (int), &int_offset);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_second_pass_kernel, 2, sizeof (cl_mem), &dst_mem);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_second_pass_kernel, 3, sizeof (int), &dst_offset);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_second_pass_kernel, 4, sizeof (int), &dst_pitch);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_second_pass_kernel, 5, sizeof (int), &output_width);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_second_pass_kernel, 6, sizeof (int), &output_width);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_second_pass_kernel, 7, sizeof (int), &output_height);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_second_pass_kernel, 8, sizeof (int), &output_width);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_second_pass_kernel, 9, sizeof (int), &yoffset);
    VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",
        ,
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel( cq, cl_data.vp8_filter_block2d_second_pass_kernel, 1, NULL, &global, &local , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);,
    );
}

void vp8_sixtap_single_pass(
    cl_command_queue cq,
    cl_kernel kernel,
    size_t local,
    size_t global,
    cl_mem src_mem,
    cl_mem dst_mem,
    unsigned char *src_base,
    int src_offset,
    size_t src_len,
    int src_pixels_per_line,
    int xoffset,
    int yoffset,
    unsigned char *dst_base,
    int dst_offset,
    int dst_pitch,
    size_t dst_len
){
    int err;

#if !STATIC_MEM
    cl_mem int_mem;
#endif

    int free_src = 0, free_dst = 0;

    if (local > global){
        local = global;
    }

    /* Make space for kernel input/output data.
     * Initialize the buffer as well if needed.
     */
    if (src_mem == NULL){
        VP8_CL_CREATE_BUF( cq, src_mem,, sizeof (unsigned char) * src_len, src_base-2,,);
        src_offset = 2;
        free_src = 1;
    } else {
        src_offset -= 2*src_pixels_per_line;
    }

    if (dst_mem == NULL){
        VP8_CL_CREATE_BUF( cq, dst_mem,, sizeof (unsigned char) * dst_len + dst_offset, dst_base,, );
        free_dst = 1;
    }

#if !STATIC_MEM
    CL_CREATE_BUF( cq, int_mem,, sizeof(cl_int)*FData_height*FData_width, NULL,, );
#endif

    err =  clSetKernelArg(kernel, 0, sizeof (cl_mem), &src_mem);
    err |= clSetKernelArg(kernel, 1, sizeof (int), &src_offset);
    err |= clSetKernelArg(kernel, 2, sizeof (cl_int), &src_pixels_per_line);
    err |= clSetKernelArg(kernel, 3, sizeof (cl_int), &xoffset);
    err |= clSetKernelArg(kernel, 4, sizeof (cl_int), &yoffset);
    err |= clSetKernelArg(kernel, 5, sizeof (cl_mem), &dst_mem);
    err |= clSetKernelArg(kernel, 6, sizeof (cl_int), &dst_offset);
    err |= clSetKernelArg(kernel, 7, sizeof (int), &dst_pitch);
    VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",
        ,
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel( cq, kernel, 1, NULL, &global, &local , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);,
    );

    if (free_src == 1)
        clReleaseMemObject(src_mem);

    if (free_dst == 1){
        /* Read back the result data from the device */
        err = clEnqueueReadBuffer(cq, dst_mem, CL_FALSE, 0, sizeof (unsigned char) * dst_len + dst_offset, dst_base, 0, NULL, NULL);
        VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
            "Error: Failed to read output array!\n",
            ,
        );
        clReleaseMemObject(dst_mem);
    }
}

void vp8_sixtap_run_cl(
    cl_command_queue cq,
    cl_mem src_mem,
    cl_mem dst_mem,
    unsigned char *src_base,
    int src_offset,
    size_t src_len,
    int src_pixels_per_line,
    int xoffset,
    int yoffset,
    unsigned char *dst_base,
    int dst_offset,
    int dst_pitch,
    size_t dst_len,
    unsigned int FData_height,
    unsigned int FData_width,
    unsigned int output_height,
    unsigned int output_width,
    int int_offset
)
{
    int err;

#if !STATIC_MEM
    cl_mem int_mem;
#endif

    int free_src = 0, free_dst = 0;

    /* Make space for kernel input/output data.
     * Initialize the buffer as well if needed.
     */
    if (src_mem == NULL){
        VP8_CL_CREATE_BUF( cq, src_mem,, sizeof (unsigned char) * src_len, src_base-2,,);
        src_offset = 2;
        free_src = 1;
    } else {
        src_offset -= 2*src_pixels_per_line;
    }

    if (dst_mem == NULL){
        VP8_CL_CREATE_BUF( cq, dst_mem,, sizeof (unsigned char) * dst_len + dst_offset, dst_base,, );
        free_dst = 1;
    }

#if !STATIC_MEM
    CL_CREATE_BUF( cq, int_mem,, sizeof(cl_int)*FData_height*FData_width, NULL,, );
#endif

    vp8_filter_block2d_first_pass_cl(
        cq, src_mem, src_offset, int_mem, src_pixels_per_line,
        FData_height, FData_width, xoffset
    );

    vp8_filter_block2d_second_pass_cl(cq,int_mem,int_offset,dst_mem,dst_offset,dst_pitch,
            output_height,output_width,yoffset);

    if (free_src == 1)
        clReleaseMemObject(src_mem);

    if (free_dst == 1){
        /* Read back the result data from the device */
        err = clEnqueueReadBuffer(cq, dst_mem, CL_FALSE, 0, sizeof (unsigned char) * dst_len + dst_offset, dst_base, 0, NULL, NULL);
        VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
            "Error: Failed to read output array!\n",
            ,
        );
        clReleaseMemObject(dst_mem);
    }

#if !STATIC_MEM
    clReleaseMemObject(int_mem);
#endif
}

void vp8_sixtap_predict4x4_cl
(
    cl_command_queue cq,
    unsigned char *src_base,
    cl_mem src_mem,
    int src_offset,
    int src_pixels_per_line,
    int xoffset,
    int yoffset,
    unsigned char *dst_base,
    cl_mem dst_mem,
    int dst_offset,
    int dst_pitch
) {

    int output_width=4, output_height=4, FData_height=9, FData_width=4;

    //Size of output to transfer
    int dst_len = DST_LEN(dst_pitch,output_height,output_width);
    int src_len = SIXTAP_SRC_LEN(FData_width,FData_height,src_pixels_per_line);

#if TWO_PASS_SIXTAP
    int int_offset = 8;
    unsigned char *src_ptr = src_base + src_offset;

    vp8_sixtap_run_cl(cq, src_mem, dst_mem,
            (src_ptr-2*src_pixels_per_line),src_offset, src_len,
            src_pixels_per_line, xoffset,yoffset,dst_base,dst_offset,
            dst_pitch,dst_len,FData_height,FData_width,output_height,
            output_width,int_offset
    );
#else
    vp8_sixtap_single_pass(
            cq,
            cl_data.vp8_sixtap_predict_kernel,
            cl_data.vp8_sixtap_predict_kernel_size,
            FData_height*FData_width,
            src_mem,
            dst_mem,
            src_base,
            src_offset,
            src_len,
            src_pixels_per_line,
            xoffset,
            yoffset,
            dst_base,
            dst_offset,
            dst_pitch,
            dst_len
    );
#endif


    return;
}

void vp8_sixtap_predict8x8_cl
(
    cl_command_queue cq,
    unsigned char *src_base,
    cl_mem src_mem,
    int src_offset,
    int src_pixels_per_line,
    int xoffset,
    int yoffset,
    unsigned char *dst_base,
    cl_mem dst_mem,
    int dst_offset,
    int dst_pitch
) {
    int output_width=8, output_height=8, FData_height=13, FData_width=8;

    //Size of output to transfer
    int dst_len = DST_LEN(dst_pitch,output_height,output_width);
    int src_len = SIXTAP_SRC_LEN(FData_width,FData_height,src_pixels_per_line);

#if TWO_PASS_SIXTAP
    int int_offset = 16;
    unsigned char *src_ptr = src_base + src_offset;

    vp8_sixtap_run_cl(cq, src_mem, dst_mem,
            (src_ptr-2*src_pixels_per_line),src_offset, src_len,
            src_pixels_per_line, xoffset,yoffset,dst_base,dst_offset,
            dst_pitch,dst_len,FData_height,FData_width,output_height,
            output_width,int_offset
    );
#else
    vp8_sixtap_single_pass(
            cq,
            cl_data.vp8_sixtap_predict8x8_kernel,
            cl_data.vp8_sixtap_predict8x8_kernel_size,
            FData_height*FData_width,
            src_mem,
            dst_mem,
            src_base,
            src_offset,
            src_len,
            src_pixels_per_line,
            xoffset,
            yoffset,
            dst_base,
            dst_offset,
            dst_pitch,
            dst_len
    );
#endif

    return;
}

void vp8_sixtap_predict8x4_cl
(
    cl_command_queue cq,
    unsigned char *src_base,
    cl_mem src_mem,
    int src_offset,
    int src_pixels_per_line,
    int xoffset,
    int yoffset,
    unsigned char *dst_base,
    cl_mem dst_mem,
    int dst_offset,
    int dst_pitch
) {

    int output_width=8, output_height=4, FData_height=9, FData_width=8;

    //Size of output to transfer
    int dst_len = DST_LEN(dst_pitch,output_height,output_width);
    int src_len = SIXTAP_SRC_LEN(FData_width,FData_height,src_pixels_per_line);

#if TWO_PASS_SIXTAP
    int int_offset = 16;
    unsigned char *src_ptr = src_base + src_offset;
    
    vp8_sixtap_run_cl(cq, src_mem, dst_mem,
            (src_ptr-2*src_pixels_per_line),src_offset, src_len,
            src_pixels_per_line, xoffset,yoffset,dst_base,dst_offset,
            dst_pitch,dst_len,FData_height,FData_width,output_height,
            output_width,int_offset
    );
#else
    vp8_sixtap_single_pass(
            cq,
            cl_data.vp8_sixtap_predict8x4_kernel,
            cl_data.vp8_sixtap_predict8x4_kernel_size,
            FData_height*FData_width,
            src_mem,
            dst_mem,
            src_base,
            src_offset,
            src_len,
            src_pixels_per_line,
            xoffset,
            yoffset,
            dst_base,
            dst_offset,
            dst_pitch,
            dst_len
    );
#endif

    return;
}

void vp8_sixtap_predict16x16_cl
(
    cl_command_queue cq,
    unsigned char *src_base,
    cl_mem src_mem,
    int src_offset,
    int src_pixels_per_line,
    int xoffset,
    int yoffset,
    unsigned char *dst_base,
    cl_mem dst_mem,
    int dst_offset,
    int dst_pitch
) {

    int output_width=16, output_height=16, FData_height=21, FData_width=16;

    //Size of output to transfer
    int dst_len = DST_LEN(dst_pitch,output_height,output_width);
    int src_len = SIXTAP_SRC_LEN(FData_width,FData_height,src_pixels_per_line);

#if TWO_PASS_SIXTAP
    int int_offset = 32;
    unsigned char *src_ptr = src_base + src_offset;

    vp8_sixtap_run_cl(cq, src_mem, dst_mem,
            (src_ptr-2*src_pixels_per_line),src_offset, src_len,
            src_pixels_per_line, xoffset,yoffset,dst_base,dst_offset,
            dst_pitch,dst_len,FData_height,FData_width,output_height,
            output_width,int_offset
    );
#else
    vp8_sixtap_single_pass(
            cq,
            cl_data.vp8_sixtap_predict16x16_kernel,
            cl_data.vp8_sixtap_predict16x16_kernel_size,
            FData_height*FData_width,
            src_mem,
            dst_mem,
            src_base,
            src_offset,
            src_len,
            src_pixels_per_line,
            xoffset,
            yoffset,
            dst_base,
            dst_offset,
            dst_pitch,
            dst_len
    );
#endif

    return;

}



void vp8_filter_block2d_bil_first_pass_cl(
    cl_command_queue cq,
    unsigned char *src_base,
    cl_mem src_mem,
    int src_offset,
    cl_mem int_mem,
    int src_pixels_per_line,
    int height,
    int width,
    int xoffset
)
{
    int err;
    size_t global = width*height;
    int free_src = 0;

    if (src_mem == NULL){
        int src_len = BIL_SRC_LEN(width,height,src_pixels_per_line);

        /*Make space for kernel input/output data. Initialize the buffer as well if needed. */
        VP8_CL_CREATE_BUF(cq, src_mem, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
            sizeof (unsigned char) * src_len, src_base+src_offset,,
        );
        src_offset = 0; //Set to zero as long as src_mem starts at base+offset
        free_src = 1;
    }

    err =  clSetKernelArg(cl_data.vp8_filter_block2d_bil_first_pass_kernel, 0, sizeof (cl_mem), &src_mem);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_bil_first_pass_kernel, 1, sizeof (int), &src_offset);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_bil_first_pass_kernel, 2, sizeof (cl_mem), &int_mem);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_bil_first_pass_kernel, 3, sizeof (int), &src_pixels_per_line);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_bil_first_pass_kernel, 4, sizeof (int), &height);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_bil_first_pass_kernel, 5, sizeof (int), &width);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_bil_first_pass_kernel, 6, sizeof (int), &xoffset);
    VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",
        ,
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel( cq, cl_data.vp8_filter_block2d_bil_first_pass_kernel, 1, NULL, &global, NULL , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);,
    );

    if (free_src == 1)
        clReleaseMemObject(src_mem);
}


void vp8_filter_block2d_bil_second_pass_cl(
    cl_command_queue cq,
    cl_mem int_mem,
    unsigned char *dst_base,
    cl_mem dst_mem,
    int dst_offset,
    int dst_pitch,
    int height,
    int width,
    int yoffset
)
{
    int err;
    size_t global = width*height;

    //Size of output data
    int dst_len = DST_LEN(dst_pitch,height,width);

    int free_dst = 0;
    if (dst_mem == NULL){
        VP8_CL_CREATE_BUF(cq, dst_mem, CL_MEM_WRITE_ONLY|CL_MEM_COPY_HOST_PTR,
            sizeof (unsigned char) * dst_len + dst_offset, dst_base,,
        );
        free_dst = 1;
    }

    err =  clSetKernelArg(cl_data.vp8_filter_block2d_bil_second_pass_kernel, 0, sizeof (cl_mem), &int_mem);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_bil_second_pass_kernel, 1, sizeof (cl_mem), &dst_mem);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_bil_second_pass_kernel, 2, sizeof (int), &dst_offset);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_bil_second_pass_kernel, 3, sizeof (int), &dst_pitch);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_bil_second_pass_kernel, 4, sizeof (int), &height);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_bil_second_pass_kernel, 5, sizeof (int), &width);
    err |= clSetKernelArg(cl_data.vp8_filter_block2d_bil_second_pass_kernel, 6, sizeof (int), &yoffset);
    VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
        "Error: Failed to set kernel arguments!\n",
        ,
    );

    /* Execute the kernel */
    err = clEnqueueNDRangeKernel( cq, cl_data.vp8_filter_block2d_bil_second_pass_kernel, 1, NULL, &global, NULL , 0, NULL, NULL);
    VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
        "Error: Failed to execute kernel!\n",
        printf("err = %d\n",err);,
    );

    if (free_dst == 1){
        /* Read back the result data from the device */
        err = clEnqueueReadBuffer(cq, dst_mem, CL_FALSE, 0, sizeof (unsigned char) * dst_len + dst_offset, dst_base, 0, NULL, NULL);
        VP8_CL_CHECK_SUCCESS( cq, err != CL_SUCCESS,
            "Error: Failed to read output array!\n",
            ,
        );
        clReleaseMemObject(dst_mem);
    }

}

void vp8_bilinear_predict4x4_cl
(
    cl_command_queue cq,
    unsigned char *src_base,
    cl_mem src_mem,
    int src_offset,
    int src_pixels_per_line,
    int xoffset,
    int yoffset,
    unsigned char *dst_base,
    cl_mem dst_mem,
    int dst_offset,
    int dst_pitch
) {

    const int height = 4, width = 4;

#if !STATIC_MEM
    int err;
    cl_mem int_mem = NULL;
    VP8_CL_CREATE_BUF(NULL, int_mem, NULL, sizeof(cl_int)*21*16, NULL, ,);
#endif
    
    /* First filter 1-D horizontally... */
    vp8_filter_block2d_bil_first_pass_cl(cq, src_base, src_mem, src_offset, int_mem, src_pixels_per_line, height + 1, width, xoffset);

    /* then 1-D vertically... */
    vp8_filter_block2d_bil_second_pass_cl(cq, int_mem, dst_base, dst_mem, dst_offset, dst_pitch, height, width, yoffset);

#if !STATIC_MEM
    clReleaseMemObject(int_mem);
#endif

}

void vp8_bilinear_predict8x8_cl
(
    cl_command_queue cq,
    unsigned char *src_base,
    cl_mem src_mem,
    int src_offset,
    int src_pixels_per_line,
    int xoffset,
    int yoffset,
    unsigned char *dst_base,
    cl_mem dst_mem,
    int dst_offset,
    int dst_pitch
) {

    const int height = 8, width = 8;

#if !STATIC_MEM
    int err;
    cl_mem int_mem = NULL;
    VP8_CL_CREATE_BUF(NULL, int_mem, NULL, sizeof(cl_int)*21*16, NULL, ,);
#endif

    /* First filter 1-D horizontally... */
    vp8_filter_block2d_bil_first_pass_cl(cq, src_base, src_mem, src_offset, int_mem, src_pixels_per_line, height + 1, width, xoffset);

    /* then 1-D vertically... */
    vp8_filter_block2d_bil_second_pass_cl(cq, int_mem, dst_base, dst_mem, dst_offset, dst_pitch, height, width, yoffset);

#if !STATIC_MEM
    clReleaseMemObject(int_mem);
#endif
    
}

void vp8_bilinear_predict8x4_cl
(
    cl_command_queue cq,
    unsigned char *src_base,
    cl_mem src_mem,
    int src_offset,
    int src_pixels_per_line,
    int xoffset,
    int yoffset,
    unsigned char *dst_base,
    cl_mem dst_mem,
    int dst_offset,
    int dst_pitch
) {

    const int height = 4, width = 8;

#if !STATIC_MEM
    int err;
    cl_mem int_mem = NULL;
    VP8_CL_CREATE_BUF(NULL, int_mem, NULL, sizeof(cl_int)*21*16, NULL, ,);
#endif

    /* First filter 1-D horizontally... */
    vp8_filter_block2d_bil_first_pass_cl(cq, src_base, src_mem, src_offset, int_mem, src_pixels_per_line, height + 1, width, xoffset);

    /* then 1-D vertically... */
    vp8_filter_block2d_bil_second_pass_cl(cq, int_mem, dst_base, dst_mem, dst_offset, dst_pitch, height, width, yoffset);

#if !STATIC_MEM
    clReleaseMemObject(int_mem);
#endif

}

void vp8_bilinear_predict16x16_cl
(
    cl_command_queue cq,
    unsigned char *src_base,
    cl_mem src_mem,
    int src_offset,
    int src_pixels_per_line,
    int xoffset,
    int yoffset,
    unsigned char *dst_base,
    cl_mem dst_mem,
    int dst_offset,
    int dst_pitch
) {

    const int height = 16, width = 16;

#if !STATIC_MEM
    int err;
    cl_mem int_mem = NULL;
    VP8_CL_CREATE_BUF(NULL, int_mem, NULL, sizeof(cl_int)*21*16, NULL, ,);
#endif

    /* First filter 1-D horizontally... */
    vp8_filter_block2d_bil_first_pass_cl(cq, src_base, src_mem, src_offset, int_mem, src_pixels_per_line, height + 1, width, xoffset);

    /* then 1-D vertically... */
    vp8_filter_block2d_bil_second_pass_cl(cq, int_mem, dst_base, dst_mem, dst_offset, dst_pitch, height, width, yoffset);

#if !STATIC_MEM
    clReleaseMemObject(int_mem);
#endif

}
