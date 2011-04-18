/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp8_opencl.h"

#include <stdio.h>

CL_FUNCTIONS cl;
void *dll = NULL;
int cl_loaded = VP8_CL_NOT_INITIALIZED;

int close_cl(){
    int ret = dlclose(dll);

    if (ret != 0)
        fprintf(stderr, "Error closing OpenCL library: %s", dlerror());

    return ret;
}

int load_cl(char *lib_name){

    //printf("Loading OpenCL library\n");
    dll = dlopen(lib_name, RTLD_NOW|RTLD_LOCAL);
    if (dll != NULL){
        //printf("Found CL library\n");
    } else {
        //printf("Didn't find CL library\n");
        return VP8_CL_TRIED_BUT_FAILED;
    }

    CL_LOAD_FN("clGetPlatformIDs", cl.getPlatformIDs);
    CL_LOAD_FN("clGetPlatformInfo", cl.getPlatformInfo);
    CL_LOAD_FN("clGetDeviceIDs", cl.getDeviceIDs);
    CL_LOAD_FN("clGetDeviceInfo", cl.getDeviceInfo);
    CL_LOAD_FN("clCreateContext", cl.createContext);
//    CL_LOAD_FN("clCreateContextFromType", cl.createContextFromType);
//    CL_LOAD_FN("clRetainContext", cl.retainContext);
    CL_LOAD_FN("clReleaseContext", cl.releaseContext);
//    CL_LOAD_FN("clGetContextInfo", cl.getContextInfo);
    CL_LOAD_FN("clCreateCommandQueue", cl.createCommandQueue);
//    CL_LOAD_FN("clRetainCommandQueue", cl.retainCommandQueue);
    CL_LOAD_FN("clReleaseCommandQueue", cl.releaseCommandQueue);
//    CL_LOAD_FN("clGetCommandQueueInfo", cl.getCommandQueue);
    CL_LOAD_FN("clCreateBuffer", cl.createBuffer);
//    CL_LOAD_FN("clCreateImage2D", cl.createImage2D);
//    CL_LOAD_FN("clCreateImage3D", cl.createImage3D);
//    CL_LOAD_FN("clRetainMemObject", cl.retainMemObject);
    CL_LOAD_FN("clReleaseMemObject", cl.releaseMemObject);
//    CL_LOAD_FN("clGetSupportedImageFormats", cl.getSupportedImageFormats);
//    CL_LOAD_FN("clGetMemObjectInfo", cl.getMemObjectInfo);
//    CL_LOAD_FN("clGetImageInfo", cl.getImageInfo);
//    CL_LOAD_FN("clCreateSampler", cl.createSampler);
//    CL_LOAD_FN("clRetainSampler", cl.retainSampler);
//    CL_LOAD_FN("clReleaseSampler", cl.releaseSampler);
//    CL_LOAD_FN("clGetSamplerInfo", cl.getSamplerInfo);
    CL_LOAD_FN("clCreateProgramWithSource", cl.createProgramWithSource);
//    CL_LOAD_FN("clCreateProgramWithBinary", cl.createProgramWithBinary);
//    CL_LOAD_FN("clRetainProgram", cl.retainProgram);
    CL_LOAD_FN("clReleaseProgram", cl.releaseProgram);
    CL_LOAD_FN("clBuildProgram", cl.buildProgram);
//    CL_LOAD_FN("clUnloadCompiler", cl.unloadCompiler);
    CL_LOAD_FN("clGetProgramInfo", cl.getProgramInfo);
    CL_LOAD_FN("clGetProgramBuildInfo", cl.getProgramBuildInfo);
    CL_LOAD_FN("clCreateKernel", cl.createKernel);
//    CL_LOAD_FN("clCreateKernelsInProgram", cl.createKernelsInProgram);
//    CL_LOAD_FN("clRetainKernel", cl.retainKernel);
    CL_LOAD_FN("clReleaseKernel", cl.releaseKernel);
    CL_LOAD_FN("clSetKernelArg", cl.setKernelArg);
//    CL_LOAD_FN("clGetKernelInfo", cl.getKernelInfo);
    CL_LOAD_FN("clGetKernelWorkGroupInfo", cl.getKernelWorkGroupInfo);
//    CL_LOAD_FN("clWaitForEvents", cl.waitForEvents);
//    CL_LOAD_FN("clGetEventInfo", cl.getEventInfo);
//    CL_LOAD_FN("clRetainEvent", cl.retainEvent);
//    CL_LOAD_FN("clReleaseEvent", cl.releaseEvent);
//    CL_LOAD_FN("clGetEventProfilingInfo", cl.getEventProfilingInfo);
    CL_LOAD_FN("clFlush", cl.flush);
    CL_LOAD_FN("clFinish", cl.finish);
    CL_LOAD_FN("clEnqueueReadBuffer", cl.enqueueReadBuffer);
    CL_LOAD_FN("clEnqueueWriteBuffer", cl.enqueueWriteBuffer);
    CL_LOAD_FN("clEnqueueCopyBuffer", cl.enqueueCopyBuffer);
//    CL_LOAD_FN("clEnqueueReadImage", cl.enqueueReadImage);
//    CL_LOAD_FN("clEnqueueWriteImage", cl.enqueueWriteImage);
//    CL_LOAD_FN("clEnqueueCopyImage", cl.enqueueCopyImage);
//    CL_LOAD_FN("clEnqueueCopyImageToBuffer", cl.enqueueCopyImageToBuffer);
//    CL_LOAD_FN("clEnqueueCopyBufferToImage", cl.enqueueCopyBufferToImage);
//    CL_LOAD_FN("clEnqueueMapBuffer", cl.enqueueMapBuffer);
//    CL_LOAD_FN("clEnqueueMapImage", cl.enqueueMapImage);
//    CL_LOAD_FN("clEnqueueUnmapMemObject", cl.enqueueUnmapMemObject);
    CL_LOAD_FN("clEnqueueNDRangeKernel", cl.enqueueNDRAngeKernel);
//    CL_LOAD_FN("clEnqueueTask", cl.enqueueTask);
//    CL_LOAD_FN("clEnqueueNativeKernel", cl.enqueueNativeKernel);
//    CL_LOAD_FN("clEnqueueMarker", cl.enqueueMarker);
//    CL_LOAD_FN("clEnqueueWaitForEvents", cl.enqueueWaitForEvents);
    CL_LOAD_FN("clEnqueueBarrier", cl.enqueueBarrier);
//    CL_LOAD_FN("clGetExtensionFunctionAddress", cl.getExtensionFunctionAddress);

    return CL_SUCCESS;
}
