/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef DYNAMIC_CL_H
#define	DYNAMIC_CL_H

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif
    
#include <dlfcn.h>

int load_cl(char *lib_name);
int close_cl();

extern int cl_loaded;

typedef cl_int(*fn_clGetPlatformIDs_t)(cl_uint, cl_platform_id *, cl_uint *);
typedef cl_int(*fn_clGetPlatformInfo_t)(cl_platform_id, cl_platform_info, size_t, void *, size_t *);
typedef cl_int(*fn_clGetDeviceIDs_t)(cl_platform_id, cl_device_type, cl_uint, cl_device_id *, cl_uint *);
typedef cl_int(*fn_clGetDeviceInfo_t)(cl_device_id, cl_device_info, size_t, void *, size_t *);
typedef cl_context(*fn_clCreateContext_t)(const cl_context_properties *, cl_uint, const cl_device_id *, void (*pfn_notify)(const char *, const void *, size_t, void *), void *, cl_int *);
typedef cl_context(*fn_clCreateContextFromType_t)(const cl_context_properties *, cl_device_type, void (*pfn_notify)(const char *, const void *, size_t, void *), void *, cl_int *);
typedef cl_int(*fn_clRetainContext_t)(cl_context);
typedef cl_int(*fn_clReleaseContext_t)(cl_context);
typedef cl_int(*fn_clGetContextInfo_t)(cl_context, cl_context_info, size_t, void *, size_t *);
typedef cl_command_queue(*fn_clCreateCommandQueue_t)(cl_context, cl_device_id, cl_command_queue_properties, cl_int *);
typedef cl_int(*fn_clRetainCommandQueue_t)(cl_command_queue);
typedef cl_int(*fn_clReleaseCommandQueue_t)(cl_command_queue);
typedef cl_int(*fn_clGetCommandQueueInfo_t)(cl_command_queue, cl_command_queue_info, size_t, void *, size_t *);
typedef cl_mem(*fn_clCreateBuffer_t)(cl_context, cl_mem_flags, size_t, void *, cl_int *);
typedef cl_mem(*fn_clCreateImage2D_t)(cl_context, cl_mem_flags, const cl_image_format *, size_t, size_t, size_t, void *, cl_int *);
typedef cl_mem(*fn_clCreateImage3D_t)(cl_context, cl_mem_flags, const cl_image_format *, size_t, size_t, size_t, size_t, size_t, void *, cl_int *);
typedef cl_int(*fn_clRetainMemObject_t)(cl_mem);
typedef cl_int(*fn_clReleaseMemObject_t)(cl_mem);
typedef cl_int(*fn_clGetSupportedImageFormats_t)(cl_context, cl_mem_flags, cl_mem_object_type, cl_uint, cl_image_format *, cl_uint *);
typedef cl_int(*fn_clGetMemObjectInfo_t)(cl_mem, cl_mem_info, size_t, void *, size_t *);
typedef cl_int(*fn_clGetImageInfo_t)(cl_mem, cl_image_info, size_t, void *, size_t *);
typedef cl_sampler(*fn_clCreateSampler_t)(cl_context, cl_bool, cl_addressing_mode, cl_filter_mode, cl_int *);
typedef cl_int(*fn_clRetainSampler_t)(cl_sampler);
typedef cl_int(*fn_clReleaseSampler_t)(cl_sampler);
typedef cl_int(*fn_clGetSamplerInfo_t)(cl_sampler, cl_sampler_info, size_t, void *, size_t *);
typedef cl_program(*fn_clCreateProgramWithSource_t)(cl_context, cl_uint, const char **, const size_t *, cl_int *);
typedef cl_program(*fn_clCreateProgramWithBinary_t)(cl_context, cl_uint, const cl_device_id *, const size_t *, const unsigned char **, cl_int *, cl_int *);
typedef cl_int(*fn_clRetainProgram_t)(cl_program);
typedef cl_int(*fn_clReleaseProgram_t)(cl_program);
typedef cl_int(*fn_clBuildProgram_t)(cl_program, cl_uint, const cl_device_id *, const char *,  void (*pfn_notify)(cl_program,void*), void *);
typedef cl_int(*fn_clUnloadCompiler_t)(void);
typedef cl_int(*fn_clGetProgramInfo_t)(cl_program, cl_program_info, size_t, void *, size_t *);
typedef cl_int(*fn_clGetProgramBuildInfo_t)(cl_program, cl_device_id, cl_program_build_info, size_t, void *, size_t *);
typedef cl_kernel(*fn_clCreateKernel_t)(cl_program, const char *, cl_int *);
typedef cl_int(*fn_clCreateKernelsInProgram_t)(cl_program, cl_uint, cl_kernel *, cl_uint *);
typedef cl_int(*fn_clRetainKernel_t)(cl_kernel);
typedef cl_int(*fn_clReleaseKernel_t)(cl_kernel);
typedef cl_int(*fn_clSetKernelArg_t)(cl_kernel, cl_uint, size_t, const void *);
typedef cl_int(*fn_clGetKernelInfo_t)(cl_kernel, cl_kernel_info, size_t, void *, size_t *);
typedef cl_int(*fn_clGetKernelWorkGroupInfo_t)(cl_kernel, cl_device_id, cl_kernel_work_group_info, size_t, void *, size_t *);
typedef cl_int(*fn_clWaitForEvents_t)(cl_uint, const cl_event *);
typedef cl_int(*fn_clGetEventInfo_t)(cl_event, cl_event_info, size_t, void *, size_t *);
typedef cl_int(*fn_clRetainEvent_t)(cl_event);
typedef cl_int(*fn_clReleaseEvent_t)(cl_event);
typedef cl_int(*fn_clGetEventProfilingInfo_t)(cl_event, cl_profiling_info, size_t, void *, size_t *);
typedef cl_int(*fn_clFlush_t)(cl_command_queue);
typedef cl_int(*fn_clFinish_t)(cl_command_queue);
typedef cl_int(*fn_clEnqueueReadBuffer_t)(cl_command_queue, cl_mem, cl_bool, size_t, size_t, void *, cl_uint, const cl_event *, cl_event *);
typedef cl_int(*fn_clEnqueueWriteBuffer_t)(cl_command_queue,  cl_mem,  cl_bool,  size_t,  size_t,  const void *,  cl_uint,  const cl_event *,  cl_event *);
typedef cl_int(*fn_clEnqueueCopyBuffer_t)(cl_command_queue,  cl_mem, cl_mem, size_t, size_t, size_t, cl_uint, const cl_event *, cl_event *);
typedef cl_int(*fn_clEnqueueReadImage_t)(cl_command_queue, cl_mem, cl_bool, const size_t *, const size_t *, size_t, size_t, void *, cl_uint, const cl_event *, cl_event *);
typedef cl_int(*fn_clEnqueueWriteImage_t)(cl_command_queue, cl_mem, cl_bool, const size_t *, const size_t *, size_t, size_t, const void *, cl_uint, const cl_event *, cl_event *);
typedef cl_int(*fn_clEnqueueCopyImage_t)(cl_command_queue, cl_mem, cl_mem, const size_t *, const size_t *, const size_t *, cl_uint, const cl_event *, cl_event *);
typedef cl_int(*fn_clEnqueueCopyImageToBuffer_t)(cl_command_queue, cl_mem, cl_mem, const size_t *, const size_t *, size_t, cl_uint, const cl_event *, cl_event *);
typedef cl_int(*fn_clEnqueueCopyBufferToImage_t)(cl_command_queue, cl_mem, cl_mem, size_t, const size_t *, const size_t *, cl_uint, const cl_event *, cl_event *);
typedef void*(*fn_clEnqueueMapBuffer_t)(cl_command_queue, cl_mem, cl_bool, cl_map_flags, size_t, size_t, cl_uint, const cl_event *, cl_event *, cl_int *);
typedef void*(*fn_clEnqueueMapImage_t)(cl_command_queue, cl_mem, cl_bool, cl_map_flags, const size_t *, const size_t *, size_t *, size_t *, cl_uint, const cl_event *, cl_event *, cl_int *);
typedef cl_int(*fn_clEnqueueUnmapMemObject_t)(cl_command_queue, cl_mem, void *, cl_uint, const cl_event *, cl_event *);
typedef cl_int(*fn_clEnqueueNDRangeKernel_t)(cl_command_queue, cl_kernel, cl_uint, const size_t *, const size_t *, const size_t *, cl_uint, const cl_event *, cl_event *);
typedef cl_int(*fn_clEnqueueTask_t)(cl_command_queue, cl_kernel, cl_uint, const cl_event *, cl_event *);
typedef cl_int(*fn_clEnqueueNativeKernel_t)(cl_command_queue,					 void (*user_func)(void *), void *, size_t, cl_uint, const cl_mem *, const void **, cl_uint, const cl_event *, cl_event *);
typedef cl_int(*fn_clEnqueueMarker_t)(cl_command_queue, cl_event *);
typedef cl_int(*fn_clEnqueueWaitForEvents_t)(cl_command_queue, cl_uint, const cl_event *);
typedef cl_int(*fn_clEnqueueBarrier_t)(cl_command_queue);
typedef void*(*fn_clGetExtensionFunctionAddress_t)(const char *);

typedef struct CL_FUNCTIONS {
    fn_clGetPlatformIDs_t getPlatformIDs;
    fn_clGetPlatformInfo_t getPlatformInfo;
    fn_clGetDeviceIDs_t getDeviceIDs;
    fn_clGetDeviceInfo_t getDeviceInfo;
    fn_clCreateContext_t createContext;
    fn_clCreateContextFromType_t createContextFromType;
    fn_clRetainContext_t retainContext;
    fn_clReleaseContext_t releaseContext;
    fn_clGetContextInfo_t getContextInfo;
    fn_clCreateCommandQueue_t createCommandQueue;
    fn_clRetainCommandQueue_t retainCommandQueue;
    fn_clReleaseCommandQueue_t releaseCommandQueue;
    fn_clGetCommandQueueInfo_t getCommandQueue;
    fn_clCreateBuffer_t createBuffer;
    fn_clCreateImage2D_t createImage2D;
    fn_clCreateImage3D_t createImage3D;
    fn_clRetainMemObject_t retainMemObject;
    fn_clReleaseMemObject_t releaseMemObject;
    fn_clGetSupportedImageFormats_t getSupportedImageFormats;
    fn_clGetMemObjectInfo_t getMemObjectInfo;
    fn_clGetImageInfo_t getImageInfo;
    fn_clCreateSampler_t createSampler;
    fn_clRetainSampler_t retainSampler;
    fn_clReleaseSampler_t releaseSampler;
    fn_clGetSamplerInfo_t getSamplerInfo;
    fn_clCreateProgramWithSource_t createProgramWithSource;
    fn_clCreateProgramWithBinary_t createProgramWithBinary;
    fn_clRetainProgram_t retainProgram;
    fn_clReleaseProgram_t releaseProgram;
    fn_clBuildProgram_t buildProgram;
    fn_clUnloadCompiler_t unloadCompiler;
    fn_clGetProgramInfo_t getProgramInfo;
    fn_clGetProgramBuildInfo_t getProgramBuildInfo;
    fn_clCreateKernel_t createKernel;
    fn_clCreateKernelsInProgram_t createKernelsInProgram;
    fn_clRetainKernel_t retainKernel;
    fn_clReleaseKernel_t releaseKernel;
    fn_clSetKernelArg_t setKernelArg;
    fn_clGetKernelInfo_t getKernelInfo;
    fn_clGetKernelWorkGroupInfo_t getKernelWorkGroupInfo;
    fn_clWaitForEvents_t waitForEvents;
    fn_clGetEventInfo_t getEventInfo;
    fn_clRetainEvent_t retainEvent;
    fn_clReleaseEvent_t releaseEvent;
    fn_clGetEventProfilingInfo_t getEventProfilingInfo;
    fn_clFlush_t flush;
    fn_clFinish_t finish;
    fn_clEnqueueReadBuffer_t enqueueReadBuffer;
    fn_clEnqueueWriteBuffer_t enqueueWriteBuffer;
    fn_clEnqueueCopyBuffer_t enqueueCopyBuffer;
    fn_clEnqueueReadImage_t enqueueReadImage;
    fn_clEnqueueWriteImage_t enqueueWriteImage;
    fn_clEnqueueCopyImage_t enqueueCopyImage;
    fn_clEnqueueCopyImageToBuffer_t enqueueCopyImageToBuffer;
    fn_clEnqueueCopyBufferToImage_t enqueueCopyBufferToImage;
    fn_clEnqueueMapBuffer_t enqueueMapBuffer;
    fn_clEnqueueMapImage_t enqueueMapImage;
    fn_clEnqueueUnmapMemObject_t enqueueUnmapMemObject;
    fn_clEnqueueNDRangeKernel_t enqueueNDRAngeKernel;
    fn_clEnqueueTask_t enqueueTask;
    fn_clEnqueueNativeKernel_t enqueueNativeKernel;
    fn_clEnqueueMarker_t enqueueMarker;
    fn_clEnqueueWaitForEvents_t enqueueWaitForEvents;
    fn_clEnqueueBarrier_t enqueueBarrier;
    fn_clGetExtensionFunctionAddress_t getExtensionFunctionAddress;
} CL_FUNCTIONS;

extern CL_FUNCTIONS cl;

#define clGetPlatformIDs cl.getPlatformIDs
#define clGetPlatformInfo cl.getPlatformInfo
#define clGetDeviceIDs cl.getDeviceIDs
#define clGetDeviceInfo cl.getDeviceInfo
#define clCreateContext cl.createContext
#define clCreateContextFromType cl.createContextFromType
#define clRetainContext cl.retainContext
#define clReleaseContext cl.releaseContext
#define clGetContextInfo cl.getContextInfo
#define clCreateCommandQueue cl.createCommandQueue
#define clRetainCommandQueue cl.retainCommandQueue
#define clReleaseCommandQueue cl.releaseCommandQueue
#define clGetCommandQueueInfo cl.getCommandQueue
#define clCreateBuffer cl.createBuffer
#define clCreateSubBuffer cl.createSubBuffer
#define clCreateImage2D cl.createImage2D
#define clCreateImage3D cl.createImage3D
#define clRetainMemObject cl.retainMemObject
#define clReleaseMemObject cl.releaseMemObject
#define clGetSupportedImageFormats cl.getSupportedImageFormats
#define clGetMemObjectInfo cl.getMemObjectInfo
#define clGetImageInfo cl.getImageInfo
#define clSetMemObjectDestructorCallback cl.setMemObjectDestructorCallback
#define clCreateSampler cl.createSampler
#define clRetainSampler cl.retainSampler
#define clReleaseSampler cl.releaseSampler
#define clGetSamplerInfo cl.getSamplerInfo
#define clCreateProgramWithSource cl.createProgramWithSource
#define clCreateProgramWithBinary cl.createProgramWithBinary
#define clRetainProgram cl.retainProgram
#define clReleaseProgram cl.releaseProgram
#define clBuildProgram cl.buildProgram
#define clUnloadCompiler cl.unloadCompiler
#define clGetProgramInfo cl.getProgramInfo
#define clGetProgramBuildInfo cl.getProgramBuildInfo
#define clCreateKernel cl.createKernel
#define clCreateKernelsInProgram cl.createKernelsInProgram
#define clRetainKernel cl.retainKernel
#define clReleaseKernel cl.releaseKernel
#define clSetKernelArg cl.setKernelArg
#define clGetKernelInfo cl.getKernelInfo
#define clGetKernelWorkGroupInfo cl.getKernelWorkGroupInfo
#define clWaitForEvents cl.waitForEvents
#define clGetEventInfo cl.getEventInfo
#define clCreateUserEvent cl.createUserEvent
#define clRetainEvent cl.retainEvent
#define clReleaseEvent cl.releaseEvent
#define clSetUserEventStatus cl.setUserEventStatus
#define clSetEventCallback cl.setEventCallback
#define clGetEventProfilingInfo cl.getEventProfilingInfo
#define clFlush cl.flush
#define clFinish cl.finish
#define clEnqueueReadBuffer cl.enqueueReadBuffer
#define clEnqueueReadBufferRect cl.enqueueReadBufferRect
#define clEnqueueWriteBuffer cl.enqueueWriteBuffer
#define clEnqueueWriteBufferRect cl.enqueueWriteBufferRect
#define clEnqueueCopyBuffer cl.enqueueCopyBuffer
#define clEnqueueCopyBufferRect cl.enqueueCopyBufferRect
#define clEnqueueReadImage cl.enqueueReadImage
#define clEnqueueWriteImage cl.enqueueWriteImage
#define clEnqueueCopyImage cl.enqueueCopyImage
#define clEnqueueCopyImageToBuffer cl.enqueueCopyImageToBuffer
#define clEnqueueCopyBufferToImage cl.enqueueCopyBufferToImage
#define clEnqueueMapBuffer cl.enqueueMapBuffer
#define clEnqueueMapImage cl.enqueueMapImage
#define clEnqueueUnmapMemObject cl.enqueueUnmapMemObject
#define clEnqueueNDRangeKernel cl.enqueueNDRAngeKernel
#define clEnqueueTask cl.enqueueTask
#define clEnqueueNativeKernel cl.enqueueNativeKernel
#define clEnqueueMarker cl.enqueueMarker
#define clEnqueueWaitForEvents cl.enqueueWaitForEvents
#define clEnqueueBarrier cl.enqueueBarrier
#define clGetExtensionFunctionAddress cl.getExtensionFunctionAddress

#define CL_LOAD_FN(name, ref) \
    ref = dlsym(dll,name); \
    if (ref == NULL){ \
        dlclose(dll); \
        return CL_INVALID_PLATFORM; \
    }


#ifdef	__cplusplus
}
#endif

#endif	/* DYNAMIC_CL_H */
