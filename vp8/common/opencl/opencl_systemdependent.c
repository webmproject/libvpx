/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_ports/config.h"
#include "../subpixel.h"
#include "subpixel_cl.h"
#include "../onyxc_int.h"
#include "vp8_opencl.h"

#if HAVE_DLOPEN
#include "dynamic_cl.h"
#endif

void vp8_arch_opencl_common_init(VP8_COMMON *ctx)
{

#if HAVE_DLOPEN

#if WIN32 //Windows .dll has no lib prefix and no extension
    	cl_loaded = load_cl("OpenCL");
#else   //But *nix needs full name
    	cl_loaded = load_cl("libOpenCL.so");
#endif

        if (cl_loaded == CL_SUCCESS)
            cl_initialized = cl_common_init();
        else
            cl_initialized = VP8_CL_TRIED_BUT_FAILED;

#else //!HAVE_DLOPEN (e.g. Apple)
        cl_initialized = cl_common_init();
#endif

}
