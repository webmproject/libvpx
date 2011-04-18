/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include "vp8/decoder/onyxd_int.h"

#include "vp8/common/opencl/vp8_opencl.h"
#include "vp8_decode_cl.h"

void vp8_arch_opencl_decode_init(VP8D_COMP *pbi)
{

    if (cl_initialized == CL_SUCCESS){
        cl_decode_init();
    }

}
