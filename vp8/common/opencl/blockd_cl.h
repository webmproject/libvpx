/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef BLOCKD_OPENCL_H
#define BLOCKD_OPENCL_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "vp8_opencl.h"
#include "../blockd.h"

#define DIFF 0x0001
#define PREDICTOR 0x0002
#define QCOEFF 0x0004
#define DQCOEFF 0x0008
#define EOBS 0x0010
#define DEQUANT 0x0020
#define PRE_BUF 0x0040
#define DST_BUF 0x0080
    
#define BLOCK_COPY_ALL 0xffff

/*
#define BLOCK_MEM_SIZE 6
enum {
    DIFF_MEM = 0,
    PRED_MEM = 1,
    QCOEFF_MEM = 2,
    DQCOEFF_MEM = 3,
    EOBS_MEM = 4,
    DEQUANT_MEM = 5
} BLOCK_MEM_TYPES;


struct cl_block_mem{
    cl_mem gpu_mem;
    size_t size;
    void *host_mem;
};

typedef struct cl_block_mem block_mem;
*/
    
extern int vp8_cl_block_finish(BLOCKD *b, int flags);
extern int vp8_cl_block_prep(BLOCKD *b, int flags);

extern int vp8_cl_mb_prep(MACROBLOCKD *x, int flags);
extern int vp8_cl_mb_finish(MACROBLOCKD *x, int flags);

#ifdef	__cplusplus
}
#endif

#endif