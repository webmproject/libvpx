/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP8_DECODFRAME_CL_H
#define VP8_DECODFRAME_CL_H

#ifdef  __cplusplus
extern "C" {
#endif

#include "../onyxd_int.h"
#include "vp8/common/blockd.h"

//Implemented in decodframe_cl.c
extern void mb_init_dequantizer_cl(MACROBLOCKD *xd);
extern void vp8_decode_frame_cl_finish(VP8D_COMP *pbi);
extern void vp8_decode_macroblock_cl(VP8D_COMP *pbi, MACROBLOCKD *xd, int eobtotal);


#ifdef  __cplusplus
}
#endif

#endif  /* VP8_DECODFRAME_CL_H */
