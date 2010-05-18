/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#ifndef __INC_RDOPT_H
#define __INC_RDOPT_H
void vp8_initialize_rd_consts(VP8_COMP *cpi, int Qvalue);
int vp8_rd_pick_intra4x4mby_modes(VP8_COMP *cpi, MACROBLOCK *mb, int *rate, int *rate_to, int *distortion);
int vp8_rd_pick_intra16x16mby_mode(VP8_COMP *cpi, MACROBLOCK *x, int *returnrate, int *rate_to, int *returndistortion);
int vp8_rd_pick_intra_mbuv_mode(VP8_COMP *cpi, MACROBLOCK *x, int *rate, int *rate_to, int *distortion);
extern int vp8_rd_pick_inter_mode(VP8_COMP *cpi, MACROBLOCK *x, int recon_yoffset, int recon_uvoffset, int *returnrate, int *returndistortion, int *returnintra);


#endif
