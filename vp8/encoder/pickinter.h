/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_PICKINTER_H
#define __INC_PICKINTER_H
#include "vpx_ports/config.h"
#include "vp8/common/onyxc_int.h"

#define RD_ESTIMATE(RM,DM,R,D) ( ((128+(R)*(RM)) >> 8) + (DM)*(D) )
extern int vp8_pick_intra4x4mby_modes(const VP8_ENCODER_RTCD *, MACROBLOCK *mb, int *Rate, int *Distortion);
extern int vp8_pick_intra_mbuv_mode(MACROBLOCK *mb);
extern int vp8_pick_inter_mode(VP8_COMP *cpi, MACROBLOCK *x, int recon_yoffset, int recon_uvoffset, int *returnrate, int *returndistortion, int *returnintra);
#endif
