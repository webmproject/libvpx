/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_ENCODEMV_H
#define __INC_ENCODEMV_H

#include "onyx_int.h"

void vp9_write_nmvprobs(VP8_COMP* const, int usehp, vp8_writer* const);
void vp9_encode_nmv(vp8_writer* const w, const MV* const mv,
                    const MV* const ref, const nmv_context* const mvctx);
void vp9_encode_nmv_fp(vp8_writer* const w, const MV* const mv,
                       const MV* const ref, const nmv_context *mvctx,
                       int usehp);
void vp9_build_nmv_cost_table(int *mvjoint,
                              int *mvcost[2],
                              const nmv_context *mvctx,
                              int usehp,
                              int mvc_flag_v,
                              int mvc_flag_h);

#endif
