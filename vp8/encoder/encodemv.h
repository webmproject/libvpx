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

#if CONFIG_NEWMVENTROPY
void vp8_write_nmvprobs(VP8_COMP* const, int usehp, vp8_writer* const);
void vp8_encode_nmv(vp8_writer* const w, const MV* const mv,
                    const MV* const ref, const nmv_context* const mvctx);
void vp8_encode_nmv_fp(vp8_writer* const w, const MV* const mv,
                       const MV* const ref, const nmv_context *mvctx,
                       int usehp);
void vp8_build_nmv_cost_table(int *mvjoint,
                              int *mvcost[2],
                              const nmv_context *mvctx,
                              int usehp,
                              int mvc_flag_v,
                              int mvc_flag_h);
#else  /* CONFIG_NEWMVENTROPY */
void vp8_write_mvprobs(VP8_COMP* const, vp8_writer* const);
void vp8_encode_motion_vector(vp8_writer* const, const MV* const,
                              const MV_CONTEXT* const);
void vp8_build_component_cost_table(int *mvcost[2],
                                    const MV_CONTEXT*,
                                    const int mvc_flag[2]);
void vp8_write_mvprobs_hp(VP8_COMP* const, vp8_writer* const);
void vp8_encode_motion_vector_hp(vp8_writer* const, const MV* const,
                                 const MV_CONTEXT_HP* const);
void vp8_build_component_cost_table_hp(int *mvcost[2],
                                       const MV_CONTEXT_HP*,
                                       const int mvc_flag[2]);
#endif  /* CONFIG_NEWMVENTROPY */

#endif
