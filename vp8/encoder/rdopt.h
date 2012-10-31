/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_RDOPT_H
#define __INC_RDOPT_H

#define RDCOST(RM,DM,R,D) ( ((128+((int64_t)R)*(RM)) >> 8) + ((int64_t)DM)*(D) )
#define RDCOST_8x8(RM,DM,R,D) ( ((128+((int64_t)R)*(RM)) >> 8) + ((int64_t)DM)*(D) )

extern void vp9_initialize_rd_consts(VP9_COMP *cpi, int Qvalue);

extern void vp9_rd_pick_inter_mode(VP9_COMP *cpi, MACROBLOCK *x,
                                   int recon_yoffset, int recon_uvoffset,
                                   int *returnrate, int *returndistortion,
                                   int64_t *returnintra);

extern void vp9_rd_pick_intra_mode(VP9_COMP *cpi, MACROBLOCK *x,
                                   int *r, int *d);

extern void vp9_rd_pick_intra_mode_sb(VP9_COMP *cpi, MACROBLOCK *x,
                                      int *r, int *d);

extern void vp9_mv_pred(VP9_COMP *cpi, MACROBLOCKD *xd,
                        const MODE_INFO *here, int_mv *mvp,
                        int refframe, int *ref_frame_sign_bias,
                        int *sr, int near_sadidx[]);

extern void vp9_init_me_luts();

extern void vp9_set_mbmode_and_mvs(MACROBLOCK *x,
                                   MB_PREDICTION_MODE mb, int_mv *mv);

#endif
