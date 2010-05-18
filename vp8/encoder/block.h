/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#ifndef __INC_BLOCK_H
#define __INC_BLOCK_H

#include "onyx.h"
#include "blockd.h"
#include "entropymv.h"
#include "entropy.h"
#include "vpx_ports/mem.h"

// motion search site
typedef struct
{
    MV mv;
    int offset;
} search_site;

typedef struct
{
    // 16 Y blocks, 4 U blocks, 4 V blocks each with 16 entries
    short *src_diff;
    short *coeff;

    // 16 Y blocks, 4 U blocks, 4 V blocks each with 16 entries
    short(*quant)[4];
    short(*zbin)[4];
    short(*zrun_zbin_boost);
    short(*round)[4];

    // Zbin Over Quant value
    short zbin_extra;

    unsigned char **base_src;
    int src;
    int src_stride;

//  MV  enc_mv;
    int force_empty;

} BLOCK;

typedef struct
{
    DECLARE_ALIGNED(16, short, src_diff[400]);       // 16x16 Y 8x8 U 8x8 V 4x4 2nd Y
    DECLARE_ALIGNED(16, short, coeff[400]);     // 16x16 Y 8x8 U 8x8 V 4x4 2nd Y

    // 16 Y blocks, 4 U blocks, 4 V blocks, 1 DC 2nd order block each with 16 entries
    BLOCK block[25];

    YV12_BUFFER_CONFIG src;

    MACROBLOCKD e_mbd;

    search_site *ss;
    int ss_count;
    int searches_per_step;

    int errorperbit;
    int sadperbit16;
    int sadperbit4;
    int errthresh;
    int rddiv;
    int rdmult;

    int mvcosts[2][MVvals+1];
    int *mvcost[2];
    int mvsadcosts[2][MVvals+1];
    int *mvsadcost[2];
    int mbmode_cost[2][MB_MODE_COUNT];
    int intra_uv_mode_cost[2][MB_MODE_COUNT];
    unsigned int bmode_costs[10][10][10];
    unsigned int inter_bmode_costs[B_MODE_COUNT];

    // These define limits to motion vector components to prevent them from extending outside the UMV borders
    int mv_col_min;
    int mv_col_max;
    int mv_row_min;
    int mv_row_max;

    int vector_range;    // Used to monitor limiting range of recent vectors to guide search.
    int skip;

    int encode_breakout;

    unsigned char *active_ptr;
    MV_CONTEXT *mvc;

    unsigned int token_costs[BLOCK_TYPES] [COEF_BANDS] [PREV_COEF_CONTEXTS] [vp8_coef_tokens];
    int optimize;

    void (*vp8_short_fdct4x4)(short *input, short *output, int pitch);
    void (*vp8_short_fdct8x4)(short *input, short *output, int pitch);
    void (*short_fdct4x4rd)(short *input, short *output, int pitch);
    void (*short_fdct8x4rd)(short *input, short *output, int pitch);
    void (*vp8_short_fdct4x4_ptr)(short *input, short *output, int pitch);
    void (*short_walsh4x4)(short *input, short *output, int pitch);

    void (*quantize_b)(BLOCK *b, BLOCKD *d);
    void (*quantize_brd)(BLOCK *b, BLOCKD *d);



} MACROBLOCK;


#endif
