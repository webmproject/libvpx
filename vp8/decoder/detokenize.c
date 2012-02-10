/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp8/common/type_aliases.h"
#include "vp8/common/blockd.h"
#include "onyxd_int.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem.h"
#include "detokenize.h"

#include "vp8/common/seg_common.h"

#define BOOL_DATA UINT8

#define OCB_X PREV_COEF_CONTEXTS * ENTROPY_NODES
DECLARE_ALIGNED(16, static const unsigned char, coef_bands_x[16]) =
{
    0 * OCB_X, 1 * OCB_X, 2 * OCB_X, 3 * OCB_X,
    6 * OCB_X, 4 * OCB_X, 5 * OCB_X, 6 * OCB_X,
    6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X,
    6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 7 * OCB_X
};
#if CONFIG_T8X8
DECLARE_ALIGNED(64, static const unsigned char, coef_bands_x_8x8[64]) = {
  0 * OCB_X, 1 * OCB_X, 2 * OCB_X, 3 * OCB_X, 5 * OCB_X, 4 * OCB_X, 4 * OCB_X, 5 * OCB_X,
  5 * OCB_X, 3 * OCB_X, 6 * OCB_X, 3 * OCB_X, 5 * OCB_X, 4 * OCB_X, 6 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 5 * OCB_X, 5 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
};
#endif
#define EOB_CONTEXT_NODE            0
#define ZERO_CONTEXT_NODE           1
#define ONE_CONTEXT_NODE            2
#define LOW_VAL_CONTEXT_NODE        3
#define TWO_CONTEXT_NODE            4
#define THREE_CONTEXT_NODE          5
#define HIGH_LOW_CONTEXT_NODE       6
#define CAT_ONE_CONTEXT_NODE        7
#define CAT_THREEFOUR_CONTEXT_NODE  8
#define CAT_THREE_CONTEXT_NODE      9
#define CAT_FIVE_CONTEXT_NODE       10

#define CAT1_MIN_VAL    5
#define CAT2_MIN_VAL    7
#define CAT3_MIN_VAL   11
#define CAT4_MIN_VAL   19
#define CAT5_MIN_VAL   35
#define CAT6_MIN_VAL   67
#define CAT1_PROB0    159
#define CAT2_PROB0    145
#define CAT2_PROB1    165

#define CAT3_PROB0 140
#define CAT3_PROB1 148
#define CAT3_PROB2 173

#define CAT4_PROB0 135
#define CAT4_PROB1 140
#define CAT4_PROB2 155
#define CAT4_PROB3 176

#define CAT5_PROB0 130
#define CAT5_PROB1 134
#define CAT5_PROB2 141
#define CAT5_PROB3 157
#define CAT5_PROB4 180

static const unsigned char cat6_prob[14] =
{ 129, 130, 133, 140, 153, 177, 196, 230, 243, 249, 252, 254, 254, 0 };

void vp8_reset_mb_tokens_context(MACROBLOCKD *x)
{
    /* Clear entropy contexts for Y2 blocks */
    if (x->mode_info_context->mbmi.mode != B_PRED &&
        x->mode_info_context->mbmi.mode != I8X8_PRED &&
        x->mode_info_context->mbmi.mode != SPLITMV)
    {
        vpx_memset(x->above_context, 0, sizeof(ENTROPY_CONTEXT_PLANES));
        vpx_memset(x->left_context, 0, sizeof(ENTROPY_CONTEXT_PLANES));
    }
    else
    {
        vpx_memset(x->above_context, 0, sizeof(ENTROPY_CONTEXT_PLANES)-1);
        vpx_memset(x->left_context, 0, sizeof(ENTROPY_CONTEXT_PLANES)-1);
    }
}

DECLARE_ALIGNED(16, extern const unsigned char, vp8_norm[256]);
#define FILL \
    if(count < 0) \
        VP8DX_BOOL_DECODER_FILL(count, value, bufptr, bufend);

#define NORMALIZE \
    /*if(range < 0x80)*/                            \
    { \
        shift = vp8_norm[range]; \
        range <<= shift; \
        value <<= shift; \
        count -= shift; \
    }

#define DECODE_AND_APPLYSIGN(value_to_sign) \
    split = (range + 1) >> 1; \
    bigsplit = (VP8_BD_VALUE)split << (VP8_BD_VALUE_SIZE - 8); \
    FILL \
    if ( value < bigsplit ) \
    { \
        range = split; \
        v= value_to_sign; \
    } \
    else \
    { \
        range = range-split; \
        value = value-bigsplit; \
        v = -value_to_sign; \
    } \
    range +=range;                   \
    value +=value;                   \
    count--;

#define DECODE_AND_BRANCH_IF_ZERO(probability,branch) \
    { \
        split = 1 +  ((( probability*(range-1) ) )>> 8); \
        bigsplit = (VP8_BD_VALUE)split << (VP8_BD_VALUE_SIZE - 8); \
        FILL \
        if ( value < bigsplit ) \
        { \
            range = split; \
            NORMALIZE \
            goto branch; \
        } \
        value -= bigsplit; \
        range = range - split; \
        NORMALIZE \
    }

#define DECODE_AND_LOOP_IF_ZERO(probability,branch) \
    { \
        split = 1 + ((( probability*(range-1) ) ) >> 8); \
        bigsplit = (VP8_BD_VALUE)split << (VP8_BD_VALUE_SIZE - 8); \
        FILL \
        if ( value < bigsplit ) \
        { \
            range = split; \
            NORMALIZE \
            Prob = coef_probs; \
            if(c<15) {\
            ++c; \
            Prob += coef_bands_x[c]; \
            goto branch; \
            } goto BLOCK_FINISHED; /*for malformed input */\
        } \
        value -= bigsplit; \
        range = range - split; \
        NORMALIZE \
    }
#if CONFIG_T8X8
#define DECODE_AND_LOOP_IF_ZERO_8x8_2(probability,branch) \
    { \
        split = 1 + ((( probability*(range-1) ) ) >> 8); \
        bigsplit = (VP8_BD_VALUE)split << (VP8_BD_VALUE_SIZE - 8); \
        FILL \
        if ( value < bigsplit ) \
        { \
            range = split; \
            NORMALIZE \
            Prob = coef_probs; \
            if(c<3) {\
            ++c; \
            Prob += coef_bands_x[c]; \
            goto branch; \
            } goto BLOCK_FINISHED_8x8; /*for malformed input */\
        } \
        value -= bigsplit; \
        range = range - split; \
        NORMALIZE \
    }
#define DECODE_AND_LOOP_IF_ZERO_8X8(probability,branch) \
    { \
        split = 1 + ((( probability*(range-1) ) ) >> 8); \
        bigsplit = (VP8_BD_VALUE)split << (VP8_BD_VALUE_SIZE - 8); \
        FILL \
        if ( value < bigsplit ) \
        { \
            range = split; \
            NORMALIZE \
            Prob = coef_probs; \
            if(c<63) {\
            ++c; \
            Prob += coef_bands_x_8x8[c]; \
            goto branch; \
            } goto BLOCK_FINISHED_8x8; /*for malformed input */\
        } \
        value -= bigsplit; \
        range = range - split; \
        NORMALIZE \
    }
#endif
#define DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT(val) \
    DECODE_AND_APPLYSIGN(val) \
    Prob = coef_probs + (ENTROPY_NODES*2); \
    if(c < 15){\
        qcoeff_ptr [ scan[c] ] = (INT16) v; \
        ++c; \
        goto DO_WHILE; }\
    qcoeff_ptr [ 15 ] = (INT16) v; \
    goto BLOCK_FINISHED;

#if CONFIG_T8X8
#define DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8_2(val) \
    DECODE_AND_APPLYSIGN(val) \
    Prob = coef_probs + (ENTROPY_NODES*2); \
    if(c < 3){\
        qcoeff_ptr [ scan[c] ] = (INT16) v; \
        ++c; \
        goto DO_WHILE_8x8; }\
    qcoeff_ptr [ scan[3] ] = (INT16) v; \
    goto BLOCK_FINISHED_8x8;
#define DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8(val) \
    DECODE_AND_APPLYSIGN(val) \
    Prob = coef_probs + (ENTROPY_NODES*2); \
    if(c < 63){\
        qcoeff_ptr [ scan[c] ] = (INT16) v; \
        ++c; \
        goto DO_WHILE_8x8; }\
    qcoeff_ptr [ scan[63] ] = (INT16) v; \
    goto BLOCK_FINISHED_8x8;
#endif

#define DECODE_EXTRABIT_AND_ADJUST_VAL(prob, bits_count)\
    split = 1 +  (((range-1) * prob) >> 8); \
    bigsplit = (VP8_BD_VALUE)split << (VP8_BD_VALUE_SIZE - 8); \
    FILL \
    if(value >= bigsplit)\
    {\
        range = range-split;\
        value = value-bigsplit;\
        val += ((UINT16)1<<bits_count);\
    }\
    else\
    {\
        range = split;\
    }\
    NORMALIZE

#if CONFIG_T8X8
int vp8_decode_mb_tokens_8x8(VP8D_COMP *dx, MACROBLOCKD *x)
{
    ENTROPY_CONTEXT *A = (ENTROPY_CONTEXT *)x->above_context;
    ENTROPY_CONTEXT *L = (ENTROPY_CONTEXT *)x->left_context;
    const VP8_COMMON *const oc = & dx->common;

    BOOL_DECODER *bc = x->current_bc;

    char *eobs = x->eobs;

    ENTROPY_CONTEXT *a, *a1;
    ENTROPY_CONTEXT *l, *l1;
    int i;

    int eobtotal = 0;

    register int count;

    const BOOL_DATA *bufptr;
    const BOOL_DATA *bufend;
    register unsigned int range;
    VP8_BD_VALUE value;
    const int *scan;//
    register unsigned int shift;
    UINT32 split;
    VP8_BD_VALUE bigsplit;
    INT16 *qcoeff_ptr;

    const vp8_prob *coef_probs;//
    int type;
    int stop;
    INT16 val, bits_count;
    INT16 c;
    INT16 v;
    const vp8_prob *Prob;//

    int seg_eob;
    int segment_id = x->mode_info_context->mbmi.segment_id;

    type = 3;
    i = 0;
    stop = 16;

    scan = vp8_default_zig_zag1d_8x8;
    qcoeff_ptr = &x->qcoeff[0];

    if (x->mode_info_context->mbmi.mode != B_PRED && x->mode_info_context->mbmi.mode != SPLITMV)
    {
        i = 24;
        stop = 24;
        type = 1;
        qcoeff_ptr += 24*16;
        eobtotal -= 4;
        scan = vp8_default_zig_zag1d;
    }

    bufend  = bc->user_buffer_end;
    bufptr  = bc->user_buffer;
    value   = bc->value;
    count   = bc->count;
    range   = bc->range;

    coef_probs = oc->fc.coef_probs_8x8 [type] [ 0 ] [0];

BLOCK_LOOP_8x8:
    a = A + vp8_block2above[i];
    l = L + vp8_block2left[i];

    if(i < 16)
    {
       a1 = A + vp8_block2above[i+1];
       l1 = L + vp8_block2left[i+4];
    }
    else if(i<24)
    {
      a1 = A + vp8_block2above[i+1];
      l1 = L + vp8_block2left[i+2];

    }
    c = (INT16)(!type);

//    Dest = ((A)!=0) + ((B)!=0);
    if(i==24)
    {
      VP8_COMBINEENTROPYCONTEXTS(v, *a, *l);
      if ( segfeature_active( x, segment_id, SEG_LVL_EOB ) )
      {
          seg_eob = get_segdata( x, segment_id, SEG_LVL_EOB );
      }
      else
          seg_eob = 64;
    }
    else
    {
      VP8_COMBINEENTROPYCONTEXTS_8x8(v, *a, *l, *a1, *l1);
      if ( segfeature_active( x, segment_id, SEG_LVL_EOB ) )
      {
          seg_eob = get_segdata( x, segment_id, SEG_LVL_EOB );
      }
      else
          seg_eob = 64;
    }

    Prob = coef_probs;
    Prob += v * ENTROPY_NODES;

DO_WHILE_8x8:
//#if CONFIG_SEGFEATURES
    if ( c == seg_eob )
        goto BLOCK_FINISHED_8x8;

    if(i==24)
      Prob += coef_bands_x[c];
    else
      Prob += coef_bands_x_8x8[c];
    DECODE_AND_BRANCH_IF_ZERO(Prob[EOB_CONTEXT_NODE], BLOCK_FINISHED_8x8);

CHECK_0_8x8_:
    if (i==24)
    {
      DECODE_AND_LOOP_IF_ZERO_8x8_2(Prob[ZERO_CONTEXT_NODE], CHECK_0_8x8_);
    }
    else
    {
      DECODE_AND_LOOP_IF_ZERO_8X8(Prob[ZERO_CONTEXT_NODE], CHECK_0_8x8_);
    }
    DECODE_AND_BRANCH_IF_ZERO(Prob[ONE_CONTEXT_NODE], ONE_CONTEXT_NODE_0_8x8_);
    DECODE_AND_BRANCH_IF_ZERO(Prob[LOW_VAL_CONTEXT_NODE],
                                LOW_VAL_CONTEXT_NODE_0_8x8_);
    DECODE_AND_BRANCH_IF_ZERO(Prob[HIGH_LOW_CONTEXT_NODE],
                                HIGH_LOW_CONTEXT_NODE_0_8x8_);
    DECODE_AND_BRANCH_IF_ZERO(Prob[CAT_THREEFOUR_CONTEXT_NODE],
                                CAT_THREEFOUR_CONTEXT_NODE_0_8x8_);
    DECODE_AND_BRANCH_IF_ZERO(Prob[CAT_FIVE_CONTEXT_NODE],
                                CAT_FIVE_CONTEXT_NODE_0_8x8_);
    val = CAT6_MIN_VAL;
    bits_count = 12;
    do
    {
        DECODE_EXTRABIT_AND_ADJUST_VAL(cat6_prob[bits_count], bits_count);
        bits_count -- ;
    }
    while (bits_count >= 0);
    if(i==24)
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8_2(val);
    }
    else
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8(val);
    }

CAT_FIVE_CONTEXT_NODE_0_8x8_:
    val = CAT5_MIN_VAL;
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT5_PROB4, 4);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT5_PROB3, 3);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT5_PROB2, 2);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT5_PROB1, 1);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT5_PROB0, 0);
    if(i==24)
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8_2(val);
    }
    else
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8(val);
    }

CAT_THREEFOUR_CONTEXT_NODE_0_8x8_:
    DECODE_AND_BRANCH_IF_ZERO(Prob[CAT_THREE_CONTEXT_NODE],
                            CAT_THREE_CONTEXT_NODE_0_8x8_);
    val = CAT4_MIN_VAL;
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT4_PROB3, 3);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT4_PROB2, 2);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT4_PROB1, 1);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT4_PROB0, 0);
    if(i==24)
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8_2(val);
    }
    else
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8(val);
    }

CAT_THREE_CONTEXT_NODE_0_8x8_:
    val = CAT3_MIN_VAL;
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT3_PROB2, 2);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT3_PROB1, 1);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT3_PROB0, 0);
    if(i==24)
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8_2(val);
    }
    else
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8(val);
    }

HIGH_LOW_CONTEXT_NODE_0_8x8_:
    DECODE_AND_BRANCH_IF_ZERO(Prob[CAT_ONE_CONTEXT_NODE],
                            CAT_ONE_CONTEXT_NODE_0_8x8_);
    val = CAT2_MIN_VAL;
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT2_PROB1, 1);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT2_PROB0, 0);
    if(i==24)
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8_2(val);
    }
    else
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8(val);
    }

CAT_ONE_CONTEXT_NODE_0_8x8_:
    val = CAT1_MIN_VAL;
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT1_PROB0, 0);
    if(i==24)
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8_2(val);
    }
    else
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8(val);
    }

LOW_VAL_CONTEXT_NODE_0_8x8_:
    DECODE_AND_BRANCH_IF_ZERO(Prob[TWO_CONTEXT_NODE],
                                TWO_CONTEXT_NODE_0_8x8_);
    DECODE_AND_BRANCH_IF_ZERO(Prob[THREE_CONTEXT_NODE],
                                THREE_CONTEXT_NODE_0_8x8_);
    if(i==24)
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8_2(4);
    }
    else
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8(4);
    }


THREE_CONTEXT_NODE_0_8x8_:
    if(i==24)
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8_2(3);
    }
    else
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8(3);
    }


TWO_CONTEXT_NODE_0_8x8_:
    if(i==24)
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8_2(2);
    }
    else
    {
        DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT_8x8(2);
    }


ONE_CONTEXT_NODE_0_8x8_:
    DECODE_AND_APPLYSIGN(1);
    Prob = coef_probs + ENTROPY_NODES;

    if (i==24)
    {
      if (c < 3)//15
      {
        qcoeff_ptr [ scan[c] ] = (INT16) v;
        ++c;
        goto DO_WHILE_8x8;
      }
    }
    else
    {
      if (c < 63)
      {
        qcoeff_ptr [ scan[c] ] = (INT16) v;
        ++c;
        goto DO_WHILE_8x8;
      }
    }

   if(i==24)
       qcoeff_ptr [ scan[3] ] = (INT16) v;//15
   else
       qcoeff_ptr [ scan[63] ] = (INT16) v;


BLOCK_FINISHED_8x8:
    *a = *l = ((eobs[i] = c) != !type);   // any nonzero data?
    /*if (i!=24) {
      *(A + vp8_block2above[i+1]) = *(A + vp8_block2above[i+2]) = *(A + vp8_block2above[i+3]) = *a;
      *(L + vp8_block2left[i+1]) = *(L + vp8_block2left[i+2]) = *(L + vp8_block2left[i+3]) = *l;
    }*/

    if (i!=24)
    {
      if(i==0)
      {
        *(A + vp8_block2above[1]) = *(A + vp8_block2above[4]) = *(A + vp8_block2above[5]) = *a;
        *(L + vp8_block2left[1]) = *(L + vp8_block2left[4]) = *(L + vp8_block2left[5]) = *l;
      }
      else if(i==4)
      {
        *(A + vp8_block2above[2]) = *(A + vp8_block2above[3]) = *(A + vp8_block2above[6]) = *(A + vp8_block2above[7]) = *a;
        *(L + vp8_block2left[2]) = *(L + vp8_block2left[3]) = *(L + vp8_block2left[6]) = *(L + vp8_block2left[7]) = *l;
        *(A + vp8_block2above[4]) = *(A + vp8_block2above[1]);
        *(L + vp8_block2left[4]) = *(L + vp8_block2left[1]);
      }
      else if(i==8)
      {
        *(A + vp8_block2above[9]) = *(A + vp8_block2above[12]) = *(A + vp8_block2above[13]) = *a;
        *(L + vp8_block2left[9]) = *(L + vp8_block2left[12]) = *(L + vp8_block2left[13]) = *l;

      }
      else if(i==12)
      {
        *(A + vp8_block2above[10]) = *(A + vp8_block2above[11]) = *(A + vp8_block2above[14]) = *(A + vp8_block2above[15]) = *a;
        *(L + vp8_block2left[10]) = *(L + vp8_block2left[11]) = *(L + vp8_block2left[14]) = *(L + vp8_block2left[15]) = *l;
        *(A + vp8_block2above[12]) = *(A + vp8_block2above[8]);
        *(L + vp8_block2left[12]) = *(L + vp8_block2left[8]);

      }
      else
      {
        *(A + vp8_block2above[i+1]) = *(A + vp8_block2above[i+2]) = *(A + vp8_block2above[i+3]) = *a;
        *(L + vp8_block2left[i+1]) = *(L + vp8_block2left[i+2]) = *(L + vp8_block2left[i+3]) = *l;

      }
    }

    eobtotal += c;
    qcoeff_ptr += (i==24 ? 16 : 64);

    i+=4;

    if (i < stop)
        goto BLOCK_LOOP_8x8;

    if (i > 24)
    {
        type = 0;
        i = 0;
        stop = 16;
        coef_probs = oc->fc.coef_probs_8x8 [type] [ 0 ] [0];
        qcoeff_ptr -= (24*16 + 16);
        scan = vp8_default_zig_zag1d_8x8;
        goto BLOCK_LOOP_8x8;
    }

    if (i == 16)
    {
        type = 2;
        coef_probs = oc->fc.coef_probs_8x8 [type] [ 0 ] [0];
        stop = 24;
        goto BLOCK_LOOP_8x8;
    }

    FILL
    bc->user_buffer = bufptr;
    bc->value = value;
    bc->count = count;
    bc->range = range;

    return eobtotal;

}
#endif
int vp8_decode_mb_tokens(VP8D_COMP *dx, MACROBLOCKD *xd)
{
    ENTROPY_CONTEXT *A = (ENTROPY_CONTEXT *)xd->above_context;
    ENTROPY_CONTEXT *L = (ENTROPY_CONTEXT *)xd->left_context;
    const FRAME_CONTEXT * const fc = &dx->common.fc;

    BOOL_DECODER *bc = xd->current_bc;

    char *eobs = xd->eobs;

    ENTROPY_CONTEXT *a;
    ENTROPY_CONTEXT *l;
    int i;

    int eobtotal = 0;

    register int count;

    const BOOL_DATA *bufptr;
    const BOOL_DATA *bufend;
    register unsigned int range;
    VP8_BD_VALUE value;
    const int *scan;
    register unsigned int shift;
    UINT32 split;
    VP8_BD_VALUE bigsplit;
    INT16 *qcoeff_ptr;

    const vp8_prob *coef_probs;
    int type;
    int stop;
    INT16 val, bits_count;
    INT16 c;
    INT16 v;
    const vp8_prob *Prob;

    int seg_eob = 16;
    int segment_id = xd->mode_info_context->mbmi.segment_id;

    if ( segfeature_active( xd, segment_id, SEG_LVL_EOB ) )
    {
        seg_eob = get_segdata( xd, segment_id, SEG_LVL_EOB );
    }

    type = 3;
    i = 0;
    stop = 16;

    scan = vp8_default_zig_zag1d;
    qcoeff_ptr = &xd->qcoeff[0];
    if (xd->mode_info_context->mbmi.mode != B_PRED &&
        xd->mode_info_context->mbmi.mode != I8X8_PRED &&
        xd->mode_info_context->mbmi.mode != SPLITMV)
    {
        i = 24;
        stop = 24;
        type = 1;
        qcoeff_ptr += 24*16;
        eobtotal -= 16;
    }

    bufend  = bc->user_buffer_end;
    bufptr  = bc->user_buffer;
    value   = bc->value;
    count   = bc->count;
    range   = bc->range;


    coef_probs = fc->coef_probs [type] [ 0 ] [0];

BLOCK_LOOP:
    a = A + vp8_block2above[i];
    l = L + vp8_block2left[i];

    c = (INT16)(!type);

    /*Dest = ((A)!=0) + ((B)!=0);*/
    VP8_COMBINEENTROPYCONTEXTS(v, *a, *l);
    Prob = coef_probs;
    Prob += v * ENTROPY_NODES;

DO_WHILE:
    if ( c == seg_eob )
        goto BLOCK_FINISHED;

    Prob += coef_bands_x[c];
    DECODE_AND_BRANCH_IF_ZERO(Prob[EOB_CONTEXT_NODE], BLOCK_FINISHED);

CHECK_0_:
    DECODE_AND_LOOP_IF_ZERO(Prob[ZERO_CONTEXT_NODE], CHECK_0_);
    DECODE_AND_BRANCH_IF_ZERO(Prob[ONE_CONTEXT_NODE], ONE_CONTEXT_NODE_0_);
    DECODE_AND_BRANCH_IF_ZERO(Prob[LOW_VAL_CONTEXT_NODE],
                              LOW_VAL_CONTEXT_NODE_0_);
    DECODE_AND_BRANCH_IF_ZERO(Prob[HIGH_LOW_CONTEXT_NODE],
                              HIGH_LOW_CONTEXT_NODE_0_);
    DECODE_AND_BRANCH_IF_ZERO(Prob[CAT_THREEFOUR_CONTEXT_NODE],
                              CAT_THREEFOUR_CONTEXT_NODE_0_);
    DECODE_AND_BRANCH_IF_ZERO(Prob[CAT_FIVE_CONTEXT_NODE],
                              CAT_FIVE_CONTEXT_NODE_0_);

    val = CAT6_MIN_VAL;
    bits_count = 12;

    do
    {
        DECODE_EXTRABIT_AND_ADJUST_VAL(cat6_prob[bits_count], bits_count);
        bits_count -- ;
    }
    while (bits_count >= 0);

    DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT(val);

CAT_FIVE_CONTEXT_NODE_0_:
    val = CAT5_MIN_VAL;
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT5_PROB4, 4);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT5_PROB3, 3);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT5_PROB2, 2);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT5_PROB1, 1);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT5_PROB0, 0);
    DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT(val);

CAT_THREEFOUR_CONTEXT_NODE_0_:
    DECODE_AND_BRANCH_IF_ZERO(Prob[CAT_THREE_CONTEXT_NODE],
                              CAT_THREE_CONTEXT_NODE_0_);
    val = CAT4_MIN_VAL;
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT4_PROB3, 3);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT4_PROB2, 2);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT4_PROB1, 1);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT4_PROB0, 0);
    DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT(val);

CAT_THREE_CONTEXT_NODE_0_:
    val = CAT3_MIN_VAL;
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT3_PROB2, 2);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT3_PROB1, 1);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT3_PROB0, 0);
    DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT(val);

HIGH_LOW_CONTEXT_NODE_0_:
    DECODE_AND_BRANCH_IF_ZERO(Prob[CAT_ONE_CONTEXT_NODE],
                              CAT_ONE_CONTEXT_NODE_0_);

    val = CAT2_MIN_VAL;
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT2_PROB1, 1);
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT2_PROB0, 0);
    DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT(val);

CAT_ONE_CONTEXT_NODE_0_:
    val = CAT1_MIN_VAL;
    DECODE_EXTRABIT_AND_ADJUST_VAL(CAT1_PROB0, 0);
    DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT(val);

LOW_VAL_CONTEXT_NODE_0_:
    DECODE_AND_BRANCH_IF_ZERO(Prob[TWO_CONTEXT_NODE], TWO_CONTEXT_NODE_0_);
    DECODE_AND_BRANCH_IF_ZERO(Prob[THREE_CONTEXT_NODE], THREE_CONTEXT_NODE_0_);
    DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT(4);

THREE_CONTEXT_NODE_0_:
    DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT(3);

TWO_CONTEXT_NODE_0_:
    DECODE_SIGN_WRITE_COEFF_AND_CHECK_EXIT(2);

ONE_CONTEXT_NODE_0_:
    DECODE_AND_APPLYSIGN(1);
    Prob = coef_probs + ENTROPY_NODES;

    if (c < 15)
    {
        qcoeff_ptr [ scan[c] ] = (INT16) v;
        ++c;
        goto DO_WHILE;
    }

    qcoeff_ptr [ 15 ] = (INT16) v;
BLOCK_FINISHED:
    *a = *l = ((eobs[i] = c) != !type);   /* any nonzero data? */
    eobtotal += c;
    qcoeff_ptr += 16;

    i++;

    if (i < stop)
        goto BLOCK_LOOP;

    if (i == 25)
    {
        type = 0;
        i = 0;
        stop = 16;
        coef_probs = fc->coef_probs [type] [ 0 ] [0];
        qcoeff_ptr -= (24*16 + 16);
        goto BLOCK_LOOP;
    }

    if (i == 16)
    {
        type = 2;
        coef_probs = fc->coef_probs [type] [ 0 ] [0];
        stop = 24;
        goto BLOCK_LOOP;
    }

    FILL
    bc->user_buffer = bufptr;
    bc->value = value;
    bc->count = count;
    bc->range = range;

    return eobtotal;

}
