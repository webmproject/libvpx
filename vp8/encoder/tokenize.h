/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#ifndef tokenize_h
#define tokenize_h

#include "entropy.h"
#include "block.h"

void vp8_tokenize_initialize();

typedef struct
{
    int Token;
    int Extra;
    const vp8_prob *context_tree;
    int skip_eob_node;
    int section;
} TOKENEXTRA;

int rd_cost_mby(MACROBLOCKD *);

#ifdef ENTROPY_STATS
void init_context_counters();
void print_context_counters();

extern _int64 context_counters[BLOCK_TYPES] [COEF_BANDS] [PREV_COEF_CONTEXTS] [vp8_coef_tokens];
#endif


#endif  /* tokenize_h */
