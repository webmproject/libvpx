/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include <stddef.h>

#include "treewriter.h"
#include "tokenize.h"
#include "onyx_int.h"

#define ct_assert(name,cond) \
    static void assert_##name(void) UNUSED;\
    static void assert_##name(void) {switch(0){case 0:case !!(cond):;}}

#define DEFINE(sym, val) int sym = val;

/*
#define BLANK() asm volatile("\n->" : : )
*/

/*
 * int main(void)
 * {
 */

//pack tokens
DEFINE(vp8_writer_lowvalue,                     offsetof(vp8_writer, lowvalue));
DEFINE(vp8_writer_range,                        offsetof(vp8_writer, range));
DEFINE(vp8_writer_value,                        offsetof(vp8_writer, value));
DEFINE(vp8_writer_count,                        offsetof(vp8_writer, count));
DEFINE(vp8_writer_pos,                          offsetof(vp8_writer, pos));
DEFINE(vp8_writer_buffer,                       offsetof(vp8_writer, buffer));

DEFINE(tokenextra_token,                        offsetof(TOKENEXTRA, Token));
DEFINE(tokenextra_extra,                        offsetof(TOKENEXTRA, Extra));
DEFINE(tokenextra_context_tree,                 offsetof(TOKENEXTRA, context_tree));
DEFINE(tokenextra_skip_eob_node,                offsetof(TOKENEXTRA, skip_eob_node));
DEFINE(TOKENEXTRA_SZ,                           sizeof(TOKENEXTRA));

DEFINE(vp8_extra_bit_struct_sz,                 sizeof(vp8_extra_bit_struct));

DEFINE(vp8_token_value,                         offsetof(vp8_token, value));
DEFINE(vp8_token_len,                           offsetof(vp8_token, Len));

DEFINE(vp8_extra_bit_struct_tree,               offsetof(vp8_extra_bit_struct, tree));
DEFINE(vp8_extra_bit_struct_prob,               offsetof(vp8_extra_bit_struct, prob));
DEFINE(vp8_extra_bit_struct_len,                offsetof(vp8_extra_bit_struct, Len));
DEFINE(vp8_extra_bit_struct_base_val,           offsetof(vp8_extra_bit_struct, base_val));

DEFINE(vp8_comp_tplist,                         offsetof(VP8_COMP, tplist));
DEFINE(vp8_comp_common,                         offsetof(VP8_COMP, common));
DEFINE(vp8_comp_bc2,                            offsetof(VP8_COMP, bc2));

DEFINE(tokenlist_start,                         offsetof(TOKENLIST, start));
DEFINE(tokenlist_stop,                          offsetof(TOKENLIST, stop));
DEFINE(TOKENLIST_SZ,                            sizeof(TOKENLIST));

DEFINE(vp8_common_mb_rows,                      offsetof(VP8_COMMON, mb_rows));

// offsets from BLOCK structure
DEFINE(vp8_block_coeff,                         offsetof(BLOCK, coeff));
DEFINE(vp8_block_quant_fast,                    offsetof(BLOCK, quant_fast));
DEFINE(vp8_block_round,                         offsetof(BLOCK, round));

// offsets from BLOCKD structure
DEFINE(vp8_blockd_qcoeff,                       offsetof(BLOCKD, qcoeff));
DEFINE(vp8_blockd_dqcoeff,                      offsetof(BLOCKD, dqcoeff));
DEFINE(vp8_blockd_dequant,                      offsetof(BLOCKD, dequant));
DEFINE(vp8_blockd_eob,                          offsetof(BLOCKD, eob));

// These two sizes are used in vp8cx_pack_tokens.  They are hard coded
// so if the size changes this will have to be adjusted.
#if HAVE_ARMV5TE
ct_assert(TOKENEXTRA_SZ, sizeof(TOKENEXTRA) == 8)
ct_assert(vp8_extra_bit_struct_sz, sizeof(vp8_extra_bit_struct) == 16)
#endif

//add asserts for any offset that is not supported by assembly code
//add asserts for any size that is not supported by assembly code
/*
 * return 0;
 * }
 */
