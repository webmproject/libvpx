/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_ENTROPYMV_H
#define __INC_ENTROPYMV_H

#include "treecoder.h"
#include "vpx_config.h"
#include "blockd.h"

enum {
  mv_max  = 1023,              /* max absolute value of a MV component */
  MVvals = (2 * mv_max) + 1,   /* # possible values "" */
  mvlong_width = 10,       /* Large MVs have 9 bit magnitudes */
  mvnum_short = 8,         /* magnitudes 0 through 7 */
  mvnum_short_bits = 3,         /* number of bits for short mvs */

  mvfp_max  = 255,              /* max absolute value of a full pixel MV component */
  MVfpvals = (2 * mvfp_max) + 1, /* # possible full pixel MV values */

  /* probability offsets for coding each MV component */

  mvpis_short = 0,         /* short (<= 7) vs long (>= 8) */
  MVPsign,                /* sign for non-zero */
  MVPshort,               /* 8 short values = 7-position tree */

  MVPbits = MVPshort + mvnum_short - 1, /* mvlong_width long value bits */
  MVPcount = MVPbits + mvlong_width    /* (with independent probabilities) */
};

typedef struct mv_context {
  vp8_prob prob[MVPcount];  /* often come in row, col pairs */
} MV_CONTEXT;

extern const MV_CONTEXT vp8_mv_update_probs[2], vp8_default_mv_context[2];

enum {
  mv_max_hp  = 2047,              /* max absolute value of a MV component */
  MVvals_hp = (2 * mv_max_hp) + 1,   /* # possible values "" */
  mvlong_width_hp = 11,       /* Large MVs have 9 bit magnitudes */
  mvnum_short_hp = 16,         /* magnitudes 0 through 15 */
  mvnum_short_bits_hp = 4,         /* number of bits for short mvs */

  mvfp_max_hp  = 255,              /* max absolute value of a full pixel MV component */
  MVfpvals_hp = (2 * mvfp_max_hp) + 1, /* # possible full pixel MV values */

  /* probability offsets for coding each MV component */

  mvpis_short_hp = 0,         /* short (<= 7) vs long (>= 8) */
  MVPsign_hp,                /* sign for non-zero */
  MVPshort_hp,               /* 8 short values = 7-position tree */

  MVPbits_hp = MVPshort_hp + mvnum_short_hp - 1, /* mvlong_width long value bits */
  MVPcount_hp = MVPbits_hp + mvlong_width_hp    /* (with independent probabilities) */
};

typedef struct mv_context_hp {
  vp8_prob prob[MVPcount_hp];  /* often come in row, col pairs */
} MV_CONTEXT_HP;

extern const MV_CONTEXT_HP vp8_mv_update_probs_hp[2], vp8_default_mv_context_hp[2];

extern const vp8_tree_index vp8_small_mvtree[];
extern struct vp8_token_struct vp8_small_mvencodings [8];
extern const vp8_tree_index vp8_small_mvtree_hp[];
extern struct vp8_token_struct vp8_small_mvencodings_hp [16];

void vp8_entropy_mv_init();
struct VP8Common;
void vp8_adapt_mv_probs(struct VP8Common *cm);

#endif
