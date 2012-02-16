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

#if CONFIG_HIGH_PRECISION_MV
#define MV_SHIFT 0
#else
#define MV_SHIFT 1
#endif

enum
{
#if CONFIG_HIGH_PRECISION_MV
    mv_max  = 2047,              /* max absolute value of a MV component */
    MVvals = (2 * mv_max) + 1,   /* # possible values "" */
    mvlong_width = 11,       /* Large MVs have 9 bit magnitudes */
    mvnum_short = 16,         /* magnitudes 0 through 15 */
    mvnum_short_bits = 4,         /* number of bits for short mvs */
#else
    mv_max  = 1023,              /* max absolute value of a MV component */
    MVvals = (2 * mv_max) + 1,   /* # possible values "" */
    mvlong_width = 10,       /* Large MVs have 9 bit magnitudes */
    mvnum_short = 8,         /* magnitudes 0 through 7 */
    mvnum_short_bits = 3,         /* number of bits for short mvs */
#endif

    mvfp_max  = 255,              /* max absolute value of a full pixel MV component */
    MVfpvals = (2 * mvfp_max) + 1, /* # possible full pixel MV values */

    /* probability offsets for coding each MV component */

    mvpis_short = 0,         /* short (<= 7) vs long (>= 8) */
    MVPsign,                /* sign for non-zero */
    MVPshort,               /* 8 short values = 7-position tree */

    MVPbits = MVPshort + mvnum_short - 1, /* mvlong_width long value bits */
    MVPcount = MVPbits + mvlong_width    /* (with independent probabilities) */
};

typedef struct mv_context
{
    vp8_prob prob[MVPcount];  /* often come in row, col pairs */
} MV_CONTEXT;

extern const MV_CONTEXT vp8_mv_update_probs[2], vp8_default_mv_context[2];

#endif
