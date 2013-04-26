/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_ENTROPY_H_
#define VP9_COMMON_VP9_ENTROPY_H_

#include "vpx/vpx_integer.h"
#include "vp9/common/vp9_treecoder.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_common.h"

extern const int vp9_i8x8_block[4];

/* Coefficient token alphabet */

#define ZERO_TOKEN              0       /* 0         Extra Bits 0+0 */
#define ONE_TOKEN               1       /* 1         Extra Bits 0+1 */
#define TWO_TOKEN               2       /* 2         Extra Bits 0+1 */
#define THREE_TOKEN             3       /* 3         Extra Bits 0+1 */
#define FOUR_TOKEN              4       /* 4         Extra Bits 0+1 */
#define DCT_VAL_CATEGORY1       5       /* 5-6       Extra Bits 1+1 */
#define DCT_VAL_CATEGORY2       6       /* 7-10      Extra Bits 2+1 */
#define DCT_VAL_CATEGORY3       7       /* 11-18     Extra Bits 3+1 */
#define DCT_VAL_CATEGORY4       8       /* 19-34     Extra Bits 4+1 */
#define DCT_VAL_CATEGORY5       9       /* 35-66     Extra Bits 5+1 */
#define DCT_VAL_CATEGORY6       10      /* 67+       Extra Bits 14+1 */
#define DCT_EOB_TOKEN           11      /* EOB       Extra Bits 0+0 */
#define MAX_ENTROPY_TOKENS      12
#define ENTROPY_NODES           11
#define EOSB_TOKEN              127     /* Not signalled, encoder only */

#define INTER_MODE_CONTEXTS     7

extern const vp9_tree_index vp9_coef_tree[];

extern struct vp9_token vp9_coef_encodings[MAX_ENTROPY_TOKENS];

typedef struct {
  vp9_tree_p tree;
  const vp9_prob *prob;
  int len;
  int base_val;
} vp9_extra_bit;

extern vp9_extra_bit vp9_extra_bits[12];    /* indexed by token value */

#define PROB_UPDATE_BASELINE_COST   7

#define MAX_PROB                255
#define DCT_MAX_VALUE           16384

/* Coefficients are predicted via a 3-dimensional probability table. */

/* Outside dimension.  0 = Y with DC, 1 = UV */
#define BLOCK_TYPES 2
#define REF_TYPES 2  // intra=0, inter=1

/* Middle dimension reflects the coefficient position within the transform. */
#define COEF_BANDS 6

/* Inside dimension is measure of nearby complexity, that reflects the energy
   of nearby coefficients are nonzero.  For the first coefficient (DC, unless
   block type is 0), we look at the (already encoded) blocks above and to the
   left of the current block.  The context index is then the number (0,1,or 2)
   of these blocks having nonzero coefficients.
   After decoding a coefficient, the measure is determined by the size of the
   most recently decoded coefficient.
   Note that the intuitive meaning of this measure changes as coefficients
   are decoded, e.g., prior to the first token, a zero means that my neighbors
   are empty while, after the first token, because of the use of end-of-block,
   a zero means we just decoded a zero and hence guarantees that a non-zero
   coefficient will appear later in this block.  However, this shift
   in meaning is perfectly OK because our context depends also on the
   coefficient band (and since zigzag positions 0, 1, and 2 are in
   distinct bands). */

/*# define DC_TOKEN_CONTEXTS        3*/ /* 00, 0!0, !0!0 */
#define PREV_COEF_CONTEXTS          6

typedef unsigned int vp9_coeff_count[REF_TYPES][COEF_BANDS][PREV_COEF_CONTEXTS]
                                    [MAX_ENTROPY_TOKENS];
typedef unsigned int vp9_coeff_stats[REF_TYPES][COEF_BANDS][PREV_COEF_CONTEXTS]
                                    [ENTROPY_NODES][2];
typedef vp9_prob vp9_coeff_probs[REF_TYPES][COEF_BANDS][PREV_COEF_CONTEXTS]
                                [ENTROPY_NODES];

#define SUBEXP_PARAM                4   /* Subexponential code parameter */
#define MODULUS_PARAM               13  /* Modulus parameter */

struct VP9Common;
void vp9_default_coef_probs(struct VP9Common *);
extern DECLARE_ALIGNED(16, const int, vp9_default_zig_zag1d_4x4[16]);

extern DECLARE_ALIGNED(16, const int, vp9_col_scan_4x4[16]);
extern DECLARE_ALIGNED(16, const int, vp9_row_scan_4x4[16]);

extern DECLARE_ALIGNED(64, const int, vp9_default_zig_zag1d_8x8[64]);

extern DECLARE_ALIGNED(16, const int, vp9_col_scan_8x8[64]);
extern DECLARE_ALIGNED(16, const int, vp9_row_scan_8x8[64]);

extern DECLARE_ALIGNED(16, const int, vp9_default_zig_zag1d_16x16[256]);

extern DECLARE_ALIGNED(16, const int, vp9_col_scan_16x16[256]);
extern DECLARE_ALIGNED(16, const int, vp9_row_scan_16x16[256]);

extern DECLARE_ALIGNED(16, const int, vp9_default_zig_zag1d_32x32[1024]);

void vp9_coef_tree_initialize(void);
void vp9_adapt_coef_probs(struct VP9Common *);

static INLINE void vp9_reset_sb_tokens_context(MACROBLOCKD* const xd,
                                               BLOCK_SIZE_TYPE bsize) {
  /* Clear entropy contexts */
  const int bw = 1 << (b_width_log2(bsize) - 2);
  const int bh = 1 << (b_height_log2(bsize) - 2);
  vpx_memset(xd->above_context, 0, sizeof(ENTROPY_CONTEXT_PLANES) * bw);
  vpx_memset(xd->left_context, 0, sizeof(ENTROPY_CONTEXT_PLANES) * bh);
}

extern const int vp9_coef_bands8x8[64];
extern const int vp9_coef_bands4x4[16];

static int get_coef_band(const int *scan, TX_SIZE tx_size, int coef_index) {
  if (tx_size == TX_4X4) {
    return vp9_coef_bands4x4[scan[coef_index]];
  } else {
    const int pos = scan[coef_index];
    const int sz = 1 << (2 + tx_size);
    const int x = pos & (sz - 1), y = pos >> (2 + tx_size);
    if (x >= 8 || y >= 8)
      return 5;
    else
      return vp9_coef_bands8x8[y * 8 + x];
  }
}
extern int vp9_get_coef_context(const int *scan, const int *neighbors,
                                int nb_pad, uint8_t *token_cache, int c, int l);
const int *vp9_get_coef_neighbors_handle(const int *scan, int *pad);

#if CONFIG_MODELCOEFPROB
#define COEFPROB_BITS               8
#define COEFPROB_MODELS             (1 << COEFPROB_BITS)

// 2 => EOB and Zero nodes are unconstrained, rest are modeled
// 3 => EOB, Zero and One nodes are unconstrained, rest are modeled
#define UNCONSTRAINED_NODES         3   // Choose one of 2 or 3

// whether forward updates are model-based
#define MODEL_BASED_UPDATE          1
// if model-based how many nodes are unconstrained
#define UNCONSTRAINED_UPDATE_NODES  3
// whether backward updates are model-based
#define MODEL_BASED_ADAPT           0
#define UNCONSTRAINED_ADAPT_NODES   3

typedef vp9_prob vp9_coeff_probs_model[REF_TYPES][COEF_BANDS]
                                      [PREV_COEF_CONTEXTS][2];
extern const vp9_prob vp9_modelcoefprobs[COEFPROB_MODELS][ENTROPY_NODES - 1];
void vp9_get_model_distribution(vp9_prob model, vp9_prob *tree_probs,
                                int b, int r);
#endif  // CONFIG_MODELCOEFPROB

#if CONFIG_CODE_ZEROGROUP

#define ZPC_STATS

typedef enum {
  HORIZONTAL = 0,
  DIAGONAL,
  VERTICAL,
} OrientationType;

/* Note EOB should become part of this symbol eventually,
 * but holding off on this for now because that is a major
 * change in the rest of the codebase */

#define ZPC_ISOLATED     (MAX_ENTROPY_TOKENS + 0)    /* Isolated zero */

/* ZPC_EOORIENT: All remaining coefficients in the same orientation are 0.
 * In other words all remaining coeffs in the current subband, and all
 * children of the current subband are zero. Subbands are defined by
 * dyadic partitioning in the coeff domain */
#define ZPC_EOORIENT     (MAX_ENTROPY_TOKENS + 1)    /* End of Orientation */

/* Band limits over which the eoo bit is sent */
#define ZPC_EOO_BAND_LOWER       0
#define ZPC_EOO_BAND_UPPER       5

#define USE_ZPC_EOORIENT         1       /* 0: not used */
                                         /* 1: used */
#define ZPC_NODES                1

#define UNKNOWN_TOKEN          255       /* Not signalled, encoder only */

#define ZPC_BANDS                3       /* context bands for izr */
#define ZPC_PTOKS                3       /* context pt for zpcs */

#define coef_to_zpc_band(b)      ((b) >> 1)
#define coef_to_zpc_ptok(p)      ((p) > 2 ? 2 : (p))

typedef vp9_prob vp9_zpc_probs[REF_TYPES][ZPC_BANDS]
                              [ZPC_PTOKS][ZPC_NODES];
typedef unsigned int vp9_zpc_count[REF_TYPES][ZPC_BANDS]
                                  [ZPC_PTOKS][ZPC_NODES][2];

OrientationType vp9_get_orientation(int rc, TX_SIZE tx_size);
int vp9_use_eoo(int c, int eob, const int *scan, TX_SIZE tx_size,
                int *is_last_zero, int *is_eoo);
int vp9_is_eoo(int c, int eob, const int *scan, TX_SIZE tx_size,
               const int16_t *qcoeff_ptr, int *last_nz_pos);

#define ZPC_USEEOO_THRESH        4
#define ZPC_ZEROSSAVED_EOO       7   /* encoder only */

void vp9_adapt_zpc_probs(struct VP9Common *cm);

#endif  // CONFIG_CODE_ZEROGROUP

static INLINE const int* get_scan_4x4(TX_TYPE tx_type) {
  switch (tx_type) {
    case ADST_DCT:
      return vp9_row_scan_4x4;
    case DCT_ADST:
      return vp9_col_scan_4x4;
    default:
      return vp9_default_zig_zag1d_4x4;
  }
}

static INLINE const int* get_scan_8x8(TX_TYPE tx_type) {
  switch (tx_type) {
    case ADST_DCT:
      return vp9_row_scan_8x8;
    case DCT_ADST:
      return vp9_col_scan_8x8;
    default:
      return vp9_default_zig_zag1d_8x8;
  }
}

static INLINE const int* get_scan_16x16(TX_TYPE tx_type) {
  switch (tx_type) {
    case ADST_DCT:
      return vp9_row_scan_16x16;
    case DCT_ADST:
      return vp9_col_scan_16x16;
    default:
      return vp9_default_zig_zag1d_16x16;
  }
}

#include "vp9/common/vp9_coefupdateprobs.h"

#endif  // VP9_COMMON_VP9_ENTROPY_H_
