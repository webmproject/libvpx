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

#define DCT_EOB_MODEL_TOKEN     3      /* EOB       Extra Bits 0+0 */
extern const vp9_tree_index vp9_coefmodel_tree[];

extern struct vp9_token vp9_coef_encodings[MAX_ENTROPY_TOKENS];

typedef struct {
  vp9_tree_p tree;
  const vp9_prob *prob;
  int len;
  int base_val;
} vp9_extra_bit;

extern const vp9_extra_bit vp9_extra_bits[12];    /* indexed by token value */

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

#define PREV_COEF_CONTEXTS          6

// #define ENTROPY_STATS

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
extern DECLARE_ALIGNED(16, const int16_t, vp9_default_scan_4x4[16]);

extern DECLARE_ALIGNED(16, const int16_t, vp9_col_scan_4x4[16]);
extern DECLARE_ALIGNED(16, const int16_t, vp9_row_scan_4x4[16]);

extern DECLARE_ALIGNED(16, const int16_t, vp9_default_scan_8x8[64]);

extern DECLARE_ALIGNED(16, const int16_t, vp9_col_scan_8x8[64]);
extern DECLARE_ALIGNED(16, const int16_t, vp9_row_scan_8x8[64]);

extern DECLARE_ALIGNED(16, const int16_t, vp9_default_scan_16x16[256]);

extern DECLARE_ALIGNED(16, const int16_t, vp9_col_scan_16x16[256]);
extern DECLARE_ALIGNED(16, const int16_t, vp9_row_scan_16x16[256]);

extern DECLARE_ALIGNED(16, const int16_t, vp9_default_scan_32x32[1024]);

extern DECLARE_ALIGNED(16, int16_t, vp9_default_iscan_4x4[16]);

extern DECLARE_ALIGNED(16, int16_t, vp9_col_iscan_4x4[16]);
extern DECLARE_ALIGNED(16, int16_t, vp9_row_iscan_4x4[16]);

extern DECLARE_ALIGNED(16, int16_t, vp9_default_iscan_8x8[64]);

extern DECLARE_ALIGNED(16, int16_t, vp9_col_iscan_8x8[64]);
extern DECLARE_ALIGNED(16, int16_t, vp9_row_iscan_8x8[64]);

extern DECLARE_ALIGNED(16, int16_t, vp9_default_iscan_16x16[256]);

extern DECLARE_ALIGNED(16, int16_t, vp9_col_iscan_16x16[256]);
extern DECLARE_ALIGNED(16, int16_t, vp9_row_iscan_16x16[256]);

extern DECLARE_ALIGNED(16, int16_t, vp9_default_iscan_32x32[1024]);

#define MAX_NEIGHBORS 2

extern DECLARE_ALIGNED(16, int16_t,
                       vp9_default_scan_4x4_neighbors[17 * MAX_NEIGHBORS]);
extern DECLARE_ALIGNED(16, int16_t,
                       vp9_col_scan_4x4_neighbors[17 * MAX_NEIGHBORS]);
extern DECLARE_ALIGNED(16, int16_t,
                       vp9_row_scan_4x4_neighbors[17 * MAX_NEIGHBORS]);
extern DECLARE_ALIGNED(16, int16_t,
                       vp9_col_scan_8x8_neighbors[65 * MAX_NEIGHBORS]);
extern DECLARE_ALIGNED(16, int16_t,
                       vp9_row_scan_8x8_neighbors[65 * MAX_NEIGHBORS]);
extern DECLARE_ALIGNED(16, int16_t,
                       vp9_default_scan_8x8_neighbors[65 * MAX_NEIGHBORS]);
extern DECLARE_ALIGNED(16, int16_t,
                       vp9_col_scan_16x16_neighbors[257 * MAX_NEIGHBORS]);
extern DECLARE_ALIGNED(16, int16_t,
                       vp9_row_scan_16x16_neighbors[257 * MAX_NEIGHBORS]);
extern DECLARE_ALIGNED(16, int16_t,
                       vp9_default_scan_16x16_neighbors[257 * MAX_NEIGHBORS]);
extern DECLARE_ALIGNED(16, int16_t,
                       vp9_default_scan_32x32_neighbors[1025 * MAX_NEIGHBORS]);

void vp9_coef_tree_initialize(void);
void vp9_adapt_coef_probs(struct VP9Common *);

static INLINE void reset_skip_context(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  const int bw = 1 << b_width_log2(bsize);
  const int bh = 1 << b_height_log2(bsize);
  int i;
  for (i = 0; i < MAX_MB_PLANE; i++) {
    struct macroblockd_plane *const pd = &xd->plane[i];
    vpx_memset(pd->above_context, 0,
               sizeof(ENTROPY_CONTEXT) * (bw >> pd->subsampling_x));
    vpx_memset(pd->left_context, 0,
               sizeof(ENTROPY_CONTEXT) * (bh >> pd->subsampling_y));
  }
}

// This is the index in the scan order beyond which all coefficients for
// 8x8 transform and above are in the top band.
// For 4x4 blocks the index is less but to keep things common the lookup
// table for 4x4 is padded out to this index.
#define MAXBAND_INDEX 21

extern const uint8_t vp9_coefband_trans_8x8plus[MAXBAND_INDEX + 1];
extern const uint8_t vp9_coefband_trans_4x4[MAXBAND_INDEX + 1];


static int get_coef_band(const uint8_t * band_translate, int coef_index) {
  return (coef_index > MAXBAND_INDEX)
    ? (COEF_BANDS-1) : band_translate[coef_index];
}

static INLINE int get_coef_context(const int16_t *neighbors,
                                   uint8_t *token_cache,
                                   int c) {
  return (1 + token_cache[neighbors[MAX_NEIGHBORS * c + 0]] +
          token_cache[neighbors[MAX_NEIGHBORS * c + 1]]) >> 1;
}

const int16_t *vp9_get_coef_neighbors_handle(const int16_t *scan);


// 128 lists of probabilities are stored for the following ONE node probs:
// 1, 3, 5, 7, ..., 253, 255
// In between probabilities are interpolated linearly

#define COEFPROB_MODELS             128

#define UNCONSTRAINED_NODES         3

#define PIVOT_NODE                  2   // which node is pivot

typedef vp9_prob vp9_coeff_probs_model[REF_TYPES][COEF_BANDS]
                                      [PREV_COEF_CONTEXTS]
                                      [UNCONSTRAINED_NODES];

typedef unsigned int vp9_coeff_count_model[REF_TYPES][COEF_BANDS]
                                          [PREV_COEF_CONTEXTS]
                                          [UNCONSTRAINED_NODES + 1];
typedef unsigned int vp9_coeff_stats_model[REF_TYPES][COEF_BANDS]
                                          [PREV_COEF_CONTEXTS]
                                          [UNCONSTRAINED_NODES][2];

void vp9_model_to_full_probs(const vp9_prob *model, vp9_prob *full);

static INLINE const int16_t* get_scan_4x4(TX_TYPE tx_type) {
  switch (tx_type) {
    case ADST_DCT:
      return vp9_row_scan_4x4;
    case DCT_ADST:
      return vp9_col_scan_4x4;
    default:
      return vp9_default_scan_4x4;
  }
}

static INLINE void get_scan_nb_4x4(TX_TYPE tx_type,
                                   const int16_t **scan, const int16_t **nb) {
  switch (tx_type) {
    case ADST_DCT:
      *scan = vp9_row_scan_4x4;
      *nb = vp9_row_scan_4x4_neighbors;
      break;
    case DCT_ADST:
      *scan = vp9_col_scan_4x4;
      *nb = vp9_col_scan_4x4_neighbors;
      break;
    default:
      *scan = vp9_default_scan_4x4;
      *nb = vp9_default_scan_4x4_neighbors;
      break;
  }
}

static INLINE const int16_t* get_iscan_4x4(TX_TYPE tx_type) {
  switch (tx_type) {
    case ADST_DCT:
      return vp9_row_iscan_4x4;
    case DCT_ADST:
      return vp9_col_iscan_4x4;
    default:
      return vp9_default_iscan_4x4;
  }
}

static INLINE const int16_t* get_scan_8x8(TX_TYPE tx_type) {
  switch (tx_type) {
    case ADST_DCT:
      return vp9_row_scan_8x8;
    case DCT_ADST:
      return vp9_col_scan_8x8;
    default:
      return vp9_default_scan_8x8;
  }
}

static INLINE void get_scan_nb_8x8(TX_TYPE tx_type,
                                   const int16_t **scan, const int16_t **nb) {
  switch (tx_type) {
    case ADST_DCT:
      *scan = vp9_row_scan_8x8;
      *nb = vp9_row_scan_8x8_neighbors;
      break;
    case DCT_ADST:
      *scan = vp9_col_scan_8x8;
      *nb = vp9_col_scan_8x8_neighbors;
      break;
    default:
      *scan = vp9_default_scan_8x8;
      *nb = vp9_default_scan_8x8_neighbors;
      break;
  }
}

static INLINE const int16_t* get_iscan_8x8(TX_TYPE tx_type) {
  switch (tx_type) {
    case ADST_DCT:
      return vp9_row_iscan_8x8;
    case DCT_ADST:
      return vp9_col_iscan_8x8;
    default:
      return vp9_default_iscan_8x8;
  }
}

static INLINE const int16_t* get_scan_16x16(TX_TYPE tx_type) {
  switch (tx_type) {
    case ADST_DCT:
      return vp9_row_scan_16x16;
    case DCT_ADST:
      return vp9_col_scan_16x16;
    default:
      return vp9_default_scan_16x16;
  }
}

static INLINE void get_scan_nb_16x16(TX_TYPE tx_type,
                                     const int16_t **scan, const int16_t **nb) {
  switch (tx_type) {
    case ADST_DCT:
      *scan = vp9_row_scan_16x16;
      *nb = vp9_row_scan_16x16_neighbors;
      break;
    case DCT_ADST:
      *scan = vp9_col_scan_16x16;
      *nb = vp9_col_scan_16x16_neighbors;
      break;
    default:
      *scan = vp9_default_scan_16x16;
      *nb = vp9_default_scan_16x16_neighbors;
      break;
  }
}

static INLINE const int16_t* get_iscan_16x16(TX_TYPE tx_type) {
  switch (tx_type) {
    case ADST_DCT:
      return vp9_row_iscan_16x16;
    case DCT_ADST:
      return vp9_col_iscan_16x16;
    default:
      return vp9_default_iscan_16x16;
  }
}

enum { VP9_COEF_UPDATE_PROB = 252 };

#endif  // VP9_COMMON_VP9_ENTROPY_H_
