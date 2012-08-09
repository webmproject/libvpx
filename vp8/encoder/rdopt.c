/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "vp8/common/pragmas.h"

#include "tokenize.h"
#include "treewriter.h"
#include "onyx_int.h"
#include "modecosts.h"
#include "encodeintra.h"
#include "vp8/common/entropymode.h"
#include "vp8/common/reconinter.h"
#include "vp8/common/reconintra.h"
#include "vp8/common/reconintra4x4.h"
#include "vp8/common/findnearmv.h"
#include "vp8/common/quant_common.h"
#include "encodemb.h"
#include "quantize.h"
#include "vp8/common/idct.h"
#include "vp8/common/g_common.h"
#include "variance.h"
#include "mcomp.h"
#include "rdopt.h"
#include "ratectrl.h"
#include "vpx_mem/vpx_mem.h"
#include "dct.h"
#include "vp8/common/systemdependent.h"

#include "vp8/common/seg_common.h"
#include "vp8/common/pred_common.h"

#if CONFIG_RUNTIME_CPU_DETECT
#define IF_RTCD(x)  (x)
#else
#define IF_RTCD(x)  NULL
#endif

extern void vp8cx_mb_init_quantizer(VP8_COMP *cpi, MACROBLOCK *x);
extern void vp8_update_zbin_extra(VP8_COMP *cpi, MACROBLOCK *x);

#if CONFIG_HYBRIDTRANSFORM
extern void vp8_ht_quantize_b(BLOCK *b, BLOCKD *d);
#endif

#define XMVCOST (x->e_mbd.allow_high_precision_mv?x->mvcost_hp:x->mvcost)

#define MAXF(a,b)            (((a) > (b)) ? (a) : (b))

#define INVALID_MV 0x80008000

#if CONFIG_SWITCHABLE_INTERP
/* Factor to weigh the rate for switchable interp filters */
#define SWITCHABLE_INTERP_RATE_FACTOR 1
#endif

static const int auto_speed_thresh[17] = {
  1000,
  200,
  150,
  130,
  150,
  125,
  120,
  115,
  115,
  115,
  115,
  115,
  115,
  115,
  115,
  115,
  105
};

#if CONFIG_PRED_FILTER
const MODE_DEFINITION vp8_mode_order[MAX_MODES] = {
  {ZEROMV,    LAST_FRAME,   0,  0},
  {ZEROMV,    LAST_FRAME,   0,  1},
  {DC_PRED,   INTRA_FRAME,  0,  0},

  {NEARESTMV, LAST_FRAME,   0,  0},
  {NEARESTMV, LAST_FRAME,   0,  1},
  {NEARMV,    LAST_FRAME,   0,  0},
  {NEARMV,    LAST_FRAME,   0,  1},

  {ZEROMV,    GOLDEN_FRAME, 0,  0},
  {ZEROMV,    GOLDEN_FRAME, 0,  1},
  {NEARESTMV, GOLDEN_FRAME, 0,  0},
  {NEARESTMV, GOLDEN_FRAME, 0,  1},

  {ZEROMV,    ALTREF_FRAME, 0,  0},
  {ZEROMV,    ALTREF_FRAME, 0,  1},
  {NEARESTMV, ALTREF_FRAME, 0,  0},
  {NEARESTMV, ALTREF_FRAME, 0,  1},

  {NEARMV,    GOLDEN_FRAME, 0,  0},
  {NEARMV,    GOLDEN_FRAME, 0,  1},
  {NEARMV,    ALTREF_FRAME, 0,  0},
  {NEARMV,    ALTREF_FRAME, 0,  1},

  {V_PRED,    INTRA_FRAME,  0,  0},
  {H_PRED,    INTRA_FRAME,  0,  0},
  {D45_PRED,  INTRA_FRAME,  0,  0},
  {D135_PRED, INTRA_FRAME,  0,  0},
  {D117_PRED, INTRA_FRAME,  0,  0},
  {D153_PRED, INTRA_FRAME,  0,  0},
  {D27_PRED,  INTRA_FRAME,  0,  0},
  {D63_PRED,  INTRA_FRAME,  0,  0},

  {TM_PRED,   INTRA_FRAME,  0,  0},

  {NEWMV,     LAST_FRAME,   0,  0},
  {NEWMV,     LAST_FRAME,   0,  1},
  {NEWMV,     GOLDEN_FRAME, 0,  0},
  {NEWMV,     GOLDEN_FRAME, 0,  1},
  {NEWMV,     ALTREF_FRAME, 0,  0},
  {NEWMV,     ALTREF_FRAME, 0,  1},

  {SPLITMV,   LAST_FRAME,   0,  0},
  {SPLITMV,   GOLDEN_FRAME, 0,  0},
  {SPLITMV,   ALTREF_FRAME, 0,  0},

  {B_PRED,    INTRA_FRAME,  0,  0},
  {I8X8_PRED, INTRA_FRAME,  0,  0},

  /* compound prediction modes */
  {ZEROMV,    LAST_FRAME,   GOLDEN_FRAME, 0},
  {NEARESTMV, LAST_FRAME,   GOLDEN_FRAME, 0},
  {NEARMV,    LAST_FRAME,   GOLDEN_FRAME, 0},

  {ZEROMV,    ALTREF_FRAME, LAST_FRAME,   0},
  {NEARESTMV, ALTREF_FRAME, LAST_FRAME,   0},
  {NEARMV,    ALTREF_FRAME, LAST_FRAME,   0},

  {ZEROMV,    GOLDEN_FRAME, ALTREF_FRAME, 0},
  {NEARESTMV, GOLDEN_FRAME, ALTREF_FRAME, 0},
  {NEARMV,    GOLDEN_FRAME, ALTREF_FRAME, 0},

  {NEWMV,     LAST_FRAME,   GOLDEN_FRAME, 0},
  {NEWMV,     ALTREF_FRAME, LAST_FRAME,   0},
  {NEWMV,     GOLDEN_FRAME, ALTREF_FRAME, 0},

  {SPLITMV,   LAST_FRAME,   GOLDEN_FRAME, 0},
  {SPLITMV,   ALTREF_FRAME, LAST_FRAME,   0},
  {SPLITMV,   GOLDEN_FRAME, ALTREF_FRAME, 0}
};
#else
const MODE_DEFINITION vp8_mode_order[MAX_MODES] = {
  {ZEROMV,    LAST_FRAME,   0},
  {DC_PRED,   INTRA_FRAME,  0},

  {NEARESTMV, LAST_FRAME,   0},
  {NEARMV,    LAST_FRAME,   0},

  {ZEROMV,    GOLDEN_FRAME, 0},
  {NEARESTMV, GOLDEN_FRAME, 0},

  {ZEROMV,    ALTREF_FRAME, 0},
  {NEARESTMV, ALTREF_FRAME, 0},

  {NEARMV,    GOLDEN_FRAME, 0},
  {NEARMV,    ALTREF_FRAME, 0},

  {V_PRED,    INTRA_FRAME,  0},
  {H_PRED,    INTRA_FRAME,  0},
  {D45_PRED,  INTRA_FRAME,  0},
  {D135_PRED, INTRA_FRAME,  0},
  {D117_PRED, INTRA_FRAME,  0},
  {D153_PRED, INTRA_FRAME,  0},
  {D27_PRED,  INTRA_FRAME,  0},
  {D63_PRED,  INTRA_FRAME,  0},

  {TM_PRED,   INTRA_FRAME,  0},

  {NEWMV,     LAST_FRAME,   0},
  {NEWMV,     GOLDEN_FRAME, 0},
  {NEWMV,     ALTREF_FRAME, 0},

  {SPLITMV,   LAST_FRAME,   0},
  {SPLITMV,   GOLDEN_FRAME, 0},
  {SPLITMV,   ALTREF_FRAME, 0},

  {B_PRED,    INTRA_FRAME,  0},
  {I8X8_PRED, INTRA_FRAME,  0},

  /* compound prediction modes */
  {ZEROMV,    LAST_FRAME,   GOLDEN_FRAME},
  {NEARESTMV, LAST_FRAME,   GOLDEN_FRAME},
  {NEARMV,    LAST_FRAME,   GOLDEN_FRAME},

  {ZEROMV,    ALTREF_FRAME, LAST_FRAME},
  {NEARESTMV, ALTREF_FRAME, LAST_FRAME},
  {NEARMV,    ALTREF_FRAME, LAST_FRAME},

  {ZEROMV,    GOLDEN_FRAME, ALTREF_FRAME},
  {NEARESTMV, GOLDEN_FRAME, ALTREF_FRAME},
  {NEARMV,    GOLDEN_FRAME, ALTREF_FRAME},

  {NEWMV,     LAST_FRAME,   GOLDEN_FRAME},
  {NEWMV,     ALTREF_FRAME, LAST_FRAME  },
  {NEWMV,     GOLDEN_FRAME, ALTREF_FRAME},

  {SPLITMV,   LAST_FRAME,   GOLDEN_FRAME},
  {SPLITMV,   ALTREF_FRAME, LAST_FRAME  },
  {SPLITMV,   GOLDEN_FRAME, ALTREF_FRAME}
};
#endif

static void fill_token_costs(
  unsigned int (*c)[COEF_BANDS] [PREV_COEF_CONTEXTS] [MAX_ENTROPY_TOKENS],
  const vp8_prob(*p)[COEF_BANDS] [PREV_COEF_CONTEXTS] [ENTROPY_NODES],
  int block_type_counts) {
  int i, j, k;

  for (i = 0; i < block_type_counts; i++)
    for (j = 0; j < COEF_BANDS; j++)
      for (k = 0; k < PREV_COEF_CONTEXTS; k++) {
        if (k == 0 && ((j > 0 && i > 0) || (j > 1 && i == 0)))
          vp8_cost_tokens_skip((int *)(c [i][j][k]),
                               p [i][j][k],
                               vp8_coef_tree);
        else
          vp8_cost_tokens((int *)(c [i][j][k]),
                          p [i][j][k],
                          vp8_coef_tree);
      }
}


static int rd_iifactor [ 32 ] =  {    4,   4,   3,   2,   1,   0,   0,   0,
                                      0,   0,   0,   0,   0,   0,   0,   0,
                                      0,   0,   0,   0,   0,   0,   0,   0,
                                      0,   0,   0,   0,   0,   0,   0,   0,
                                 };

// 3* dc_qlookup[Q]*dc_qlookup[Q];

/* values are now correlated to quantizer */
static int sad_per_bit16lut[QINDEX_RANGE];
static int sad_per_bit4lut[QINDEX_RANGE];

void vp8_init_me_luts() {
  int i;

  // Initialize the sad lut tables using a formulaic calculation for now
  // This is to make it easier to resolve the impact of experimental changes
  // to the quantizer tables.
  for (i = 0; i < QINDEX_RANGE; i++) {
    sad_per_bit16lut[i] =
      (int)((0.0418 * vp8_convert_qindex_to_q(i)) + 2.4107);
    sad_per_bit4lut[i] = (int)((0.063 * vp8_convert_qindex_to_q(i)) + 2.742);
  }
}

int compute_rd_mult(int qindex) {
  int q;

  q = vp8_dc_quant(qindex, 0);
  return (11 * q * q) >> 6;
}

void vp8cx_initialize_me_consts(VP8_COMP *cpi, int QIndex) {
  cpi->mb.sadperbit16 =  sad_per_bit16lut[QIndex];
  cpi->mb.sadperbit4  =  sad_per_bit4lut[QIndex];
}


void vp8_initialize_rd_consts(VP8_COMP *cpi, int QIndex) {
  int q, i;

  vp8_clear_system_state();  // __asm emms;

  // Further tests required to see if optimum is different
  // for key frames, golden frames and arf frames.
  // if (cpi->common.refresh_golden_frame ||
  //     cpi->common.refresh_alt_ref_frame)
  QIndex = (QIndex < 0) ? 0 : ((QIndex > MAXQ) ? MAXQ : QIndex);

  cpi->RDMULT = compute_rd_mult(QIndex);

  // Extend rate multiplier along side quantizer zbin increases
  if (cpi->zbin_over_quant  > 0) {
    double oq_factor;

    // Experimental code using the same basic equation as used for Q above
    // The units of cpi->zbin_over_quant are 1/128 of Q bin size
    oq_factor = 1.0 + ((double)0.0015625 * cpi->zbin_over_quant);
    cpi->RDMULT = (int)((double)cpi->RDMULT * oq_factor * oq_factor);
  }

  if (cpi->pass == 2 && (cpi->common.frame_type != KEY_FRAME)) {
    if (cpi->twopass.next_iiratio > 31)
      cpi->RDMULT += (cpi->RDMULT * rd_iifactor[31]) >> 4;
    else
      cpi->RDMULT +=
        (cpi->RDMULT * rd_iifactor[cpi->twopass.next_iiratio]) >> 4;
  }

  if (cpi->RDMULT < 7)
    cpi->RDMULT = 7;

  cpi->mb.errorperbit = (cpi->RDMULT / 110);
  cpi->mb.errorperbit += (cpi->mb.errorperbit == 0);

  vp8_set_speed_features(cpi);

  q = (int)pow(vp8_dc_quant(QIndex, 0) >> 2, 1.25);
  q = q << 2;
  cpi->RDMULT = cpi->RDMULT << 4;

  if (q < 8)
    q = 8;

  if (cpi->RDMULT > 1000) {
    cpi->RDDIV = 1;
    cpi->RDMULT /= 100;

    for (i = 0; i < MAX_MODES; i++) {
      if (cpi->sf.thresh_mult[i] < INT_MAX) {
        cpi->rd_threshes[i] = cpi->sf.thresh_mult[i] * q / 100;
      } else {
        cpi->rd_threshes[i] = INT_MAX;
      }

      cpi->rd_baseline_thresh[i] = cpi->rd_threshes[i];
    }
  } else {
    cpi->RDDIV = 100;

    for (i = 0; i < MAX_MODES; i++) {
      if (cpi->sf.thresh_mult[i] < (INT_MAX / q)) {
        cpi->rd_threshes[i] = cpi->sf.thresh_mult[i] * q;
      } else {
        cpi->rd_threshes[i] = INT_MAX;
      }

      cpi->rd_baseline_thresh[i] = cpi->rd_threshes[i];
    }
  }

  fill_token_costs(
    cpi->mb.token_costs[TX_4X4],
    (const vp8_prob( *)[8][PREV_COEF_CONTEXTS][11]) cpi->common.fc.coef_probs,
    BLOCK_TYPES);

  fill_token_costs(
    cpi->mb.token_costs[TX_8X8],
    (const vp8_prob( *)[8][PREV_COEF_CONTEXTS][11]) cpi->common.fc.coef_probs_8x8,
    BLOCK_TYPES_8X8);

#if CONFIG_TX16X16
  fill_token_costs(
    cpi->mb.token_costs[TX_16X16],
    (const vp8_prob(*)[8][PREV_COEF_CONTEXTS][11]) cpi->common.fc.coef_probs_16x16,
    BLOCK_TYPES_16X16);
#endif

  /*rough estimate for costing*/
  cpi->common.kf_ymode_probs_index = cpi->common.base_qindex >> 4;
  vp8_init_mode_costs(cpi);

}

void vp8_auto_select_speed(VP8_COMP *cpi) {
  int milliseconds_for_compress = (int)(1000000 / cpi->oxcf.frame_rate);

  milliseconds_for_compress = milliseconds_for_compress * (16 - cpi->oxcf.cpu_used) / 16;

#if 0

  if (0) {
    FILE *f;

    f = fopen("speed.stt", "a");
    fprintf(f, " %8ld %10ld %10ld %10ld\n",
            cpi->common.current_video_frame, cpi->Speed, milliseconds_for_compress, cpi->avg_pick_mode_time);
    fclose(f);
  }

#endif

  /*
  // this is done during parameter valid check
  if( cpi->oxcf.cpu_used > 16)
      cpi->oxcf.cpu_used = 16;
  if( cpi->oxcf.cpu_used < -16)
      cpi->oxcf.cpu_used = -16;
  */

  if (cpi->avg_pick_mode_time < milliseconds_for_compress && (cpi->avg_encode_time - cpi->avg_pick_mode_time) < milliseconds_for_compress) {
    if (cpi->avg_pick_mode_time == 0) {
      cpi->Speed = 4;
    } else {
      if (milliseconds_for_compress * 100 < cpi->avg_encode_time * 95) {
        cpi->Speed          += 2;
        cpi->avg_pick_mode_time = 0;
        cpi->avg_encode_time = 0;

        if (cpi->Speed > 16) {
          cpi->Speed = 16;
        }
      }

      if (milliseconds_for_compress * 100 > cpi->avg_encode_time * auto_speed_thresh[cpi->Speed]) {
        cpi->Speed          -= 1;
        cpi->avg_pick_mode_time = 0;
        cpi->avg_encode_time = 0;

        // In real-time mode, cpi->speed is in [4, 16].
        if (cpi->Speed < 4) {      // if ( cpi->Speed < 0 )
          cpi->Speed = 4;        // cpi->Speed = 0;
        }
      }
    }
  } else {
    cpi->Speed += 4;

    if (cpi->Speed > 16)
      cpi->Speed = 16;


    cpi->avg_pick_mode_time = 0;
    cpi->avg_encode_time = 0;
  }
}

int vp8_block_error_c(short *coeff, short *dqcoeff, int block_size) {
  int i, error = 0;

  for (i = 0; i < block_size; i++) {
    int this_diff = coeff[i] - dqcoeff[i];
    error += this_diff * this_diff;
  }

  return error;
}

int vp8_mbblock_error_c(MACROBLOCK *mb, int dc) {
  BLOCK  *be;
  BLOCKD *bd;
  int i, j;
  int berror, error = 0;

  for (i = 0; i < 16; i++) {
    be = &mb->block[i];
    bd = &mb->e_mbd.block[i];

    berror = 0;

    for (j = dc; j < 16; j++) {
      int this_diff = be->coeff[j] - bd->dqcoeff[j];
      berror += this_diff * this_diff;
    }

    error += berror;
  }

  return error;
}

int vp8_mbuverror_c(MACROBLOCK *mb) {
  BLOCK  *be;
  BLOCKD *bd;

  int i, error = 0;

  for (i = 16; i < 24; i++) {
    be = &mb->block[i];
    bd = &mb->e_mbd.block[i];

    error += vp8_block_error_c(be->coeff, bd->dqcoeff, 16);
  }

  return error;
}

int VP8_UVSSE(MACROBLOCK *x, const vp8_variance_rtcd_vtable_t *rtcd) {
  unsigned char *uptr, *vptr;
  unsigned char *upred_ptr = (*(x->block[16].base_src) + x->block[16].src);
  unsigned char *vpred_ptr = (*(x->block[20].base_src) + x->block[20].src);
  int uv_stride = x->block[16].src_stride;

  unsigned int sse1 = 0;
  unsigned int sse2 = 0;
  int mv_row = x->e_mbd.mode_info_context->mbmi.mv.as_mv.row;
  int mv_col = x->e_mbd.mode_info_context->mbmi.mv.as_mv.col;
  int offset;
  int pre_stride = x->e_mbd.block[16].pre_stride;

  if (mv_row < 0)
    mv_row -= 1;
  else
    mv_row += 1;

  if (mv_col < 0)
    mv_col -= 1;
  else
    mv_col += 1;

  mv_row /= 2;
  mv_col /= 2;

  offset = (mv_row >> 3) * pre_stride + (mv_col >> 3);
  uptr = x->e_mbd.pre.u_buffer + offset;
  vptr = x->e_mbd.pre.v_buffer + offset;

  if ((mv_row | mv_col) & 7) {
    VARIANCE_INVOKE(rtcd, subpixvar8x8)(uptr, pre_stride,
                                        (mv_col & 7) << 1, (mv_row & 7) << 1, upred_ptr, uv_stride, &sse2);
    VARIANCE_INVOKE(rtcd, subpixvar8x8)(vptr, pre_stride,
                                        (mv_col & 7) << 1, (mv_row & 7) << 1, vpred_ptr, uv_stride, &sse1);
    sse2 += sse1;
  } else {
    VARIANCE_INVOKE(rtcd, var8x8)(uptr, pre_stride,
                                  upred_ptr, uv_stride, &sse2);
    VARIANCE_INVOKE(rtcd, var8x8)(vptr, pre_stride,
                                  vpred_ptr, uv_stride, &sse1);
    sse2 += sse1;
  }
  return sse2;

}

static int cost_coeffs_2x2(MACROBLOCK *mb,
                           BLOCKD *b, int type,
                           ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l) {
  int c = !type;              /* start at coef 0, unless Y with Y2 */
  int eob = b->eob;
  int pt;    /* surrounding block/prev coef predictor */
  int cost = 0;
  short *qcoeff_ptr = b->qcoeff;

  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
  assert(eob <= 4);

  for (; c < eob; c++) {
    int v = qcoeff_ptr[vp8_default_zig_zag1d[c]];
    int t = vp8_dct_value_tokens_ptr[v].Token;
    cost += mb->token_costs[TX_8X8][type][vp8_coef_bands[c]][pt][t];
    cost += vp8_dct_value_cost_ptr[v];
    pt = vp8_prev_token_class[t];
  }

  if (c < 4)
    cost += mb->token_costs[TX_8X8][type][vp8_coef_bands[c]]
            [pt] [DCT_EOB_TOKEN];

  pt = (c != !type); // is eob first coefficient;
  *a = *l = pt;
  return cost;
}

static int cost_coeffs(MACROBLOCK *mb, BLOCKD *b, int type,
                       ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                       int tx_type) {
  const int eob = b->eob;
  int c = !type;              /* start at coef 0, unless Y with Y2 */
  int cost = 0, default_eob, seg_eob;
  int pt;                     /* surrounding block/prev coef predictor */
  int const *scan, *band;
  short *qcoeff_ptr = b->qcoeff;

  int segment_id = mb->e_mbd.mode_info_context->mbmi.segment_id;

  switch (tx_type) {
    case TX_4X4:
      scan = vp8_default_zig_zag1d;
      band = vp8_coef_bands;
      default_eob = 16;
#if CONFIG_HYBRIDTRANSFORM
      {
        int active_ht = (mb->q_index < ACTIVE_HT) &&
                      (mb->e_mbd.mode_info_context->mbmi.mode_rdopt == B_PRED);

        if((type == PLANE_TYPE_Y_WITH_DC) && active_ht) {
          switch (b->bmi.as_mode.tx_type) {
            case ADST_DCT:
              scan = vp8_row_scan;
              break;

            case DCT_ADST:
              scan = vp8_col_scan;
              break;

            default:
              scan = vp8_default_zig_zag1d;
              break;
          }

        } else
          scan = vp8_default_zig_zag1d;
      }
#endif
      break;
    case TX_8X8:
      scan = vp8_default_zig_zag1d_8x8;
      band = vp8_coef_bands_8x8;
      default_eob = 64;
      break;
#if CONFIG_TX16X16
    case TX_16X16:
      scan = vp8_default_zig_zag1d_16x16;
      band = vp8_coef_bands_16x16;
      default_eob = 256;
      break;
#endif
    default:
      break;
  }
  if (segfeature_active(&mb->e_mbd, segment_id, SEG_LVL_EOB))
    seg_eob = get_segdata(&mb->e_mbd, segment_id, SEG_LVL_EOB);
  else
    seg_eob = default_eob;


  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);

  for (; c < eob; c++) {
    int v = qcoeff_ptr[scan[c]];
    int t = vp8_dct_value_tokens_ptr[v].Token;
    cost += mb->token_costs[tx_type][type][band[c]][pt][t];
    cost += vp8_dct_value_cost_ptr[v];
    pt = vp8_prev_token_class[t];
  }

  if (c < seg_eob)
    cost += mb->token_costs[tx_type][type][band[c]]
            [pt][DCT_EOB_TOKEN];

  pt = (c != !type); // is eob first coefficient;
  *a = *l = pt;
  return cost;
}

static int vp8_rdcost_mby(MACROBLOCK *mb) {
  int cost = 0;
  int b;
  MACROBLOCKD *x = &mb->e_mbd;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta;
  ENTROPY_CONTEXT *tl;

  vpx_memcpy(&t_above, mb->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, mb->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;

  for (b = 0; b < 16; b++)
    cost += cost_coeffs(mb, x->block + b, PLANE_TYPE_Y_NO_DC,
                        ta + vp8_block2above[b], tl + vp8_block2left[b],
                        TX_4X4);

  cost += cost_coeffs(mb, x->block + 24, PLANE_TYPE_Y2,
                      ta + vp8_block2above[24], tl + vp8_block2left[24],
                      TX_4X4);

  return cost;
}

static void macro_block_yrd(MACROBLOCK *mb,
                            int *Rate,
                            int *Distortion,
                            const VP8_ENCODER_RTCD *rtcd) {
  int b;
  MACROBLOCKD *const x = &mb->e_mbd;
  BLOCK   *const mb_y2 = mb->block + 24;
  BLOCKD *const x_y2  = x->block + 24;
  short *Y2DCPtr = mb_y2->src_diff;
  BLOCK *beptr;
  int d;

  ENCODEMB_INVOKE(&rtcd->encodemb, submby)(
    mb->src_diff,
    *(mb->block[0].base_src),
    mb->e_mbd.predictor,
    mb->block[0].src_stride);

  // Fdct and building the 2nd order block
  for (beptr = mb->block; beptr < mb->block + 16; beptr += 2) {
    mb->vp8_short_fdct8x4(beptr->src_diff, beptr->coeff, 32);
    *Y2DCPtr++ = beptr->coeff[0];
    *Y2DCPtr++ = beptr->coeff[16];
  }

  // 2nd order fdct
  mb->short_walsh4x4(mb_y2->src_diff, mb_y2->coeff, 8);

  // Quantization
  for (b = 0; b < 16; b++) {
    mb->quantize_b(&mb->block[b], &mb->e_mbd.block[b]);
  }

  // DC predication and Quantization of 2nd Order block
  mb->quantize_b(mb_y2, x_y2);

  // Distortion
  d = ENCODEMB_INVOKE(&rtcd->encodemb, mberr)(mb, 1);

  d += ENCODEMB_INVOKE(&rtcd->encodemb, berr)(mb_y2->coeff, x_y2->dqcoeff, 16);

  *Distortion = (d >> 2);
  // rate
  *Rate = vp8_rdcost_mby(mb);
}

static int vp8_rdcost_mby_8x8(MACROBLOCK *mb) {
  int cost = 0;
  int b;
  MACROBLOCKD *x = &mb->e_mbd;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta;
  ENTROPY_CONTEXT *tl;

  vpx_memcpy(&t_above, mb->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, mb->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;

  for (b = 0; b < 16; b += 4)
    cost += cost_coeffs(mb, x->block + b, PLANE_TYPE_Y_NO_DC,
                        ta + vp8_block2above_8x8[b], tl + vp8_block2left_8x8[b],
                        TX_8X8);

  cost += cost_coeffs_2x2(mb, x->block + 24, PLANE_TYPE_Y2,
                          ta + vp8_block2above[24], tl + vp8_block2left[24]);
  return cost;
}

static void macro_block_yrd_8x8(MACROBLOCK *mb,
                                int *Rate,
                                int *Distortion,
                                const VP8_ENCODER_RTCD *rtcd) {
  MACROBLOCKD *const x = &mb->e_mbd;
  BLOCK   *const mb_y2 = mb->block + 24;
  BLOCKD *const x_y2  = x->block + 24;
  int d;

  ENCODEMB_INVOKE(&rtcd->encodemb, submby)(
    mb->src_diff,
    *(mb->block[0].base_src),
    mb->e_mbd.predictor,
    mb->block[0].src_stride);

  vp8_transform_mby_8x8(mb);
  vp8_quantize_mby_8x8(mb);

  /* remove 1st order dc to properly combine 1st/2nd order distortion */
  mb->coeff[0] = 0;
  mb->coeff[64] = 0;
  mb->coeff[128] = 0;
  mb->coeff[192] = 0;
  mb->e_mbd.dqcoeff[0] = 0;
  mb->e_mbd.dqcoeff[64] = 0;
  mb->e_mbd.dqcoeff[128] = 0;
  mb->e_mbd.dqcoeff[192] = 0;

  d = ENCODEMB_INVOKE(&rtcd->encodemb, mberr)(mb, 0);
  d += ENCODEMB_INVOKE(&rtcd->encodemb, berr)(mb_y2->coeff, x_y2->dqcoeff, 16);

  *Distortion = (d >> 2);
  // rate
  *Rate = vp8_rdcost_mby_8x8(mb);
}

#if CONFIG_TX16X16
static int vp8_rdcost_mby_16x16(MACROBLOCK *mb) {
  int cost;
  MACROBLOCKD *x = &mb->e_mbd;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta, *tl;

  vpx_memcpy(&t_above, mb->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, mb->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;

  cost = cost_coeffs(mb, x->block, PLANE_TYPE_Y_WITH_DC, ta, tl, TX_16X16);
  return cost;
}

static void macro_block_yrd_16x16(MACROBLOCK *mb, int *Rate, int *Distortion,
                                  const VP8_ENCODER_RTCD *rtcd) {
  int d;

  ENCODEMB_INVOKE(&rtcd->encodemb, submby)(
    mb->src_diff,
    *(mb->block[0].base_src),
    mb->e_mbd.predictor,
    mb->block[0].src_stride);

  vp8_transform_mby_16x16(mb);
  vp8_quantize_mby_16x16(mb);
  d = ENCODEMB_INVOKE(&rtcd->encodemb, mberr)(mb, 0);

  *Distortion = (d >> 2);
  // rate
  *Rate = vp8_rdcost_mby_16x16(mb);
}
#endif

static void copy_predictor(unsigned char *dst, const unsigned char *predictor) {
  const unsigned int *p = (const unsigned int *)predictor;
  unsigned int *d = (unsigned int *)dst;
  d[0] = p[0];
  d[4] = p[4];
  d[8] = p[8];
  d[12] = p[12];
}

static void copy_predictor_8x8(unsigned char *dst, const unsigned char *predictor) {
  const unsigned int *p = (const unsigned int *)predictor;
  unsigned int *d = (unsigned int *)dst;
  d[0] = p[0];
  d[1] = p[1];
  d[4] = p[4];
  d[5] = p[5];
  d[8] = p[8];
  d[9] = p[9];
  d[12] = p[12];
  d[13] = p[13];
  d[16] = p[16];
  d[17] = p[17];
  d[20] = p[20];
  d[21] = p[21];
  d[24] = p[24];
  d[25] = p[25];
  d[28] = p[28];
  d[29] = p[29];
}

static int64_t rd_pick_intra4x4block(VP8_COMP *cpi, MACROBLOCK *x, BLOCK *be,
                                     BLOCKD *b, B_PREDICTION_MODE *best_mode,
#if CONFIG_COMP_INTRA_PRED
                                     B_PREDICTION_MODE *best_second_mode,
                                     int allow_comp,
#endif
                                     int *bmode_costs,
                                     ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                                     int *bestrate, int *bestratey,
                                     int *bestdistortion) {
  B_PREDICTION_MODE mode;

#if CONFIG_HYBRIDTRANSFORM
  int QIndex = x->q_index;
  int active_ht = (QIndex < ACTIVE_HT);
  TX_TYPE best_tx_type;
#endif

#if CONFIG_COMP_INTRA_PRED
  B_PREDICTION_MODE mode2;
#endif
  int64_t best_rd = INT64_MAX;
  int rate = 0;
  int distortion;

  ENTROPY_CONTEXT ta = *a, tempa = *a;
  ENTROPY_CONTEXT tl = *l, templ = *l;
  /*
   * The predictor buffer is a 2d buffer with a stride of 16.  Create
   * a temp buffer that meets the stride requirements, but we are only
   * interested in the left 4x4 block
   * */
  DECLARE_ALIGNED_ARRAY(16, unsigned char,  best_predictor, 16 * 4);
  DECLARE_ALIGNED_ARRAY(16, short, best_dqcoeff, 16);

  for (mode = B_DC_PRED; mode <= B_HU_PRED; mode++) {
#if CONFIG_COMP_INTRA_PRED
    for (mode2 = (allow_comp ? 0 : (B_DC_PRED - 1));
                   mode2 != (allow_comp ? (mode + 1) : 0); mode2++) {
#endif
      int64_t this_rd;
      int ratey;

      // TODO Temporarily ignore modes that need the above-right data. SB
      // encoding means this data is not available for the bottom right MB
      // Do we need to do this for mode2 also?
      if (mode == B_LD_PRED || mode == B_VL_PRED)
        continue;
      rate = bmode_costs[mode];

#if CONFIG_COMP_INTRA_PRED
      if (mode2 == (B_PREDICTION_MODE)(B_DC_PRED - 1)) {
#endif
        RECON_INVOKE(&cpi->rtcd.common->recon, intra4x4_predict)
        (b, mode, b->predictor);
#if CONFIG_COMP_INTRA_PRED
      } else {
        RECON_INVOKE(&cpi->rtcd.common->recon, comp_intra4x4_predict)
        (b, mode, mode2, b->predictor);
        rate += bmode_costs[mode2];
      }
#endif
      ENCODEMB_INVOKE(IF_RTCD(&cpi->rtcd.encodemb), subb)(be, b, 16);

#if CONFIG_HYBRIDTRANSFORM
      if (active_ht) {
        b->bmi.as_mode.test = mode;
        txfm_map(b, mode);
        vp8_fht_c(be->src_diff, be->coeff, 32, b->bmi.as_mode.tx_type, 4);
        vp8_ht_quantize_b(be, b);
      } else {
        x->vp8_short_fdct4x4(be->src_diff, be->coeff, 32);
        x->quantize_b(be, b);
      }
#else
        x->vp8_short_fdct4x4(be->src_diff, be->coeff, 32);
        x->quantize_b(be, b);
#endif

        tempa = ta;
        templ = tl;

        ratey = cost_coeffs(x, b, PLANE_TYPE_Y_WITH_DC, &tempa, &templ, TX_4X4);
        rate += ratey;
        distortion = ENCODEMB_INVOKE(IF_RTCD(&cpi->rtcd.encodemb), berr)(
            be->coeff, b->dqcoeff, 16) >> 2;

        this_rd = RDCOST(x->rdmult, x->rddiv, rate, distortion);

        if (this_rd < best_rd) {
          *bestrate = rate;
          *bestratey = ratey;
          *bestdistortion = distortion;
          best_rd = this_rd;
          *best_mode = mode;
#if CONFIG_HYBRIDTRANSFORM
          best_tx_type = b->bmi.as_mode.tx_type ;
#endif

#if CONFIG_COMP_INTRA_PRED
          *best_second_mode = mode2;
#endif
          *a = tempa;
          *l = templ;
          copy_predictor(best_predictor, b->predictor);
          vpx_memcpy(best_dqcoeff, b->dqcoeff, 32);
        }
#if CONFIG_COMP_INTRA_PRED
    }
#endif
  }
  b->bmi.as_mode.first = (B_PREDICTION_MODE)(*best_mode);
#if CONFIG_COMP_INTRA_PRED
  b->bmi.as_mode.second = (B_PREDICTION_MODE)(*best_second_mode);
#endif

#if CONFIG_HYBRIDTRANSFORM
  b->bmi.as_mode.tx_type = best_tx_type;

  // inverse transform
  if (active_ht)
    vp8_ihtllm_c(best_dqcoeff, b->diff, 32, b->bmi.as_mode.tx_type, 4);
  else
    IDCT_INVOKE(IF_RTCD(&cpi->rtcd.common->idct), idct16)(best_dqcoeff,
                                                                b->diff, 32);
#else
  IDCT_INVOKE(IF_RTCD(&cpi->rtcd.common->idct), idct16)(best_dqcoeff, b->diff, 32);
#endif

  RECON_INVOKE(IF_RTCD(&cpi->rtcd.common->recon), recon)(best_predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);

  return best_rd;
}

static int64_t rd_pick_intra4x4mby_modes(VP8_COMP *cpi, MACROBLOCK *mb, int *Rate,
                                     int *rate_y, int *Distortion, int64_t best_rd,
#if CONFIG_COMP_INTRA_PRED
                                     int allow_comp,
#endif
                                     int update_contexts) {
  int i;
  MACROBLOCKD *const xd = &mb->e_mbd;
  int cost = mb->mbmode_cost [xd->frame_type] [B_PRED];
  int distortion = 0;
  int tot_rate_y = 0;
  int64_t total_rd = 0;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta, *tl;
  int *bmode_costs;

  if (update_contexts) {
    ta = (ENTROPY_CONTEXT *)mb->e_mbd.above_context;
    tl = (ENTROPY_CONTEXT *)mb->e_mbd.left_context;
  } else {
    vpx_memcpy(&t_above, mb->e_mbd.above_context,
               sizeof(ENTROPY_CONTEXT_PLANES));
    vpx_memcpy(&t_left, mb->e_mbd.left_context,
               sizeof(ENTROPY_CONTEXT_PLANES));

    ta = (ENTROPY_CONTEXT *)&t_above;
    tl = (ENTROPY_CONTEXT *)&t_left;
  }

  // TODO(agrange)
  // vp8_intra_prediction_down_copy(xd);

  bmode_costs = mb->inter_bmode_costs;

  for (i = 0; i < 16; i++) {
    MODE_INFO *const mic = xd->mode_info_context;
    const int mis = xd->mode_info_stride;
    B_PREDICTION_MODE UNINITIALIZED_IS_SAFE(best_mode);
#if CONFIG_COMP_INTRA_PRED
    B_PREDICTION_MODE UNINITIALIZED_IS_SAFE(best_second_mode);
#endif
    int UNINITIALIZED_IS_SAFE(r), UNINITIALIZED_IS_SAFE(ry), UNINITIALIZED_IS_SAFE(d);

    if (mb->e_mbd.frame_type == KEY_FRAME) {
      const B_PREDICTION_MODE A = above_block_mode(mic, i, mis);
      const B_PREDICTION_MODE L = left_block_mode(mic, i);

      bmode_costs  = mb->bmode_costs[A][L];
    }

    total_rd += rd_pick_intra4x4block(
                  cpi, mb, mb->block + i, xd->block + i, &best_mode,
#if CONFIG_COMP_INTRA_PRED
                  & best_second_mode, allow_comp,
#endif
                  bmode_costs, ta + vp8_block2above[i],
                  tl + vp8_block2left[i], &r, &ry, &d);

    cost += r;
    distortion += d;
    tot_rate_y += ry;

    mic->bmi[i].as_mode.first = best_mode;
#if CONFIG_COMP_INTRA_PRED
    mic->bmi[i].as_mode.second = best_second_mode;
#endif

    if (total_rd >= best_rd)
      break;
  }

  if (total_rd >= best_rd)
    return INT64_MAX;

#if CONFIG_COMP_INTRA_PRED
  cost += vp8_cost_bit(128, allow_comp);
#endif
  *Rate = cost;
  *rate_y += tot_rate_y;
  *Distortion = distortion;

  return RDCOST(mb->rdmult, mb->rddiv, cost, distortion);
}


static int64_t rd_pick_intra16x16mby_mode(VP8_COMP *cpi,
                                      MACROBLOCK *x,
                                      int *Rate,
                                      int *rate_y,
                                      int *Distortion) {
  MB_PREDICTION_MODE mode;
  MB_PREDICTION_MODE UNINITIALIZED_IS_SAFE(mode_selected);
#if CONFIG_COMP_INTRA_PRED
  MB_PREDICTION_MODE mode2;
  MB_PREDICTION_MODE UNINITIALIZED_IS_SAFE(mode2_selected);
#endif
  int rate, ratey;
  int distortion;
  int64_t best_rd = INT64_MAX;
  int64_t this_rd;

  // Y Search for 16x16 intra prediction mode
  for (mode = DC_PRED; mode <= TM_PRED; mode++) {
    x->e_mbd.mode_info_context->mbmi.mode = mode;
#if CONFIG_COMP_INTRA_PRED
    for (mode2 = DC_PRED - 1; mode2 != TM_PRED + 1; mode2++) {
      x->e_mbd.mode_info_context->mbmi.second_mode = mode2;
      if (mode2 == (MB_PREDICTION_MODE)(DC_PRED - 1)) {
#endif
        RECON_INVOKE(&cpi->common.rtcd.recon, build_intra_predictors_mby)
        (&x->e_mbd);
#if CONFIG_COMP_INTRA_PRED
      } else {
        continue; // i.e. disable for now
        RECON_INVOKE(&cpi->common.rtcd.recon, build_comp_intra_predictors_mby)(&x->e_mbd);
      }
#endif

#if CONFIG_TX16X16
      macro_block_yrd_16x16(x, &ratey, &distortion, IF_RTCD(&cpi->rtcd));
#else
      macro_block_yrd_8x8(x, &ratey, &distortion, IF_RTCD(&cpi->rtcd));
#endif
      // FIXME add compoundmode cost
      // FIXME add rate for mode2
      rate = ratey + x->mbmode_cost[x->e_mbd.frame_type]
             [x->e_mbd.mode_info_context->mbmi.mode];

      this_rd = RDCOST(x->rdmult, x->rddiv, rate, distortion);

      if (this_rd < best_rd) {
        mode_selected = mode;
#if CONFIG_COMP_INTRA_PRED
        mode2_selected = mode2;
#endif
        best_rd = this_rd;
        *Rate = rate;
        *rate_y = ratey;
        *Distortion = distortion;
      }
#if CONFIG_COMP_INTRA_PRED
    }
#endif
  }

  x->e_mbd.mode_info_context->mbmi.mode = mode_selected;
#if CONFIG_COMP_INTRA_PRED
  x->e_mbd.mode_info_context->mbmi.second_mode = mode2_selected;
#endif
  return best_rd;
}


static int64_t rd_pick_intra8x8block(VP8_COMP *cpi, MACROBLOCK *x, int ib,
                                     B_PREDICTION_MODE *best_mode,
#if CONFIG_COMP_INTRA_PRED
                                     B_PREDICTION_MODE *best_second_mode,
#endif
                                     int *mode_costs,
                                     ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                                     int *bestrate, int *bestratey,
                                     int *bestdistortion) {
  MB_PREDICTION_MODE mode;
#if CONFIG_COMP_INTRA_PRED
  MB_PREDICTION_MODE mode2;
#endif
  MACROBLOCKD *xd = &x->e_mbd;
  int64_t best_rd = INT64_MAX;
  int distortion, rate = 0;
  BLOCK  *be = x->block + ib;
  BLOCKD *b = x->e_mbd.block + ib;
  ENTROPY_CONTEXT ta0, ta1, besta0 = 0, besta1 = 0;
  ENTROPY_CONTEXT tl0, tl1, bestl0 = 0, bestl1 = 0;

  /*
   * The predictor buffer is a 2d buffer with a stride of 16.  Create
   * a temp buffer that meets the stride requirements, but we are only
   * interested in the left 8x8 block
   * */
  DECLARE_ALIGNED_ARRAY(16, unsigned char,  best_predictor, 16 * 8);
  DECLARE_ALIGNED_ARRAY(16, short, best_dqcoeff, 16 * 4);

#if CONFIG_HYBRIDTRANSFORM8X8
  // perform transformation of dimension 8x8
  // note the input and output index mapping
  int idx = (ib & 0x02) ? (ib + 2) : ib;
#endif

  for (mode = DC_PRED; mode <= TM_PRED; mode++) {
#if CONFIG_COMP_INTRA_PRED
    for (mode2 = DC_PRED - 1; mode2 != TM_PRED + 1; mode2++) {
#endif
      int64_t this_rd;
      int rate_t;

      // FIXME rate for compound mode and second intrapred mode
      rate = mode_costs[mode];

#if CONFIG_COMP_INTRA_PRED
      if (mode2 == (MB_PREDICTION_MODE)(DC_PRED - 1)) {
#endif
        RECON_INVOKE(&cpi->rtcd.common->recon, intra8x8_predict)
        (b, mode, b->predictor);
#if CONFIG_COMP_INTRA_PRED
      } else {
        continue; // i.e. disable for now
        RECON_INVOKE(&cpi->rtcd.common->recon, comp_intra8x8_predict)
        (b, mode, mode2, b->predictor);
      }
#endif

      vp8_subtract_4b_c(be, b, 16);

#if CONFIG_HYBRIDTRANSFORM8X8
      txfm_map(b, pred_mode_conv(mode));
      vp8_fht_c(be->src_diff, (x->block + idx)->coeff, 32,
                b->bmi.as_mode.tx_type, 8);
      x->quantize_b_8x8(x->block + idx, xd->block + idx);

      // compute quantization mse of 8x8 block
      distortion = vp8_block_error_c((x->block + idx)->coeff,
                                     (xd->block + idx)->dqcoeff, 64)>>2;

      ta0 = *(a + vp8_block2above_8x8[idx]);
      tl0 = *(l + vp8_block2left_8x8 [idx]);

      rate_t = cost_coeffs(x, xd->block + idx, PLANE_TYPE_Y_WITH_DC,
                           &ta0, &tl0, TX_8X8);
      rate += rate_t;
      ta1 = ta0;
      tl1 = tl0;
#else
      x->vp8_short_fdct8x4(be->src_diff, be->coeff, 32);
      x->vp8_short_fdct8x4(be->src_diff + 64, be->coeff + 64, 32);

      x->quantize_b_pair(x->block + ib, x->block + ib + 1,
                         xd->block + ib, xd->block + ib + 1);
      x->quantize_b_pair(x->block + ib + 4, x->block + ib + 5,
                         xd->block + ib + 4, xd->block + ib + 5);

      distortion = ENCODEMB_INVOKE(IF_RTCD(&cpi->rtcd.encodemb), berr)
                   ((x->block + ib)->coeff, (xd->block + ib)->dqcoeff, 16) >> 2;
      distortion += ENCODEMB_INVOKE(IF_RTCD(&cpi->rtcd.encodemb), berr)
                    ((x->block + ib + 1)->coeff, (xd->block + ib + 1)->dqcoeff, 16) >> 2;
      distortion += ENCODEMB_INVOKE(IF_RTCD(&cpi->rtcd.encodemb), berr)
                    ((x->block + ib + 4)->coeff, (xd->block + ib + 4)->dqcoeff, 16) >> 2;
      distortion += ENCODEMB_INVOKE(IF_RTCD(&cpi->rtcd.encodemb), berr)
                    ((x->block + ib + 5)->coeff, (xd->block + ib + 5)->dqcoeff, 16) >> 2;

      ta0 = *(a + vp8_block2above[ib]);
      ta1 = *(a + vp8_block2above[ib + 1]);
      tl0 = *(l + vp8_block2above[ib]);
      tl1 = *(l + vp8_block2above[ib + 4]);
      rate_t = cost_coeffs(x, xd->block + ib, PLANE_TYPE_Y_WITH_DC,
                           &ta0, &tl0, TX_4X4);
      rate_t += cost_coeffs(x, xd->block + ib + 1, PLANE_TYPE_Y_WITH_DC,
                            &ta1, &tl0, TX_4X4);
      rate_t += cost_coeffs(x, xd->block + ib + 4, PLANE_TYPE_Y_WITH_DC,
                            &ta0, &tl1, TX_4X4);
      rate_t += cost_coeffs(x, xd->block + ib + 5, PLANE_TYPE_Y_WITH_DC,
                            &ta1, &tl1, TX_4X4);
      rate += rate_t;
#endif

      this_rd = RDCOST(x->rdmult, x->rddiv, rate, distortion);
      if (this_rd < best_rd) {
        *bestrate = rate;
        *bestratey = rate_t;
        *bestdistortion = distortion;
        besta0 = ta0;
        besta1 = ta1;
        bestl0 = tl0;
        bestl1 = tl1;
        best_rd = this_rd;
        *best_mode = mode;
#if CONFIG_COMP_INTRA_PRED
        *best_second_mode = mode2;
#endif
        copy_predictor_8x8(best_predictor, b->predictor);
        vpx_memcpy(best_dqcoeff, b->dqcoeff, 64);
        vpx_memcpy(best_dqcoeff + 32, b->dqcoeff + 64, 64);
#if CONFIG_COMP_INTRA_PRED
      }
#endif
    }
  }
  b->bmi.as_mode.first = (*best_mode);
#if CONFIG_COMP_INTRA_PRED
  b->bmi.as_mode.second = (*best_second_mode);
#endif
  vp8_encode_intra8x8(IF_RTCD(&cpi->rtcd), x, ib);

#if CONFIG_HYBRIDTRANSFORM8X8
  *(a + vp8_block2above_8x8[idx])     = besta0;
  *(a + vp8_block2above_8x8[idx] + 1) = besta1;
  *(l + vp8_block2left_8x8 [idx])     = bestl0;
  *(l + vp8_block2left_8x8 [idx] + 1) = bestl1;
#else
  *(a + vp8_block2above[ib])   = besta0;
  *(a + vp8_block2above[ib + 1]) = besta1;
  *(l + vp8_block2above[ib])   = bestl0;
  *(l + vp8_block2above[ib + 4]) = bestl1;
#endif
  return best_rd;
}

const int vp8_i8x8_block[4] = {0, 2, 8, 10};
int64_t rd_pick_intra8x8mby_modes(VP8_COMP *cpi, MACROBLOCK *mb,
                                  int *Rate, int *rate_y,
                                  int *Distortion, int64_t best_rd) {
  MACROBLOCKD *const xd = &mb->e_mbd;
  int i, ib;
  int cost = mb->mbmode_cost [xd->frame_type] [I8X8_PRED];
  int distortion = 0;
  int tot_rate_y = 0;
  long long total_rd = 0;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta, *tl;
  int *i8x8mode_costs;

  vpx_memcpy(&t_above, mb->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, mb->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;

  i8x8mode_costs  = mb->i8x8_mode_costs;

  for (i = 0; i < 4; i++) {
    MODE_INFO *const mic = xd->mode_info_context;
    B_PREDICTION_MODE UNINITIALIZED_IS_SAFE(best_mode);
#if CONFIG_COMP_INTRA_PRED
    B_PREDICTION_MODE UNINITIALIZED_IS_SAFE(best_second_mode);
#endif
    int UNINITIALIZED_IS_SAFE(r), UNINITIALIZED_IS_SAFE(ry), UNINITIALIZED_IS_SAFE(d);

    ib = vp8_i8x8_block[i];
    total_rd += rd_pick_intra8x8block(
                  cpi, mb, ib, &best_mode,
#if CONFIG_COMP_INTRA_PRED
                  & best_second_mode,
#endif
                  i8x8mode_costs, ta, tl, &r, &ry, &d);
    cost += r;
    distortion += d;
    tot_rate_y += ry;
    mic->bmi[ib].as_mode.first = best_mode;
#if CONFIG_COMP_INTRA_PRED
    mic->bmi[ib].as_mode.second = best_second_mode;
#endif
  }
  *Rate = cost;
  *rate_y += tot_rate_y;
  *Distortion = distortion;
  return RDCOST(mb->rdmult, mb->rddiv, cost, distortion);
}

static int rd_cost_mbuv(MACROBLOCK *mb) {
  int b;
  int cost = 0;
  MACROBLOCKD *x = &mb->e_mbd;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta, *tl;

  vpx_memcpy(&t_above, mb->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, mb->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;

  for (b = 16; b < 24; b++)
    cost += cost_coeffs(mb, x->block + b, PLANE_TYPE_UV,
                        ta + vp8_block2above[b], tl + vp8_block2left[b],
                        TX_4X4);

  return cost;
}


static int64_t rd_inter16x16_uv(VP8_COMP *cpi, MACROBLOCK *x, int *rate,
                            int *distortion, int fullpixel) {
  ENCODEMB_INVOKE(IF_RTCD(&cpi->rtcd.encodemb), submbuv)(x->src_diff,
                                                         x->src.u_buffer, x->src.v_buffer, x->e_mbd.predictor, x->src.uv_stride);

  vp8_transform_mbuv(x);
  vp8_quantize_mbuv(x);

  *rate       = rd_cost_mbuv(x);
  *distortion = ENCODEMB_INVOKE(&cpi->rtcd.encodemb, mbuverr)(x) / 4;

  return RDCOST(x->rdmult, x->rddiv, *rate, *distortion);
}

static int rd_cost_mbuv_8x8(MACROBLOCK *mb) {
  int b;
  int cost = 0;
  MACROBLOCKD *x = &mb->e_mbd;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta, *tl;

  vpx_memcpy(&t_above, mb->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, mb->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;

  for (b = 16; b < 24; b += 4)
    cost += cost_coeffs(mb, x->block + b, PLANE_TYPE_UV,
                        ta + vp8_block2above_8x8[b],
                        tl + vp8_block2left_8x8[b], TX_8X8);

  return cost;
}


static int64_t rd_inter16x16_uv_8x8(VP8_COMP *cpi, MACROBLOCK *x, int *rate,
                                    int *distortion, int fullpixel) {
  ENCODEMB_INVOKE(IF_RTCD(&cpi->rtcd.encodemb), submbuv)(x->src_diff,
                                                         x->src.u_buffer, x->src.v_buffer, x->e_mbd.predictor, x->src.uv_stride);

  vp8_transform_mbuv_8x8(x);

  vp8_quantize_mbuv_8x8(x);

  *rate       = rd_cost_mbuv_8x8(x);
  *distortion = ENCODEMB_INVOKE(&cpi->rtcd.encodemb, mbuverr)(x) / 4;

  return RDCOST(x->rdmult, x->rddiv, *rate, *distortion);
}


static int64_t rd_inter4x4_uv(VP8_COMP *cpi, MACROBLOCK *x, int *rate,
                              int *distortion, int fullpixel) {
  vp8_build_inter4x4_predictors_mbuv(&x->e_mbd);
  ENCODEMB_INVOKE(IF_RTCD(&cpi->rtcd.encodemb), submbuv)(x->src_diff,
                                                         x->src.u_buffer, x->src.v_buffer, x->e_mbd.predictor, x->src.uv_stride);

  vp8_transform_mbuv(x);
  vp8_quantize_mbuv(x);

  *rate       = rd_cost_mbuv(x);
  *distortion = ENCODEMB_INVOKE(&cpi->rtcd.encodemb, mbuverr)(x) / 4;

  return RDCOST(x->rdmult, x->rddiv, *rate, *distortion);
}

static void rd_pick_intra_mbuv_mode(VP8_COMP *cpi,
                                    MACROBLOCK *x,
                                    int *rate,
                                    int *rate_tokenonly,
                                    int *distortion) {
  MB_PREDICTION_MODE mode;
  MB_PREDICTION_MODE UNINITIALIZED_IS_SAFE(mode_selected);
#if CONFIG_COMP_INTRA_PRED
  MB_PREDICTION_MODE mode2;
  MB_PREDICTION_MODE UNINITIALIZED_IS_SAFE(mode2_selected);
#endif
  int64_t best_rd = INT64_MAX;
  int UNINITIALIZED_IS_SAFE(d), UNINITIALIZED_IS_SAFE(r);
  int rate_to;

  for (mode = DC_PRED; mode <= TM_PRED; mode++) {
#if CONFIG_COMP_INTRA_PRED
    for (mode2 = DC_PRED - 1; mode2 != TM_PRED + 1; mode2++) {
#endif
      int rate;
      int distortion;
      int64_t this_rd;

      x->e_mbd.mode_info_context->mbmi.uv_mode = mode;
#if CONFIG_COMP_INTRA_PRED
      x->e_mbd.mode_info_context->mbmi.second_uv_mode = mode2;
      if (mode2 == (MB_PREDICTION_MODE)(DC_PRED - 1)) {
#endif
        RECON_INVOKE(&cpi->rtcd.common->recon, build_intra_predictors_mbuv)
        (&x->e_mbd);
#if CONFIG_COMP_INTRA_PRED
      } else {
        continue;
        RECON_INVOKE(&cpi->rtcd.common->recon, build_comp_intra_predictors_mbuv)
        (&x->e_mbd);
      }
#endif

      ENCODEMB_INVOKE(IF_RTCD(&cpi->rtcd.encodemb), submbuv)(x->src_diff,
                                                             x->src.u_buffer, x->src.v_buffer, x->e_mbd.predictor,
                                                             x->src.uv_stride);
      vp8_transform_mbuv(x);
      vp8_quantize_mbuv(x);

      rate_to = rd_cost_mbuv(x);
      rate = rate_to
             + x->intra_uv_mode_cost[x->e_mbd.frame_type]
             [x->e_mbd.mode_info_context->mbmi.uv_mode];

      distortion = ENCODEMB_INVOKE(&cpi->rtcd.encodemb, mbuverr)(x) / 4;

      this_rd = RDCOST(x->rdmult, x->rddiv, rate, distortion);

      if (this_rd < best_rd) {
        best_rd = this_rd;
        d = distortion;
        r = rate;
        *rate_tokenonly = rate_to;
        mode_selected = mode;
#if CONFIG_COMP_INTRA_PRED
        mode2_selected = mode2;
      }
#endif
    }
  }

  *rate = r;
  *distortion = d;

  x->e_mbd.mode_info_context->mbmi.uv_mode = mode_selected;
#if CONFIG_COMP_INTRA_PRED
  x->e_mbd.mode_info_context->mbmi.second_uv_mode = mode2_selected;
#endif
}

static void rd_pick_intra_mbuv_mode_8x8(VP8_COMP *cpi,
                                        MACROBLOCK *x,
                                        int *rate,
                                        int *rate_tokenonly,
                                        int *distortion) {
  MB_PREDICTION_MODE mode;
  MB_PREDICTION_MODE UNINITIALIZED_IS_SAFE(mode_selected);
  int64_t best_rd = INT64_MAX;
  int UNINITIALIZED_IS_SAFE(d), UNINITIALIZED_IS_SAFE(r);
  int rate_to;

  for (mode = DC_PRED; mode <= TM_PRED; mode++) {
    int rate;
    int distortion;
    int64_t this_rd;

    x->e_mbd.mode_info_context->mbmi.uv_mode = mode;
    RECON_INVOKE(&cpi->rtcd.common->recon, build_intra_predictors_mbuv)
    (&x->e_mbd);
    ENCODEMB_INVOKE(IF_RTCD(&cpi->rtcd.encodemb), submbuv)(x->src_diff,
                                                           x->src.u_buffer, x->src.v_buffer, x->e_mbd.predictor,
                                                           x->src.uv_stride);
    vp8_transform_mbuv_8x8(x);

    vp8_quantize_mbuv_8x8(x);

    rate_to = rd_cost_mbuv_8x8(x);
    rate = rate_to + x->intra_uv_mode_cost[x->e_mbd.frame_type]
           [x->e_mbd.mode_info_context->mbmi.uv_mode];

    distortion = ENCODEMB_INVOKE(&cpi->rtcd.encodemb, mbuverr)(x) / 4;
    this_rd = RDCOST(x->rdmult, x->rddiv, rate, distortion);

    if (this_rd < best_rd) {
      best_rd = this_rd;
      d = distortion;
      r = rate;
      *rate_tokenonly = rate_to;
      mode_selected = mode;
    }
  }
  *rate = r;
  *distortion = d;
  x->e_mbd.mode_info_context->mbmi.uv_mode = mode_selected;
}

int vp8_cost_mv_ref(VP8_COMP *cpi,
                    MB_PREDICTION_MODE m,
                    const int near_mv_ref_ct[4]) {
  MACROBLOCKD *xd = &cpi->mb.e_mbd;
  int segment_id = xd->mode_info_context->mbmi.segment_id;

  // If the mode coding is done entirely at the segment level
  // we should not account for it at the per mb level in rd code.
  // Note that if the segment level coding is expanded from single mode
  // to multiple mode masks as per reference frame coding we will need
  // to do something different here.
  if (!segfeature_active(xd, segment_id, SEG_LVL_MODE)) {
    VP8_COMMON *pc = &cpi->common;

    vp8_prob p [VP8_MVREFS - 1];
    assert(NEARESTMV <= m  &&  m <= SPLITMV);
    vp8_mv_ref_probs(pc, p, near_mv_ref_ct);
    return vp8_cost_token(vp8_mv_ref_tree, p,
                          vp8_mv_ref_encoding_array - NEARESTMV + m);
  } else
    return 0;
}

void vp8_set_mbmode_and_mvs(MACROBLOCK *x, MB_PREDICTION_MODE mb, int_mv *mv) {
  x->e_mbd.mode_info_context->mbmi.mode = mb;
  x->e_mbd.mode_info_context->mbmi.mv.as_int = mv->as_int;
}

static int labels2mode(MACROBLOCK *x, int const *labelings, int which_label,
                       B_PREDICTION_MODE this_mode,
                       int_mv *this_mv, int_mv *this_second_mv,
                       int_mv seg_mvs[MAX_REF_FRAMES - 1], int_mv *best_ref_mv,
                       int_mv *second_best_ref_mv, int *mvcost[2]) {
  MACROBLOCKD *const xd = & x->e_mbd;
  MODE_INFO *const mic = xd->mode_info_context;
  const int mis = xd->mode_info_stride;

  int i, cost = 0, thismvcost = 0;

  /* We have to be careful retrieving previously-encoded motion vectors.
     Ones from this macroblock have to be pulled from the BLOCKD array
     as they have not yet made it to the bmi array in our MB_MODE_INFO. */
  for (i = 0; i < 16; ++i) {
    BLOCKD *const d = xd->block + i;
    const int row = i >> 2,  col = i & 3;

    B_PREDICTION_MODE m;

    if (labelings[i] != which_label)
      continue;

    if (col  &&  labelings[i] == labelings[i - 1])
      m = LEFT4X4;
    else if (row  &&  labelings[i] == labelings[i - 4])
      m = ABOVE4X4;
    else {
      // the only time we should do costing for new motion vector or mode
      // is when we are on a new label  (jbb May 08, 2007)
      switch (m = this_mode) {
        case NEW4X4 :
          if (xd->mode_info_context->mbmi.second_ref_frame) {
            this_mv->as_int        = seg_mvs[xd->mode_info_context->mbmi.ref_frame        - 1].as_int;
            this_second_mv->as_int = seg_mvs[xd->mode_info_context->mbmi.second_ref_frame - 1].as_int;
          }

          thismvcost  = vp8_mv_bit_cost(this_mv, best_ref_mv, mvcost,
                                        102, xd->allow_high_precision_mv);
          if (xd->mode_info_context->mbmi.second_ref_frame) {
            thismvcost += vp8_mv_bit_cost(this_second_mv, second_best_ref_mv,
                                          mvcost, 102,
                                          xd->allow_high_precision_mv);
          }
          break;
        case LEFT4X4:
          this_mv->as_int = col ? d[-1].bmi.as_mv.first.as_int : left_block_mv(mic, i);
          if (xd->mode_info_context->mbmi.second_ref_frame)
            this_second_mv->as_int = col ? d[-1].bmi.as_mv.second.as_int : left_block_second_mv(mic, i);
          break;
        case ABOVE4X4:
          this_mv->as_int = row ? d[-4].bmi.as_mv.first.as_int : above_block_mv(mic, i, mis);
          if (xd->mode_info_context->mbmi.second_ref_frame)
            this_second_mv->as_int = row ? d[-4].bmi.as_mv.second.as_int : above_block_second_mv(mic, i, mis);
          break;
        case ZERO4X4:
          this_mv->as_int = 0;
          if (xd->mode_info_context->mbmi.second_ref_frame)
            this_second_mv->as_int = 0;
          break;
        default:
          break;
      }

      if (m == ABOVE4X4) { // replace above with left if same
        int_mv left_mv, left_second_mv;

        left_mv.as_int = col ? d[-1].bmi.as_mv.first.as_int :
                         left_block_mv(mic, i);
        if (xd->mode_info_context->mbmi.second_ref_frame)
          left_second_mv.as_int = col ? d[-1].bmi.as_mv.second.as_int :
                                  left_block_second_mv(mic, i);

        if (left_mv.as_int == this_mv->as_int &&
            (!xd->mode_info_context->mbmi.second_ref_frame ||
             left_second_mv.as_int == this_second_mv->as_int))
          m = LEFT4X4;
      }

      cost = x->inter_bmode_costs[ m];
    }

    d->bmi.as_mv.first.as_int = this_mv->as_int;
    if (xd->mode_info_context->mbmi.second_ref_frame)
      d->bmi.as_mv.second.as_int = this_second_mv->as_int;

    x->partition_info->bmi[i].mode = m;
    x->partition_info->bmi[i].mv.as_int = this_mv->as_int;
    if (xd->mode_info_context->mbmi.second_ref_frame)
      x->partition_info->bmi[i].second_mv.as_int = this_second_mv->as_int;
  }

  cost += thismvcost;
  return cost;
}

static int rdcost_mbsegment_y(MACROBLOCK *mb, const int *labels,
                              int which_label, ENTROPY_CONTEXT *ta,
                              ENTROPY_CONTEXT *tl) {
  int b, cost = 0;
  MACROBLOCKD *x = &mb->e_mbd;

  for (b = 0; b < 16; b++)
    if (labels[ b] == which_label)
      cost += cost_coeffs(mb, x->block + b, PLANE_TYPE_Y_WITH_DC,
                          ta + vp8_block2above[b],
                          tl + vp8_block2left[b], TX_4X4);

  return cost;

}

static unsigned int vp8_encode_inter_mb_segment(MACROBLOCK *x,
                                                int const *labels,
                                                int which_label,
                                                const VP8_ENCODER_RTCD *rtcd) {
  int i;
  unsigned int distortion = 0;

  for (i = 0; i < 16; i++) {
    if (labels[i] == which_label) {
      BLOCKD *bd = &x->e_mbd.block[i];
      BLOCK *be = &x->block[i];
      int thisdistortion;

      vp8_build_inter_predictors_b(bd, 16, x->e_mbd.subpixel_predict);
      if (x->e_mbd.mode_info_context->mbmi.second_ref_frame)
        vp8_build_2nd_inter_predictors_b(bd, 16, x->e_mbd.subpixel_predict_avg);
      ENCODEMB_INVOKE(&rtcd->encodemb, subb)(be, bd, 16);
      x->vp8_short_fdct4x4(be->src_diff, be->coeff, 32);

      // set to 0 no way to account for 2nd order DC so discount
      // be->coeff[0] = 0;
      x->quantize_b(be, bd);
      thisdistortion = ENCODEMB_INVOKE(&rtcd->encodemb, berr)(
                         be->coeff, bd->dqcoeff, 16) / 4;
      distortion += thisdistortion;
    }
  }
  return distortion;
}


static const unsigned int segmentation_to_sseshift[4] = {3, 3, 2, 0};


typedef struct {
  int_mv *ref_mv, *second_ref_mv;
  int_mv mvp;

  int64_t segment_rd;
  int segment_num;
  int r;
  int d;
  int segment_yrate;
  B_PREDICTION_MODE modes[16];
  int_mv mvs[16], second_mvs[16];
  unsigned char eobs[16];

  int mvthresh;
  int *mdcounts;

  int_mv sv_mvp[4];     // save 4 mvp from 8x8
  int sv_istep[2];  // save 2 initial step_param for 16x8/8x16

} BEST_SEG_INFO;

static __inline
int mv_check_bounds(MACROBLOCK *x, int_mv *mv) {
  int r = 0;
  r |= (mv->as_mv.row >> 3) < x->mv_row_min;
  r |= (mv->as_mv.row >> 3) > x->mv_row_max;
  r |= (mv->as_mv.col >> 3) < x->mv_col_min;
  r |= (mv->as_mv.col >> 3) > x->mv_col_max;
  return r;
}

static void rd_check_segment(VP8_COMP *cpi, MACROBLOCK *x,
                             BEST_SEG_INFO *bsi, unsigned int segmentation,
                             int_mv seg_mvs[16 /* n_blocks */][MAX_REF_FRAMES - 1]) {
  int i;
  int const *labels;
  int br = 0, bd = 0;
  B_PREDICTION_MODE this_mode;


  int label_count;
  int64_t this_segment_rd = 0;
  int label_mv_thresh;
  int rate = 0;
  int sbr = 0, sbd = 0;
  int segmentyrate = 0;

  vp8_variance_fn_ptr_t *v_fn_ptr;

  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta, *tl;
  ENTROPY_CONTEXT_PLANES t_above_b, t_left_b;
  ENTROPY_CONTEXT *ta_b, *tl_b;

  vpx_memcpy(&t_above, x->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, x->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;
  ta_b = (ENTROPY_CONTEXT *)&t_above_b;
  tl_b = (ENTROPY_CONTEXT *)&t_left_b;

  v_fn_ptr = &cpi->fn_ptr[segmentation];
  labels = vp8_mbsplits[segmentation];
  label_count = vp8_mbsplit_count[segmentation];

  // 64 makes this threshold really big effectively
  // making it so that we very rarely check mvs on
  // segments.   setting this to 1 would make mv thresh
  // roughly equal to what it is for macroblocks
  label_mv_thresh = 1 * bsi->mvthresh / label_count;

  // Segmentation method overheads
  rate = vp8_cost_token(vp8_mbsplit_tree, vp8_mbsplit_probs, vp8_mbsplit_encodings + segmentation);
  rate += vp8_cost_mv_ref(cpi, SPLITMV, bsi->mdcounts);
  this_segment_rd += RDCOST(x->rdmult, x->rddiv, rate, 0);
  br += rate;

  for (i = 0; i < label_count; i++) {
    int_mv mode_mv[B_MODE_COUNT], second_mode_mv[B_MODE_COUNT];
    int64_t best_label_rd = INT64_MAX;
    B_PREDICTION_MODE mode_selected = ZERO4X4;
    int bestlabelyrate = 0;

    // search for the best motion vector on this segment
    for (this_mode = LEFT4X4; this_mode <= NEW4X4; this_mode ++) {
      int64_t this_rd;
      int distortion;
      int labelyrate;
      ENTROPY_CONTEXT_PLANES t_above_s, t_left_s;
      ENTROPY_CONTEXT *ta_s;
      ENTROPY_CONTEXT *tl_s;

      vpx_memcpy(&t_above_s, &t_above, sizeof(ENTROPY_CONTEXT_PLANES));
      vpx_memcpy(&t_left_s, &t_left, sizeof(ENTROPY_CONTEXT_PLANES));

      ta_s = (ENTROPY_CONTEXT *)&t_above_s;
      tl_s = (ENTROPY_CONTEXT *)&t_left_s;

      // motion search for newmv (single predictor case only)
      if (!x->e_mbd.mode_info_context->mbmi.second_ref_frame && this_mode == NEW4X4) {
        int sseshift, n;
        int step_param = 0;
        int further_steps;
        int thissme, bestsme = INT_MAX;
        BLOCK *c;
        BLOCKD *e;

        // Is the best so far sufficiently good that we cant justify doing and new motion search.
        if (best_label_rd < label_mv_thresh)
          break;

        if (cpi->compressor_speed) {
          if (segmentation == BLOCK_8X16 || segmentation == BLOCK_16X8) {
            bsi->mvp.as_int = bsi->sv_mvp[i].as_int;
            if (i == 1 && segmentation == BLOCK_16X8)
              bsi->mvp.as_int = bsi->sv_mvp[2].as_int;

            step_param = bsi->sv_istep[i];
          }

          // use previous block's result as next block's MV predictor.
          if (segmentation == BLOCK_4X4 && i > 0) {
            bsi->mvp.as_int = x->e_mbd.block[i - 1].bmi.as_mv.first.as_int;
            if (i == 4 || i == 8 || i == 12)
              bsi->mvp.as_int = x->e_mbd.block[i - 4].bmi.as_mv.first.as_int;
            step_param = 2;
          }
        }

        further_steps = (MAX_MVSEARCH_STEPS - 1) - step_param;

        {
          int sadpb = x->sadperbit4;
          int_mv mvp_full;

          mvp_full.as_mv.row = bsi->mvp.as_mv.row >> 3;
          mvp_full.as_mv.col = bsi->mvp.as_mv.col >> 3;

          // find first label
          n = vp8_mbsplit_offset[segmentation][i];

          c = &x->block[n];
          e = &x->e_mbd.block[n];

          bestsme = vp8_full_pixel_diamond(cpi, x, c, e, &mvp_full, step_param,
                                           sadpb, further_steps, 0, v_fn_ptr,
                                           bsi->ref_mv, &mode_mv[NEW4X4]);

          sseshift = segmentation_to_sseshift[segmentation];

          // Should we do a full search (best quality only)
          if ((cpi->compressor_speed == 0) && (bestsme >> sseshift) > 4000) {
            /* Check if mvp_full is within the range. */
            vp8_clamp_mv(&mvp_full, x->mv_col_min, x->mv_col_max, x->mv_row_min, x->mv_row_max);

            thissme = cpi->full_search_sad(x, c, e, &mvp_full,
                                           sadpb, 16, v_fn_ptr,
                                           XMVCOST, bsi->ref_mv);

            if (thissme < bestsme) {
              bestsme = thissme;
              mode_mv[NEW4X4].as_int = e->bmi.as_mv.first.as_int;
            } else {
              // The full search result is actually worse so re-instate the previous best vector
              e->bmi.as_mv.first.as_int = mode_mv[NEW4X4].as_int;
            }
          }
        }

        if (bestsme < INT_MAX) {
          int distortion;
          unsigned int sse;
          cpi->find_fractional_mv_step(x, c, e, &mode_mv[NEW4X4],
                                       bsi->ref_mv, x->errorperbit, v_fn_ptr, XMVCOST,
                                       &distortion, &sse);

          // safe motion search result for use in compound prediction
          seg_mvs[i][x->e_mbd.mode_info_context->mbmi.ref_frame - 1].as_int = mode_mv[NEW4X4].as_int;
        }
      } /* NEW4X4 */
      else if (x->e_mbd.mode_info_context->mbmi.second_ref_frame && this_mode == NEW4X4) {
        // motion search not completed? Then skip newmv for this block with comppred
        if (seg_mvs[i][x->e_mbd.mode_info_context->mbmi.second_ref_frame - 1].as_int == INVALID_MV ||
            seg_mvs[i][x->e_mbd.mode_info_context->mbmi.ref_frame        - 1].as_int == INVALID_MV) {
          continue;
        }
      }

      rate = labels2mode(x, labels, i, this_mode, &mode_mv[this_mode],
                         &second_mode_mv[this_mode], seg_mvs[i], bsi->ref_mv, bsi->second_ref_mv, XMVCOST);

      // Trap vectors that reach beyond the UMV borders
      if (((mode_mv[this_mode].as_mv.row >> 3) < x->mv_row_min) || ((mode_mv[this_mode].as_mv.row >> 3) > x->mv_row_max) ||
          ((mode_mv[this_mode].as_mv.col >> 3) < x->mv_col_min) || ((mode_mv[this_mode].as_mv.col >> 3) > x->mv_col_max)) {
        continue;
      }
      if (x->e_mbd.mode_info_context->mbmi.second_ref_frame &&
          mv_check_bounds(x, &second_mode_mv[this_mode]))
        continue;

      distortion = vp8_encode_inter_mb_segment(
                     x, labels, i,
                     IF_RTCD(&cpi->rtcd));

      labelyrate = rdcost_mbsegment_y(x, labels, i, ta_s, tl_s);
      rate += labelyrate;

      this_rd = RDCOST(x->rdmult, x->rddiv, rate, distortion);

      if (this_rd < best_label_rd) {
        sbr = rate;
        sbd = distortion;
        bestlabelyrate = labelyrate;
        mode_selected = this_mode;
        best_label_rd = this_rd;

        vpx_memcpy(ta_b, ta_s, sizeof(ENTROPY_CONTEXT_PLANES));
        vpx_memcpy(tl_b, tl_s, sizeof(ENTROPY_CONTEXT_PLANES));

      }
    } /*for each 4x4 mode*/

    vpx_memcpy(ta, ta_b, sizeof(ENTROPY_CONTEXT_PLANES));
    vpx_memcpy(tl, tl_b, sizeof(ENTROPY_CONTEXT_PLANES));

    labels2mode(x, labels, i, mode_selected, &mode_mv[mode_selected],
                &second_mode_mv[mode_selected], seg_mvs[i], bsi->ref_mv, bsi->second_ref_mv, XMVCOST);

    br += sbr;
    bd += sbd;
    segmentyrate += bestlabelyrate;
    this_segment_rd += best_label_rd;

    if (this_segment_rd >= bsi->segment_rd) {
      break;
    }


  } /* for each label */

  if (this_segment_rd < bsi->segment_rd) {
    bsi->r = br;
    bsi->d = bd;
    bsi->segment_yrate = segmentyrate;
    bsi->segment_rd = this_segment_rd;
    bsi->segment_num = segmentation;

    // store everything needed to come back to this!!
    for (i = 0; i < 16; i++) {
      BLOCKD *bd = &x->e_mbd.block[i];

      bsi->mvs[i].as_mv = x->partition_info->bmi[i].mv.as_mv;
      if (x->e_mbd.mode_info_context->mbmi.second_ref_frame)
        bsi->second_mvs[i].as_mv = x->partition_info->bmi[i].second_mv.as_mv;
      bsi->modes[i] = x->partition_info->bmi[i].mode;
      bsi->eobs[i] = bd->eob;
    }
  }
}

static __inline
void vp8_cal_step_param(int sr, int *sp) {
  int step = 0;

  if (sr > MAX_FIRST_STEP) sr = MAX_FIRST_STEP;
  else if (sr < 1) sr = 1;

  while (sr >>= 1)
    step++;

  *sp = MAX_MVSEARCH_STEPS - 1 - step;
}

static int vp8_rd_pick_best_mbsegmentation(VP8_COMP *cpi, MACROBLOCK *x,
                                           int_mv *best_ref_mv, int_mv *second_best_ref_mv, int64_t best_rd,
                                           int *mdcounts, int *returntotrate,
                                           int *returnyrate, int *returndistortion,
                                           int mvthresh,
                                           int_mv seg_mvs[BLOCK_MAX_SEGMENTS - 1][16 /* n_blocks */][MAX_REF_FRAMES - 1]) {
  int i;
  BEST_SEG_INFO bsi;

  vpx_memset(&bsi, 0, sizeof(bsi));

  bsi.segment_rd = best_rd;
  bsi.ref_mv = best_ref_mv;
  bsi.second_ref_mv = second_best_ref_mv;
  bsi.mvp.as_int = best_ref_mv->as_int;
  bsi.mvthresh = mvthresh;
  bsi.mdcounts = mdcounts;

  for (i = 0; i < 16; i++)
    bsi.modes[i] = ZERO4X4;

  if (cpi->compressor_speed == 0) {
    /* for now, we will keep the original segmentation order
       when in best quality mode */
    rd_check_segment(cpi, x, &bsi, BLOCK_16X8, seg_mvs[BLOCK_16X8]);
    rd_check_segment(cpi, x, &bsi, BLOCK_8X16, seg_mvs[BLOCK_8X16]);
    rd_check_segment(cpi, x, &bsi, BLOCK_8X8,  seg_mvs[BLOCK_8X8]);
    rd_check_segment(cpi, x, &bsi, BLOCK_4X4,  seg_mvs[BLOCK_4X4]);
  } else {
    int sr;

    rd_check_segment(cpi, x, &bsi, BLOCK_8X8, seg_mvs[BLOCK_8X8]);


    if (bsi.segment_rd < best_rd) {
      int tmp_col_min = x->mv_col_min;
      int tmp_col_max = x->mv_col_max;
      int tmp_row_min = x->mv_row_min;
      int tmp_row_max = x->mv_row_max;

      vp8_clamp_mv_min_max(x, best_ref_mv);

      /* Get 8x8 result */
      bsi.sv_mvp[0].as_int = bsi.mvs[0].as_int;
      bsi.sv_mvp[1].as_int = bsi.mvs[2].as_int;
      bsi.sv_mvp[2].as_int = bsi.mvs[8].as_int;
      bsi.sv_mvp[3].as_int = bsi.mvs[10].as_int;

      /* Use 8x8 result as 16x8/8x16's predictor MV. Adjust search range according to the closeness of 2 MV. */
      /* block 8X16 */
      {
        sr = MAXF((abs(bsi.sv_mvp[0].as_mv.row - bsi.sv_mvp[2].as_mv.row)) >> 3, (abs(bsi.sv_mvp[0].as_mv.col - bsi.sv_mvp[2].as_mv.col)) >> 3);
        vp8_cal_step_param(sr, &bsi.sv_istep[0]);

        sr = MAXF((abs(bsi.sv_mvp[1].as_mv.row - bsi.sv_mvp[3].as_mv.row)) >> 3, (abs(bsi.sv_mvp[1].as_mv.col - bsi.sv_mvp[3].as_mv.col)) >> 3);
        vp8_cal_step_param(sr, &bsi.sv_istep[1]);

        rd_check_segment(cpi, x, &bsi, BLOCK_8X16, seg_mvs[BLOCK_8X16]);
      }

      /* block 16X8 */
      {
        sr = MAXF((abs(bsi.sv_mvp[0].as_mv.row - bsi.sv_mvp[1].as_mv.row)) >> 3, (abs(bsi.sv_mvp[0].as_mv.col - bsi.sv_mvp[1].as_mv.col)) >> 3);
        vp8_cal_step_param(sr, &bsi.sv_istep[0]);

        sr = MAXF((abs(bsi.sv_mvp[2].as_mv.row - bsi.sv_mvp[3].as_mv.row)) >> 3, (abs(bsi.sv_mvp[2].as_mv.col - bsi.sv_mvp[3].as_mv.col)) >> 3);
        vp8_cal_step_param(sr, &bsi.sv_istep[1]);

        rd_check_segment(cpi, x, &bsi, BLOCK_16X8, seg_mvs[BLOCK_16X8]);
      }

      /* If 8x8 is better than 16x8/8x16, then do 4x4 search */
      /* Not skip 4x4 if speed=0 (good quality) */
      if (cpi->sf.no_skip_block4x4_search || bsi.segment_num == BLOCK_8X8) { /* || (sv_segment_rd8x8-bsi.segment_rd) < sv_segment_rd8x8>>5) */
        bsi.mvp.as_int = bsi.sv_mvp[0].as_int;
        rd_check_segment(cpi, x, &bsi, BLOCK_4X4, seg_mvs[BLOCK_4X4]);
      }

      /* restore UMV window */
      x->mv_col_min = tmp_col_min;
      x->mv_col_max = tmp_col_max;
      x->mv_row_min = tmp_row_min;
      x->mv_row_max = tmp_row_max;
    }
  }

  /* set it to the best */
  for (i = 0; i < 16; i++) {
    BLOCKD *bd = &x->e_mbd.block[i];

    bd->bmi.as_mv.first.as_int = bsi.mvs[i].as_int;
    if (x->e_mbd.mode_info_context->mbmi.second_ref_frame)
      bd->bmi.as_mv.second.as_int = bsi.second_mvs[i].as_int;
    bd->eob = bsi.eobs[i];
  }

  *returntotrate = bsi.r;
  *returndistortion = bsi.d;
  *returnyrate = bsi.segment_yrate;

  /* save partitions */
  x->e_mbd.mode_info_context->mbmi.partitioning = bsi.segment_num;
  x->partition_info->count = vp8_mbsplit_count[bsi.segment_num];

  for (i = 0; i < x->partition_info->count; i++) {
    int j;

    j = vp8_mbsplit_offset[bsi.segment_num][i];

    x->partition_info->bmi[i].mode = bsi.modes[j];
    x->partition_info->bmi[i].mv.as_mv = bsi.mvs[j].as_mv;
    if (x->e_mbd.mode_info_context->mbmi.second_ref_frame)
      x->partition_info->bmi[i].second_mv.as_mv = bsi.second_mvs[j].as_mv;
  }
  /*
   * used to set x->e_mbd.mode_info_context->mbmi.mv.as_int
   */
  x->partition_info->bmi[15].mv.as_int = bsi.mvs[15].as_int;
  if (x->e_mbd.mode_info_context->mbmi.second_ref_frame)
    x->partition_info->bmi[15].second_mv.as_int = bsi.second_mvs[15].as_int;

  return bsi.segment_rd;
}

/* Order arr in increasing order, original position stored in idx */
static void insertsortmv(int arr[], int len) {
  int i, j, k;

  for (i = 1; i <= len - 1; i++) {
    for (j = 0; j < i; j++) {
      if (arr[j] > arr[i]) {
        int temp;

        temp = arr[i];

        for (k = i; k > j; k--)
          arr[k] = arr[k - 1];

        arr[j] = temp;
      }
    }
  }
}

static void insertsortsad(int arr[], int idx[], int len) {
  int i, j, k;

  for (i = 1; i <= len - 1; i++) {
    for (j = 0; j < i; j++) {
      if (arr[j] > arr[i]) {
        int temp, tempi;

        temp = arr[i];
        tempi = idx[i];

        for (k = i; k > j; k--) {
          arr[k] = arr[k - 1];
          idx[k] = idx[k - 1];
        }

        arr[j] = temp;
        idx[j] = tempi;
      }
    }
  }
}

// The improved MV prediction
void vp8_mv_pred(VP8_COMP *cpi, MACROBLOCKD *xd, const MODE_INFO *here,
                 int_mv *mvp, int refframe, int *ref_frame_sign_bias,
                 int *sr, int near_sadidx[]) {
  const MODE_INFO *above = here - xd->mode_info_stride;
  const MODE_INFO *left = here - 1;
  const MODE_INFO *aboveleft = above - 1;
  int_mv           near_mvs[8];
  int              near_ref[8];
  int_mv           mv;
  int              vcnt = 0;
  int              find = 0;
  int              mb_offset;

  int              mvx[8];
  int              mvy[8];
  int              i;

  mv.as_int = 0;

  if (here->mbmi.ref_frame != INTRA_FRAME) {
    near_mvs[0].as_int = near_mvs[1].as_int = near_mvs[2].as_int = near_mvs[3].as_int = near_mvs[4].as_int = near_mvs[5].as_int = near_mvs[6].as_int = near_mvs[7].as_int = 0;
    near_ref[0] = near_ref[1] = near_ref[2] = near_ref[3] = near_ref[4] = near_ref[5] = near_ref[6] = near_ref[7] = 0;

    // read in 3 nearby block's MVs from current frame as prediction candidates.
    if (above->mbmi.ref_frame != INTRA_FRAME) {
      near_mvs[vcnt].as_int = above->mbmi.mv.as_int;
      mv_bias(ref_frame_sign_bias[above->mbmi.ref_frame], refframe, &near_mvs[vcnt], ref_frame_sign_bias);
      near_ref[vcnt] =  above->mbmi.ref_frame;
    }
    vcnt++;
    if (left->mbmi.ref_frame != INTRA_FRAME) {
      near_mvs[vcnt].as_int = left->mbmi.mv.as_int;
      mv_bias(ref_frame_sign_bias[left->mbmi.ref_frame], refframe, &near_mvs[vcnt], ref_frame_sign_bias);
      near_ref[vcnt] =  left->mbmi.ref_frame;
    }
    vcnt++;
    if (aboveleft->mbmi.ref_frame != INTRA_FRAME) {
      near_mvs[vcnt].as_int = aboveleft->mbmi.mv.as_int;
      mv_bias(ref_frame_sign_bias[aboveleft->mbmi.ref_frame], refframe, &near_mvs[vcnt], ref_frame_sign_bias);
      near_ref[vcnt] =  aboveleft->mbmi.ref_frame;
    }
    vcnt++;

    // read in 5 nearby block's MVs from last frame.
    if (cpi->common.last_frame_type != KEY_FRAME) {
      mb_offset = (-xd->mb_to_top_edge / 128 + 1) * (xd->mode_info_stride + 1) + (-xd->mb_to_left_edge / 128 + 1);

      // current in last frame
      if (cpi->lf_ref_frame[mb_offset] != INTRA_FRAME) {
        near_mvs[vcnt].as_int = cpi->lfmv[mb_offset].as_int;
        mv_bias(cpi->lf_ref_frame_sign_bias[mb_offset], refframe, &near_mvs[vcnt], ref_frame_sign_bias);
        near_ref[vcnt] =  cpi->lf_ref_frame[mb_offset];
      }
      vcnt++;

      // above in last frame
      if (cpi->lf_ref_frame[mb_offset - xd->mode_info_stride - 1] != INTRA_FRAME) {
        near_mvs[vcnt].as_int = cpi->lfmv[mb_offset - xd->mode_info_stride - 1].as_int;
        mv_bias(cpi->lf_ref_frame_sign_bias[mb_offset - xd->mode_info_stride - 1], refframe, &near_mvs[vcnt], ref_frame_sign_bias);
        near_ref[vcnt] =  cpi->lf_ref_frame[mb_offset - xd->mode_info_stride - 1];
      }
      vcnt++;

      // left in last frame
      if (cpi->lf_ref_frame[mb_offset - 1] != INTRA_FRAME) {
        near_mvs[vcnt].as_int = cpi->lfmv[mb_offset - 1].as_int;
        mv_bias(cpi->lf_ref_frame_sign_bias[mb_offset - 1], refframe, &near_mvs[vcnt], ref_frame_sign_bias);
        near_ref[vcnt] =  cpi->lf_ref_frame[mb_offset - 1];
      }
      vcnt++;

      // right in last frame
      if (cpi->lf_ref_frame[mb_offset + 1] != INTRA_FRAME) {
        near_mvs[vcnt].as_int = cpi->lfmv[mb_offset + 1].as_int;
        mv_bias(cpi->lf_ref_frame_sign_bias[mb_offset + 1], refframe, &near_mvs[vcnt], ref_frame_sign_bias);
        near_ref[vcnt] =  cpi->lf_ref_frame[mb_offset + 1];
      }
      vcnt++;

      // below in last frame
      if (cpi->lf_ref_frame[mb_offset + xd->mode_info_stride + 1] != INTRA_FRAME) {
        near_mvs[vcnt].as_int = cpi->lfmv[mb_offset + xd->mode_info_stride + 1].as_int;
        mv_bias(cpi->lf_ref_frame_sign_bias[mb_offset + xd->mode_info_stride + 1], refframe, &near_mvs[vcnt], ref_frame_sign_bias);
        near_ref[vcnt] =  cpi->lf_ref_frame[mb_offset + xd->mode_info_stride + 1];
      }
      vcnt++;
    }

    for (i = 0; i < vcnt; i++) {
      if (near_ref[near_sadidx[i]] != INTRA_FRAME) {
        if (here->mbmi.ref_frame == near_ref[near_sadidx[i]]) {
          mv.as_int = near_mvs[near_sadidx[i]].as_int;
          find = 1;
          if (i < 3)
            *sr = 3;
          else
            *sr = 2;
          break;
        }
      }
    }

    if (!find) {
      for (i = 0; i < vcnt; i++) {
        mvx[i] = near_mvs[i].as_mv.row;
        mvy[i] = near_mvs[i].as_mv.col;
      }

      insertsortmv(mvx, vcnt);
      insertsortmv(mvy, vcnt);
      mv.as_mv.row = mvx[vcnt / 2];
      mv.as_mv.col = mvy[vcnt / 2];

      find = 1;
      // sr is set to 0 to allow calling function to decide the search range.
      *sr = 0;
    }
  }

  /* Set up return values */
  mvp->as_int = mv.as_int;
  vp8_clamp_mv2(mvp, xd);
}

void vp8_cal_sad(VP8_COMP *cpi, MACROBLOCKD *xd, MACROBLOCK *x, int recon_yoffset, int near_sadidx[]) {

  int near_sad[8] = {0}; // 0-cf above, 1-cf left, 2-cf aboveleft, 3-lf current, 4-lf above, 5-lf left, 6-lf right, 7-lf below
  BLOCK *b = &x->block[0];
  unsigned char *src_y_ptr = *(b->base_src);

  // calculate sad for current frame 3 nearby MBs.
  if (xd->mb_to_top_edge == 0 && xd->mb_to_left_edge == 0) {
    near_sad[0] = near_sad[1] = near_sad[2] = INT_MAX;
  } else if (xd->mb_to_top_edge == 0) {
    // only has left MB for sad calculation.
    near_sad[0] = near_sad[2] = INT_MAX;
    near_sad[1] = cpi->fn_ptr[BLOCK_16X16].sdf(src_y_ptr, b->src_stride, xd->dst.y_buffer - 16, xd->dst.y_stride, 0x7fffffff);
  } else if (xd->mb_to_left_edge == 0) {
    // only has left MB for sad calculation.
    near_sad[1] = near_sad[2] = INT_MAX;
    near_sad[0] = cpi->fn_ptr[BLOCK_16X16].sdf(src_y_ptr, b->src_stride, xd->dst.y_buffer - xd->dst.y_stride * 16, xd->dst.y_stride, 0x7fffffff);
  } else {
    near_sad[0] = cpi->fn_ptr[BLOCK_16X16].sdf(src_y_ptr, b->src_stride, xd->dst.y_buffer - xd->dst.y_stride * 16, xd->dst.y_stride, 0x7fffffff);
    near_sad[1] = cpi->fn_ptr[BLOCK_16X16].sdf(src_y_ptr, b->src_stride, xd->dst.y_buffer - 16, xd->dst.y_stride, 0x7fffffff);
    near_sad[2] = cpi->fn_ptr[BLOCK_16X16].sdf(src_y_ptr, b->src_stride, xd->dst.y_buffer - xd->dst.y_stride * 16 - 16, xd->dst.y_stride, 0x7fffffff);
  }

  if (cpi->common.last_frame_type != KEY_FRAME) {
    // calculate sad for last frame 5 nearby MBs.
    unsigned char *pre_y_buffer = cpi->common.yv12_fb[cpi->common.lst_fb_idx].y_buffer + recon_yoffset;
    int pre_y_stride = cpi->common.yv12_fb[cpi->common.lst_fb_idx].y_stride;

    if (xd->mb_to_top_edge == 0) near_sad[4] = INT_MAX;
    if (xd->mb_to_left_edge == 0) near_sad[5] = INT_MAX;
    if (xd->mb_to_right_edge == 0) near_sad[6] = INT_MAX;
    if (xd->mb_to_bottom_edge == 0) near_sad[7] = INT_MAX;

    if (near_sad[4] != INT_MAX)
      near_sad[4] = cpi->fn_ptr[BLOCK_16X16].sdf(src_y_ptr, b->src_stride, pre_y_buffer - pre_y_stride * 16, pre_y_stride, 0x7fffffff);
    if (near_sad[5] != INT_MAX)
      near_sad[5] = cpi->fn_ptr[BLOCK_16X16].sdf(src_y_ptr, b->src_stride, pre_y_buffer - 16, pre_y_stride, 0x7fffffff);
    near_sad[3] = cpi->fn_ptr[BLOCK_16X16].sdf(src_y_ptr, b->src_stride, pre_y_buffer, pre_y_stride, 0x7fffffff);
    if (near_sad[6] != INT_MAX)
      near_sad[6] = cpi->fn_ptr[BLOCK_16X16].sdf(src_y_ptr, b->src_stride, pre_y_buffer + 16, pre_y_stride, 0x7fffffff);
    if (near_sad[7] != INT_MAX)
      near_sad[7] = cpi->fn_ptr[BLOCK_16X16].sdf(src_y_ptr, b->src_stride, pre_y_buffer + pre_y_stride * 16, pre_y_stride, 0x7fffffff);
  }

  if (cpi->common.last_frame_type != KEY_FRAME) {
    insertsortsad(near_sad, near_sadidx, 8);
  } else {
    insertsortsad(near_sad, near_sadidx, 3);
  }
}

void rd_update_mvcount(VP8_COMP *cpi, MACROBLOCK *x,
                       int_mv *best_ref_mv, int_mv *second_best_ref_mv) {
  if (x->e_mbd.mode_info_context->mbmi.mode == SPLITMV) {
    int i;

    for (i = 0; i < x->partition_info->count; i++) {
      if (x->partition_info->bmi[i].mode == NEW4X4) {
        if (x->e_mbd.allow_high_precision_mv) {
          cpi->MVcount_hp[0][mv_max_hp + (x->partition_info->bmi[i].mv.as_mv.row
                                          - best_ref_mv->as_mv.row)]++;
          cpi->MVcount_hp[1][mv_max_hp + (x->partition_info->bmi[i].mv.as_mv.col
                                          - best_ref_mv->as_mv.col)]++;
          if (x->e_mbd.mode_info_context->mbmi.second_ref_frame) {
            cpi->MVcount_hp[0][mv_max_hp + (x->partition_info->bmi[i].second_mv.as_mv.row
                                            - second_best_ref_mv->as_mv.row)]++;
            cpi->MVcount_hp[1][mv_max_hp + (x->partition_info->bmi[i].second_mv.as_mv.col
                                            - second_best_ref_mv->as_mv.col)]++;
          }
        } else
        {
          cpi->MVcount[0][mv_max + ((x->partition_info->bmi[i].mv.as_mv.row
                                     - best_ref_mv->as_mv.row) >> 1)]++;
          cpi->MVcount[1][mv_max + ((x->partition_info->bmi[i].mv.as_mv.col
                                     - best_ref_mv->as_mv.col) >> 1)]++;
          if (x->e_mbd.mode_info_context->mbmi.second_ref_frame) {
            cpi->MVcount[0][mv_max + ((x->partition_info->bmi[i].second_mv.as_mv.row
                                       - second_best_ref_mv->as_mv.row) >> 1)]++;
            cpi->MVcount[1][mv_max + ((x->partition_info->bmi[i].second_mv.as_mv.col
                                       - second_best_ref_mv->as_mv.col) >> 1)]++;
          }
        }
      }
    }
  } else if (x->e_mbd.mode_info_context->mbmi.mode == NEWMV) {
    if (x->e_mbd.allow_high_precision_mv) {
      cpi->MVcount_hp[0][mv_max_hp + (x->e_mbd.mode_info_context->mbmi.mv.as_mv.row
                                      - best_ref_mv->as_mv.row)]++;
      cpi->MVcount_hp[1][mv_max_hp + (x->e_mbd.mode_info_context->mbmi.mv.as_mv.col
                                      - best_ref_mv->as_mv.col)]++;
      if (x->e_mbd.mode_info_context->mbmi.second_ref_frame) {
        cpi->MVcount_hp[0][mv_max_hp + (x->e_mbd.mode_info_context->mbmi.second_mv.as_mv.row
                                        - second_best_ref_mv->as_mv.row)]++;
        cpi->MVcount_hp[1][mv_max_hp + (x->e_mbd.mode_info_context->mbmi.second_mv.as_mv.col
                                        - second_best_ref_mv->as_mv.col)]++;
      }
    } else
    {
      cpi->MVcount[0][mv_max + ((x->e_mbd.mode_info_context->mbmi.mv.as_mv.row
                                 - best_ref_mv->as_mv.row) >> 1)]++;
      cpi->MVcount[1][mv_max + ((x->e_mbd.mode_info_context->mbmi.mv.as_mv.col
                                 - best_ref_mv->as_mv.col) >> 1)]++;
      if (x->e_mbd.mode_info_context->mbmi.second_ref_frame) {
        cpi->MVcount[0][mv_max + ((x->e_mbd.mode_info_context->mbmi.second_mv.as_mv.row
                                   - second_best_ref_mv->as_mv.row) >> 1)]++;
        cpi->MVcount[1][mv_max + ((x->e_mbd.mode_info_context->mbmi.second_mv.as_mv.col
                                   - second_best_ref_mv->as_mv.col) >> 1)]++;
      }
    }
  }
}

static void set_i8x8_block_modes(MACROBLOCK *x, int modes[2][4]) {
  int i;
  MACROBLOCKD *xd = &x->e_mbd;
  for (i = 0; i < 4; i++) {
    int ib = vp8_i8x8_block[i];
    x->e_mbd.mode_info_context->bmi[ib + 0].as_mode.first = modes[0][i];
    x->e_mbd.mode_info_context->bmi[ib + 1].as_mode.first = modes[0][i];
    x->e_mbd.mode_info_context->bmi[ib + 4].as_mode.first = modes[0][i];
    x->e_mbd.mode_info_context->bmi[ib + 5].as_mode.first = modes[0][i];
#if CONFIG_COMP_INTRA_PRED
    x->e_mbd.mode_info_context->bmi[ib + 0].as_mode.second = modes[1][i];
    x->e_mbd.mode_info_context->bmi[ib + 1].as_mode.second = modes[1][i];
    x->e_mbd.mode_info_context->bmi[ib + 4].as_mode.second = modes[1][i];
    x->e_mbd.mode_info_context->bmi[ib + 5].as_mode.second = modes[1][i];
#endif
    // printf("%d,%d,%d,%d %d,%d,%d,%d\n",
    //       modes[0][0], modes[0][1], modes[0][2], modes[0][3],
    //       modes[1][0], modes[1][1], modes[1][2], modes[1][3]);
  }

  for (i = 0; i < 16; i++) {
    xd->block[i].bmi = xd->mode_info_context->bmi[i];
  }
}

extern void calc_ref_probs(int *count, vp8_prob *probs);
static void estimate_curframe_refprobs(VP8_COMP *cpi, vp8_prob mod_refprobs[3], int pred_ref) {
  int norm_cnt[MAX_REF_FRAMES];
  const int *const rfct = cpi->count_mb_ref_frame_usage;
  int intra_count = rfct[INTRA_FRAME];
  int last_count  = rfct[LAST_FRAME];
  int gf_count    = rfct[GOLDEN_FRAME];
  int arf_count   = rfct[ALTREF_FRAME];

  // Work out modified reference frame probabilities to use where prediction
  // of the reference frame fails
  if (pred_ref == INTRA_FRAME) {
    norm_cnt[0] = 0;
    norm_cnt[1] = last_count;
    norm_cnt[2] = gf_count;
    norm_cnt[3] = arf_count;
    calc_ref_probs(norm_cnt, mod_refprobs);
    mod_refprobs[0] = 0;    // This branch implicit
  } else if (pred_ref == LAST_FRAME) {
    norm_cnt[0] = intra_count;
    norm_cnt[1] = 0;
    norm_cnt[2] = gf_count;
    norm_cnt[3] = arf_count;
    calc_ref_probs(norm_cnt, mod_refprobs);
    mod_refprobs[1] = 0;    // This branch implicit
  } else if (pred_ref == GOLDEN_FRAME) {
    norm_cnt[0] = intra_count;
    norm_cnt[1] = last_count;
    norm_cnt[2] = 0;
    norm_cnt[3] = arf_count;
    calc_ref_probs(norm_cnt, mod_refprobs);
    mod_refprobs[2] = 0;  // This branch implicit
  } else {
    norm_cnt[0] = intra_count;
    norm_cnt[1] = last_count;
    norm_cnt[2] = gf_count;
    norm_cnt[3] = 0;
    calc_ref_probs(norm_cnt, mod_refprobs);
    mod_refprobs[2] = 0;  // This branch implicit
  }
}

static __inline unsigned weighted_cost(vp8_prob *tab0, vp8_prob *tab1, int idx, int val, int weight) {
  unsigned cost0 = tab0[idx] ? vp8_cost_bit(tab0[idx], val) : 0;
  unsigned cost1 = tab1[idx] ? vp8_cost_bit(tab1[idx], val) : 0;
  // weight is 16-bit fixed point, so this basically calculates:
  // 0.5 + weight * cost1 + (1.0 - weight) * cost0
  return (0x8000 + weight * cost1 + (0x10000 - weight) * cost0) >> 16;
}

static void vp8_estimate_ref_frame_costs(VP8_COMP *cpi, int segment_id, unsigned int *ref_costs) {
  VP8_COMMON *cm = &cpi->common;
  MACROBLOCKD *xd = &cpi->mb.e_mbd;
  vp8_prob *mod_refprobs;

  unsigned int cost;
  int pred_ref;
  int pred_flag;
  int pred_ctx;
  int i;
  int tot_count;

  vp8_prob pred_prob, new_pred_prob;
  int seg_ref_active;
  int seg_ref_count = 0;
  seg_ref_active = segfeature_active(xd,
                                     segment_id,
                                     SEG_LVL_REF_FRAME);

  if (seg_ref_active) {
    seg_ref_count = check_segref(xd, segment_id, INTRA_FRAME)  +
                    check_segref(xd, segment_id, LAST_FRAME)   +
                    check_segref(xd, segment_id, GOLDEN_FRAME) +
                    check_segref(xd, segment_id, ALTREF_FRAME);
  }

  // Get the predicted reference for this mb
  pred_ref = get_pred_ref(cm, xd);

  // Get the context probability for the prediction flag (based on last frame)
  pred_prob = get_pred_prob(cm, xd, PRED_REF);

  // Predict probability for current frame based on stats so far
  pred_ctx = get_pred_context(cm, xd, PRED_REF);
  tot_count = cpi->ref_pred_count[pred_ctx][0] + cpi->ref_pred_count[pred_ctx][1];
  if (tot_count) {
    new_pred_prob =
      (cpi->ref_pred_count[pred_ctx][0] * 255 + (tot_count >> 1)) / tot_count;
    new_pred_prob += !new_pred_prob;
  } else
    new_pred_prob = 128;

  // Get the set of probabilities to use if prediction fails
  mod_refprobs = cm->mod_refprobs[pred_ref];

  // For each possible selected reference frame work out a cost.
  for (i = 0; i < MAX_REF_FRAMES; i++) {
    if (seg_ref_active && seg_ref_count == 1) {
      cost = 0;
    } else {
      pred_flag = (i == pred_ref);

      // Get the prediction for the current mb
      cost = weighted_cost(&pred_prob, &new_pred_prob, 0,
                           pred_flag, cpi->seg0_progress);
      if (cost > 1024) cost = 768; // i.e. account for 4 bits max.

      // for incorrectly predicted cases
      if (! pred_flag) {
        vp8_prob curframe_mod_refprobs[3];

        if (cpi->seg0_progress) {
          estimate_curframe_refprobs(cpi, curframe_mod_refprobs, pred_ref);
        } else {
          vpx_memset(curframe_mod_refprobs, 0, sizeof(curframe_mod_refprobs));
        }

        cost += weighted_cost(mod_refprobs, curframe_mod_refprobs, 0,
                              (i != INTRA_FRAME), cpi->seg0_progress);
        if (i != INTRA_FRAME) {
          cost += weighted_cost(mod_refprobs, curframe_mod_refprobs, 1,
                                (i != LAST_FRAME), cpi->seg0_progress);
          if (i != LAST_FRAME) {
            cost += weighted_cost(mod_refprobs, curframe_mod_refprobs, 2,
                                  (i != GOLDEN_FRAME), cpi->seg0_progress);
          }
        }
      }
    }

    ref_costs[i] = cost;
  }
}

static void store_coding_context(MACROBLOCK *x, int mb_index,
                                 int mode_index,
                                 PARTITION_INFO *partition,
                                 int_mv *ref_mv,
                                 int_mv *second_ref_mv) {
  MACROBLOCKD *xd = &x->e_mbd;

  // Take a snapshot of the coding context so it can be
  // restored if we decide to encode this way
  x->mb_context[mb_index].best_mode_index = mode_index;
  vpx_memcpy(&x->mb_context[mb_index].mic, xd->mode_info_context,
             sizeof(MODE_INFO));
  vpx_memcpy(&x->mb_context[mb_index].partition_info, partition,
             sizeof(PARTITION_INFO));
  x->mb_context[mb_index].best_ref_mv.as_int = ref_mv->as_int;
  x->mb_context[mb_index].second_best_ref_mv.as_int = second_ref_mv->as_int;

  // x->mb_context[mb_index].rddiv = x->rddiv;
  // x->mb_context[mb_index].rdmult = x->rdmult;
}

static void inter_mode_cost(VP8_COMP *cpi, MACROBLOCK *x, int this_mode,
                            int *rate2, int *distortion2, int *rate_y,
                            int *distortion, int* rate_uv, int *distortion_uv) {
  // Y cost and distortion
#if CONFIG_TX16X16
  if (this_mode == ZEROMV ||
      this_mode == NEARESTMV ||
      this_mode == NEARMV ||
      this_mode == NEWMV)
    macro_block_yrd_16x16(x, rate_y, distortion, IF_RTCD(&cpi->rtcd));
  else {
#endif
    if (cpi->common.txfm_mode == ALLOW_8X8)
      macro_block_yrd_8x8(x, rate_y, distortion, IF_RTCD(&cpi->rtcd));
    else
      macro_block_yrd(x, rate_y, distortion, IF_RTCD(&cpi->rtcd));
#if CONFIG_TX16X16
  }
#endif

  *rate2 += *rate_y;
  *distortion2 += *distortion;

  // UV cost and distortion
  if (cpi->common.txfm_mode == ALLOW_8X8
#if CONFIG_TX16X16
      || this_mode == ZEROMV ||
      this_mode == NEARESTMV ||
      this_mode == NEARMV ||
      this_mode == NEWMV
#endif
      )
    rd_inter16x16_uv_8x8(cpi, x, rate_uv, distortion_uv,
                         cpi->common.full_pixel);
  else
    rd_inter16x16_uv(cpi, x, rate_uv, distortion_uv, cpi->common.full_pixel);
  *rate2 += *rate_uv;
  *distortion2 += *distortion_uv;
}

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
void setup_buffer_inter(VP8_COMP *cpi, MACROBLOCK *x, int idx, int frame_type,
                        int recon_yoffset, int recon_uvoffset,
                        int_mv frame_nearest_mv[4], int_mv frame_near_mv[4],
                        int_mv frame_best_ref_mv[4],
#if CONFIG_NEWBESTREFMV
                        int_mv ref_mv[MAX_REF_FRAMES],
#endif
                        int frame_mdcounts[4][4],
                        unsigned char *y_buffer[4], unsigned char *u_buffer[4],
                        unsigned char *v_buffer[4]) {
  YV12_BUFFER_CONFIG *yv12 = &cpi->common.yv12_fb[idx];

  vp8_find_near_mvs(&x->e_mbd, x->e_mbd.mode_info_context,
                    x->e_mbd.prev_mode_info_context,
                    &frame_nearest_mv[frame_type], &frame_near_mv[frame_type],
                    &frame_best_ref_mv[frame_type], frame_mdcounts[frame_type],
                    frame_type, cpi->common.ref_frame_sign_bias);

  y_buffer[frame_type] = yv12->y_buffer + recon_yoffset;
  u_buffer[frame_type] = yv12->u_buffer + recon_uvoffset;
  v_buffer[frame_type] = yv12->v_buffer + recon_uvoffset;
#if CONFIG_NEWBESTREFMV
  vp8_find_best_ref_mvs(&x->e_mbd, y_buffer[frame_type],
                        yv12->y_stride, &frame_best_ref_mv[frame_type]);
  ref_mv[frame_type].as_int = frame_best_ref_mv[frame_type].as_int;
#endif
}

void vp8_rd_pick_inter_mode(VP8_COMP *cpi, MACROBLOCK *x, int recon_yoffset, int recon_uvoffset,
                            int *returnrate, int *returndistortion, int64_t *returnintra,
                            int64_t *best_single_rd_diff, int64_t *best_comp_rd_diff,
                            int64_t *best_hybrid_rd_diff) {
  VP8_COMMON *cm = &cpi->common;
  BLOCK *b = &x->block[0];
  BLOCKD *d = &x->e_mbd.block[0];
  MACROBLOCKD *xd = &x->e_mbd;
  union b_mode_info best_bmodes[16];
  MB_MODE_INFO best_mbmode;
  PARTITION_INFO best_partition;
  int_mv best_ref_mv, second_best_ref_mv;
  int_mv mode_mv[MB_MODE_COUNT];
  MB_PREDICTION_MODE this_mode;
  int i, best_mode_index = 0;
  int mode8x8[2][4];
  unsigned char segment_id = xd->mode_info_context->mbmi.segment_id;

  int mode_index;
  int mdcounts[4];
  int rate, distortion;
  int rate2, distortion2;
  int64_t best_rd = INT64_MAX, best_intra_rd = INT64_MAX;
  int64_t best_comp_rd = INT64_MAX;
  int64_t best_single_rd = INT64_MAX;
  int64_t best_hybrid_rd = INT64_MAX;
#if CONFIG_PRED_FILTER
  int64_t best_overall_rd = INT64_MAX;
#endif
  int uv_intra_rate, uv_intra_distortion, uv_intra_rate_tokenonly;
  int uv_intra_skippable = 0;
  int uv_intra_rate_8x8 = 0, uv_intra_distortion_8x8 = 0, uv_intra_rate_tokenonly_8x8 = 0;
  int uv_intra_skippable_8x8 = 0;
  int rate_y, UNINITIALIZED_IS_SAFE(rate_uv);
  int distortion_uv;
  int64_t best_yrd = INT64_MAX;
#if CONFIG_PRED_FILTER
  int best_filter_state;
#endif
#if CONFIG_NEWBESTREFMV
  int_mv ref_mv[MAX_REF_FRAMES] = {0};
#endif

#if CONFIG_SWITCHABLE_INTERP
  int switchable_filter_index = 0;
#endif

  MB_PREDICTION_MODE uv_intra_mode;
  MB_PREDICTION_MODE uv_intra_mode_8x8 = 0;

  int_mv mvp;
  int near_sadidx[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  int saddone = 0;
  int sr = 0;  // search range got from mv_pred(). It uses step_param levels. (0-7)

  int_mv frame_nearest_mv[4];
  int_mv frame_near_mv[4];
  int_mv frame_best_ref_mv[4];
  int_mv mc_search_result[4];
  int frame_mdcounts[4][4];
  unsigned char *y_buffer[4], *u_buffer[4], *v_buffer[4];

  unsigned int ref_costs[MAX_REF_FRAMES];
  int_mv seg_mvs[BLOCK_MAX_SEGMENTS - 1][16 /* n_blocks */][MAX_REF_FRAMES - 1];

  vpx_memset(&best_mbmode, 0, sizeof(best_mbmode));
  vpx_memset(&best_bmodes, 0, sizeof(best_bmodes));
  vpx_memset(&x->mb_context[xd->mb_index], 0, sizeof(PICK_MODE_CONTEXT));

  for (i = 0; i < 4; i++)
    mc_search_result[i].as_int = INVALID_MV;

  for (i = 0; i < BLOCK_MAX_SEGMENTS - 1; i++) {
    int j, k;

    for (j = 0; j < 16; j++)
      for (k = 0; k < MAX_REF_FRAMES - 1; k++)
        seg_mvs[i][j][k].as_int = INVALID_MV;
  }

  if (cpi->ref_frame_flags & VP8_LAST_FLAG) {
    setup_buffer_inter(cpi, x, cpi->common.lst_fb_idx, LAST_FRAME,
                       recon_yoffset, recon_uvoffset, frame_nearest_mv,
                       frame_near_mv, frame_best_ref_mv,
#if CONFIG_NEWBESTREFMV
                       ref_mv,
#endif
                       frame_mdcounts, y_buffer, u_buffer, v_buffer);
  }

  if (cpi->ref_frame_flags & VP8_GOLD_FLAG) {
    setup_buffer_inter(cpi, x, cpi->common.gld_fb_idx, GOLDEN_FRAME,
                       recon_yoffset, recon_uvoffset, frame_nearest_mv,
                       frame_near_mv, frame_best_ref_mv,
#if CONFIG_NEWBESTREFMV
                       ref_mv,
#endif
                       frame_mdcounts, y_buffer, u_buffer, v_buffer);
  }

  if (cpi->ref_frame_flags & VP8_ALT_FLAG) {
    setup_buffer_inter(cpi, x, cpi->common.alt_fb_idx, ALTREF_FRAME,
                       recon_yoffset, recon_uvoffset, frame_nearest_mv,
                       frame_near_mv, frame_best_ref_mv,
#if CONFIG_NEWBESTREFMV
                       ref_mv,
#endif
                       frame_mdcounts, y_buffer, u_buffer, v_buffer);
  }

  *returnintra = INT64_MAX;

  x->skip = 0;

  vpx_memset(mode_mv, 0, sizeof(mode_mv));

  x->e_mbd.mode_info_context->mbmi.ref_frame = INTRA_FRAME;

  /* Initialize zbin mode boost for uv costing */
  cpi->zbin_mode_boost = 0;
  vp8_update_zbin_extra(cpi, x);

  rd_pick_intra_mbuv_mode(cpi, x, &uv_intra_rate,
                          &uv_intra_rate_tokenonly, &uv_intra_distortion);
  uv_intra_mode = x->e_mbd.mode_info_context->mbmi.uv_mode;
  uv_intra_skippable = mbuv_is_skippable(&x->e_mbd);

  /* rough estimate for now */
  if (cpi->common.txfm_mode == ALLOW_8X8) {
    rd_pick_intra_mbuv_mode_8x8(cpi, x, &uv_intra_rate_8x8,
                                &uv_intra_rate_tokenonly_8x8,
                                &uv_intra_distortion_8x8);
    uv_intra_mode_8x8 = x->e_mbd.mode_info_context->mbmi.uv_mode;
    uv_intra_skippable_8x8 = mbuv_is_skippable_8x8(&x->e_mbd);
  }

  // Get estimates of reference frame costs for each reference frame
  // that depend on the current prediction etc.
  vp8_estimate_ref_frame_costs(cpi, segment_id, ref_costs);

#if CONFIG_SWITCHABLE_INTERP
  for (mode_index = 0; mode_index < MAX_MODES;
       mode_index += (!switchable_filter_index)) {
#else
  for (mode_index = 0; mode_index < MAX_MODES; ++mode_index) {
#endif
    int64_t this_rd = INT64_MAX;
    int disable_skip = 0;
    int other_cost = 0;
    int compmode_cost = 0;
    int mode_excluded = 0;

    // These variables hold are rolling total cost and distortion for this mode
    rate2 = 0;
    distortion2 = 0;
    rate_y = 0;
    rate_uv = 0;

    this_mode = vp8_mode_order[mode_index].mode;
    xd->mode_info_context->mbmi.mode = this_mode;
    xd->mode_info_context->mbmi.uv_mode = DC_PRED;
    xd->mode_info_context->mbmi.ref_frame =
      vp8_mode_order[mode_index].ref_frame;
    xd->mode_info_context->mbmi.second_ref_frame =
      vp8_mode_order[mode_index].second_ref_frame;
#if CONFIG_NEWBESTREFMV
    x->e_mbd.mode_info_context->mbmi.ref_mv =
      ref_mv[x->e_mbd.mode_info_context->mbmi.ref_frame];
    x->e_mbd.mode_info_context->mbmi.second_ref_mv =
      ref_mv[x->e_mbd.mode_info_context->mbmi.second_ref_frame];
#endif
#if CONFIG_PRED_FILTER
    xd->mode_info_context->mbmi.pred_filter_enabled = 0;
#endif
#if CONFIG_SWITCHABLE_INTERP
    if (cpi->common.mcomp_filter_type == SWITCHABLE &&
        this_mode >= NEARESTMV && this_mode <= SPLITMV) {
      xd->mode_info_context->mbmi.interp_filter =
          vp8_switchable_interp[switchable_filter_index++];
      if (switchable_filter_index == VP8_SWITCHABLE_FILTERS)
        switchable_filter_index = 0;
        //printf("Searching %d (%d)\n", this_mode, switchable_filter_index);
    } else {
      xd->mode_info_context->mbmi.interp_filter = cpi->common.mcomp_filter_type;
    }
    vp8_setup_interp_filters(xd, xd->mode_info_context->mbmi.interp_filter,
                             &cpi->common);
#endif

    // Test best rd so far against threshold for trying this mode.
    if (best_rd <= cpi->rd_threshes[mode_index])
      continue;

    // current coding mode under rate-distortion optimization test loop
#if CONFIG_HYBRIDTRANSFORM
    xd->mode_info_context->mbmi.mode_rdopt = this_mode;
#endif

#if CONFIG_COMP_INTRA_PRED
    xd->mode_info_context->mbmi.second_mode = (MB_PREDICTION_MODE)(DC_PRED - 1);
    xd->mode_info_context->mbmi.second_uv_mode = (MB_PREDICTION_MODE)(DC_PRED - 1);
#endif

    // If the segment reference frame feature is enabled....
    // then do nothing if the current ref frame is not allowed..
    if (segfeature_active(xd, segment_id, SEG_LVL_REF_FRAME) &&
        !check_segref(xd, segment_id,
                      xd->mode_info_context->mbmi.ref_frame)) {
      continue;
    }
    // If the segment mode feature is enabled....
    // then do nothing if the current mode is not allowed..
    else if (segfeature_active(xd, segment_id, SEG_LVL_MODE)  &&
             (this_mode !=
              get_segdata(xd, segment_id, SEG_LVL_MODE))) {
      continue;
    }
    // Disable this drop out case if either the mode or ref frame
    // segment level feature is enabled for this segment. This is to
    // prevent the possibility that the we end up unable to pick any mode.
    else if (!segfeature_active(xd, segment_id, SEG_LVL_REF_FRAME) &&
             !segfeature_active(xd, segment_id, SEG_LVL_MODE)) {
      // Only consider ZEROMV/ALTREF_FRAME for alt ref frame,
      // unless ARNR filtering is enabled in which case we want
      // an unfiltered alternative
      if (cpi->is_src_frame_alt_ref && (cpi->oxcf.arnr_max_frames == 0)) {
        if (this_mode != ZEROMV ||
            x->e_mbd.mode_info_context->mbmi.ref_frame != ALTREF_FRAME) {
          continue;
        }
      }
    }

    /* everything but intra */
    if (x->e_mbd.mode_info_context->mbmi.ref_frame) {
      int ref = x->e_mbd.mode_info_context->mbmi.ref_frame;

      x->e_mbd.pre.y_buffer = y_buffer[ref];
      x->e_mbd.pre.u_buffer = u_buffer[ref];
      x->e_mbd.pre.v_buffer = v_buffer[ref];
      mode_mv[NEARESTMV] = frame_nearest_mv[ref];
      mode_mv[NEARMV] = frame_near_mv[ref];
      best_ref_mv = frame_best_ref_mv[ref];
      vpx_memcpy(mdcounts, frame_mdcounts[ref], sizeof(mdcounts));
    }

    if (x->e_mbd.mode_info_context->mbmi.second_ref_frame) {
      int ref = x->e_mbd.mode_info_context->mbmi.second_ref_frame;

      x->e_mbd.second_pre.y_buffer = y_buffer[ref];
      x->e_mbd.second_pre.u_buffer = u_buffer[ref];
      x->e_mbd.second_pre.v_buffer = v_buffer[ref];
      second_best_ref_mv  = frame_best_ref_mv[ref];
    }

    // Experimental code. Special case for gf and arf zeromv modes.
    // Increase zbin size to suppress noise
    if (cpi->zbin_mode_boost_enabled) {
      if (vp8_mode_order[mode_index].ref_frame == INTRA_FRAME)
        cpi->zbin_mode_boost = 0;
      else {
        if (vp8_mode_order[mode_index].mode == ZEROMV) {
          if (vp8_mode_order[mode_index].ref_frame != LAST_FRAME)
            cpi->zbin_mode_boost = GF_ZEROMV_ZBIN_BOOST;
          else
            cpi->zbin_mode_boost = LF_ZEROMV_ZBIN_BOOST;
        } else if (vp8_mode_order[mode_index].mode == SPLITMV)
          cpi->zbin_mode_boost = 0;
        else
          cpi->zbin_mode_boost = MV_ZBIN_BOOST;
      }

      vp8_update_zbin_extra(cpi, x);
    }

    // Intra
    if (!x->e_mbd.mode_info_context->mbmi.ref_frame) {
      switch (this_mode) {
        default:
        case DC_PRED:
        case V_PRED:
        case H_PRED:
        case TM_PRED:
        case D45_PRED:
        case D135_PRED:
        case D117_PRED:
        case D153_PRED:
        case D27_PRED:
        case D63_PRED:
          x->e_mbd.mode_info_context->mbmi.ref_frame = INTRA_FRAME;
          // FIXME compound intra prediction
          RECON_INVOKE(&cpi->common.rtcd.recon, build_intra_predictors_mby)
              (&x->e_mbd);
#if CONFIG_TX16X16
          // FIXME: breaks lossless since 4x4 isn't allowed
          macro_block_yrd_16x16(x, &rate_y, &distortion,
                                IF_RTCD(&cpi->rtcd));
          rate2 += rate_y;
          distortion2 += distortion;
          rate2 += x->mbmode_cost[x->e_mbd.frame_type][x->e_mbd.mode_info_context->mbmi.mode];
          rate2 += uv_intra_rate_8x8;
          rate_uv = uv_intra_rate_tokenonly_8x8;
          distortion2 += uv_intra_distortion_8x8;
          distortion_uv = uv_intra_distortion_8x8;
          break;
#else
          if (cpi->common.txfm_mode == ALLOW_8X8)
            macro_block_yrd_8x8(x, &rate_y, &distortion,
                                IF_RTCD(&cpi->rtcd));
          else
            macro_block_yrd(x, &rate_y, &distortion,
                            IF_RTCD(&cpi->rtcd));
          rate2 += rate_y;
          distortion2 += distortion;
          rate2 += x->mbmode_cost[x->e_mbd.frame_type][x->e_mbd.mode_info_context->mbmi.mode];
          if (cpi->common.txfm_mode == ALLOW_8X8) {
            rate2 += uv_intra_rate_8x8;
            rate_uv = uv_intra_rate_tokenonly_8x8;
            distortion2 += uv_intra_distortion_8x8;
            distortion_uv = uv_intra_distortion_8x8;
          } else {
            rate2 += uv_intra_rate;
            rate_uv = uv_intra_rate_tokenonly;
            distortion2 += uv_intra_distortion;
            distortion_uv = uv_intra_distortion;
          }
          break;
#endif
        case B_PRED: {
          int64_t tmp_rd;

          // Note the rate value returned here includes the cost of coding
          // the BPRED mode : x->mbmode_cost[x->e_mbd.frame_type][BPRED];
          tmp_rd = rd_pick_intra4x4mby_modes(cpi, x, &rate, &rate_y, &distortion, best_yrd,
#if CONFIG_COMP_INTRA_PRED
                                             0,
#endif
                                             0);
          rate2 += rate;
          distortion2 += distortion;

          if (tmp_rd < best_yrd) {
            rate2 += uv_intra_rate;
            rate_uv = uv_intra_rate_tokenonly;
            distortion2 += uv_intra_distortion;
            distortion_uv = uv_intra_distortion;
          } else {
            this_rd = INT64_MAX;
            disable_skip = 1;
          }
        }
        break;
        case I8X8_PRED: {
          int64_t tmp_rd;
          tmp_rd = rd_pick_intra8x8mby_modes(cpi,
                                             x, &rate, &rate_y, &distortion,
                                             best_yrd);
          rate2 += rate;
          distortion2 += distortion;

          mode8x8[0][0] = x->e_mbd.mode_info_context->bmi[0].as_mode.first;
          mode8x8[0][1] = x->e_mbd.mode_info_context->bmi[2].as_mode.first;
          mode8x8[0][2] = x->e_mbd.mode_info_context->bmi[8].as_mode.first;
          mode8x8[0][3] = x->e_mbd.mode_info_context->bmi[10].as_mode.first;
#if CONFIG_COMP_INTRA_PRED
          mode8x8[1][0] = x->e_mbd.mode_info_context->bmi[0].as_mode.second;
          mode8x8[1][1] = x->e_mbd.mode_info_context->bmi[2].as_mode.second;
          mode8x8[1][2] = x->e_mbd.mode_info_context->bmi[8].as_mode.second;
          mode8x8[1][3] = x->e_mbd.mode_info_context->bmi[10].as_mode.second;
#endif

          /* TODO: uv rate maybe over-estimated here since there is UV intra
                   mode coded in I8X8_PRED prediction */
          if (tmp_rd < best_yrd) {
            rate2 += uv_intra_rate;
            rate_uv = uv_intra_rate_tokenonly;
            distortion2 += uv_intra_distortion;
            distortion_uv = uv_intra_distortion;
          } else {
            this_rd = INT64_MAX;
            disable_skip = 1;
          }
        }
        break;
      }
    }
    // Split MV. The code is very different from the other inter modes so
    // special case it.
    else if (this_mode == SPLITMV) {
      int64_t tmp_rd, this_rd_thresh;
      int is_comp_pred = x->e_mbd.mode_info_context->mbmi.second_ref_frame != 0;
      int_mv *second_ref = is_comp_pred ? &second_best_ref_mv : NULL;

      this_rd_thresh =
          (x->e_mbd.mode_info_context->mbmi.ref_frame == LAST_FRAME) ?
          cpi->rd_threshes[THR_NEWMV] : cpi->rd_threshes[THR_NEWA];
      this_rd_thresh =
          (x->e_mbd.mode_info_context->mbmi.ref_frame == GOLDEN_FRAME) ?
          cpi->rd_threshes[THR_NEWG] : this_rd_thresh;

      tmp_rd = vp8_rd_pick_best_mbsegmentation(cpi, x, &best_ref_mv,
                                               second_ref, best_yrd, mdcounts,
                                               &rate, &rate_y, &distortion,
                                               this_rd_thresh, seg_mvs);
      rate2 += rate;
      distortion2 += distortion;

#if CONFIG_SWITCHABLE_INTERP
      if (cpi->common.mcomp_filter_type == SWITCHABLE)
        rate2 += SWITCHABLE_INTERP_RATE_FACTOR * x->switchable_interp_costs
            [get_pred_context(&cpi->common, xd, PRED_SWITCHABLE_INTERP)]
            [vp8_switchable_interp_map[
            x->e_mbd.mode_info_context->mbmi.interp_filter]];
#endif
      // If even the 'Y' rd value of split is higher than best so far
      // then dont bother looking at UV
      if (tmp_rd < best_yrd) {
        rd_inter4x4_uv(cpi, x, &rate_uv, &distortion_uv, cpi->common.full_pixel);
        rate2 += rate_uv;
        distortion2 += distortion_uv;
      } else {
        this_rd = INT64_MAX;
        disable_skip = 1;
      }

      if (is_comp_pred)
        mode_excluded = cpi->common.comp_pred_mode == SINGLE_PREDICTION_ONLY;
      else
        mode_excluded = cpi->common.comp_pred_mode == COMP_PREDICTION_ONLY;

      compmode_cost =
        vp8_cost_bit(get_pred_prob(cm, xd, PRED_COMP), is_comp_pred);
      x->e_mbd.mode_info_context->mbmi.mode = this_mode;
    }
    // Single prediction inter
    else if (!x->e_mbd.mode_info_context->mbmi.second_ref_frame) {
      switch (this_mode) {
        case NEWMV: {
          int bestsme = INT_MAX;
          int further_steps, step_param = cpi->sf.first_step;
          int sadpb = x->sadperbit16;
          int_mv mvp_full;

          int tmp_col_min = x->mv_col_min;
          int tmp_col_max = x->mv_col_max;
          int tmp_row_min = x->mv_row_min;
          int tmp_row_max = x->mv_row_max;

          vp8_clamp_mv_min_max(x, &best_ref_mv);

          if (!saddone) {
            vp8_cal_sad(cpi, xd, x, recon_yoffset, &near_sadidx[0]);
            saddone = 1;
          }

          vp8_mv_pred(cpi, &x->e_mbd, x->e_mbd.mode_info_context, &mvp,
                      x->e_mbd.mode_info_context->mbmi.ref_frame,
                      cpi->common.ref_frame_sign_bias, &sr, &near_sadidx[0]);

          mvp_full.as_mv.col = mvp.as_mv.col >> 3;
          mvp_full.as_mv.row = mvp.as_mv.row >> 3;

          // adjust search range according to sr from mv prediction
          step_param = MAX(step_param, sr);

          // Further step/diamond searches as necessary
          further_steps = (cpi->sf.max_step_search_steps - 1) - step_param;

          bestsme = vp8_full_pixel_diamond(cpi, x, b, d, &mvp_full, step_param,
                                           sadpb, further_steps, 1,
                                           &cpi->fn_ptr[BLOCK_16X16],
                                           &best_ref_mv, &mode_mv[NEWMV]);
          d->bmi.as_mv.first.as_int = mode_mv[NEWMV].as_int;

          x->mv_col_min = tmp_col_min;
          x->mv_col_max = tmp_col_max;
          x->mv_row_min = tmp_row_min;
          x->mv_row_max = tmp_row_max;

          if (bestsme < INT_MAX) {
            int dis; /* TODO: use dis in distortion calculation later. */
            unsigned int sse;
            cpi->find_fractional_mv_step(x, b, d, &d->bmi.as_mv.first, &best_ref_mv,
                                         x->errorperbit,
                                         &cpi->fn_ptr[BLOCK_16X16],
                                         XMVCOST, &dis, &sse);
          }
          mc_search_result[x->e_mbd.mode_info_context->mbmi.ref_frame].as_int = d->bmi.as_mv.first.as_int;

          mode_mv[NEWMV].as_int = d->bmi.as_mv.first.as_int;

          // Add the new motion vector cost to our rolling cost variable
          rate2 += vp8_mv_bit_cost(&mode_mv[NEWMV], &best_ref_mv,
                                   XMVCOST, 96,
                                   x->e_mbd.allow_high_precision_mv);
        }

        case NEARESTMV:
        case NEARMV:
          // Clip "next_nearest" so that it does not extend to far out of image
          vp8_clamp_mv2(&mode_mv[this_mode], xd);

          // Do not bother proceeding if the vector (from newmv,nearest or near) is 0,0 as this should then be coded using the zeromv mode.
          if (((this_mode == NEARMV) || (this_mode == NEARESTMV)) &&
              (mode_mv[this_mode].as_int == 0)) {
            continue;
          }

        case ZEROMV:
          // Trap vectors that reach beyond the UMV borders
          // Note that ALL New MV, Nearest MV Near MV and Zero MV code drops through to this point
          // because of the lack of break statements in the previous two cases.
          if (mv_check_bounds(x, &mode_mv[this_mode]))
            continue;

          vp8_set_mbmode_and_mvs(x, this_mode, &mode_mv[this_mode]);

#if CONFIG_PRED_FILTER
          // Filtered prediction:
          xd->mode_info_context->mbmi.pred_filter_enabled =
            vp8_mode_order[mode_index].pred_filter_flag;
          rate2 += vp8_cost_bit(cpi->common.prob_pred_filter_off,
                                xd->mode_info_context->mbmi.pred_filter_enabled);
#endif
#if CONFIG_SWITCHABLE_INTERP
          if (cpi->common.mcomp_filter_type == SWITCHABLE)
            rate2 += SWITCHABLE_INTERP_RATE_FACTOR * x->switchable_interp_costs
                [get_pred_context(&cpi->common, xd, PRED_SWITCHABLE_INTERP)]
                [vp8_switchable_interp_map[
                x->e_mbd.mode_info_context->mbmi.interp_filter]];
#endif

          vp8_build_1st_inter16x16_predictors_mby(&x->e_mbd,
                                                  xd->predictor, 16);

          compmode_cost =
            vp8_cost_bit(get_pred_prob(cm, xd, PRED_COMP), 0);

          // Add in the Mv/mode cost
          rate2 += vp8_cost_mv_ref(cpi, this_mode, mdcounts);

          if (cpi->active_map_enabled && x->active_ptr[0] == 0)
            x->skip = 1;
          else if (x->encode_breakout) {
            unsigned int sse, var;
            int threshold = (xd->block[0].dequant[1]
                             * xd->block[0].dequant[1] >> 4);

            if (threshold < x->encode_breakout)
              threshold = x->encode_breakout;

            var = VARIANCE_INVOKE(&cpi->rtcd.variance, var16x16)
                  (*(b->base_src), b->src_stride,
                   x->e_mbd.predictor, 16, &sse);

            if (sse < threshold) {
              unsigned int q2dc = xd->block[24].dequant[0];
              /* If there is no codeable 2nd order dc
                 or a very small uniform pixel change change */
              if ((sse - var < q2dc *q2dc >> 4) ||
                  (sse / 2 > var && sse - var < 64)) {
                // Check u and v to make sure skip is ok
                int sse2 =  VP8_UVSSE(x, IF_RTCD(&cpi->rtcd.variance));
                if (sse2 * 2 < threshold) {
                  x->skip = 1;
                  distortion2 = sse + sse2;
                  rate2 = 500;

                  /* for best_yrd calculation */
                  rate_uv = 0;
                  distortion_uv = sse2;

                  disable_skip = 1;
                  this_rd = RDCOST(x->rdmult, x->rddiv, rate2, distortion2);

                  break;
                }
              }
            }
          }

          vp8_build_1st_inter16x16_predictors_mbuv(&x->e_mbd,
                                                   &xd->predictor[256],
                                                   &xd->predictor[320], 8);
          inter_mode_cost(cpi, x, this_mode, &rate2, &distortion2,
                          &rate_y, &distortion, &rate_uv, &distortion_uv);
          mode_excluded = cpi->common.comp_pred_mode == COMP_PREDICTION_ONLY;
          break;

        default:
          break;
      }
    } else { /* x->e_mbd.mode_info_context->mbmi.second_ref_frame != 0 */
      int ref1 = x->e_mbd.mode_info_context->mbmi.ref_frame;
      int ref2 = x->e_mbd.mode_info_context->mbmi.second_ref_frame;

      mode_excluded = cpi->common.comp_pred_mode == SINGLE_PREDICTION_ONLY;
      switch (this_mode) {
        case NEWMV:
          if (mc_search_result[ref1].as_int == INVALID_MV ||
              mc_search_result[ref2].as_int == INVALID_MV)
            continue;
          x->e_mbd.mode_info_context->mbmi.mv.as_int        = mc_search_result[ref1].as_int;
          x->e_mbd.mode_info_context->mbmi.second_mv.as_int = mc_search_result[ref2].as_int;
          rate2 += vp8_mv_bit_cost(&mc_search_result[ref1],
                                   &frame_best_ref_mv[ref1],
                                   XMVCOST, 96,
                                   x->e_mbd.allow_high_precision_mv);
          rate2 += vp8_mv_bit_cost(&mc_search_result[ref2],
                                   &frame_best_ref_mv[ref2],
                                   XMVCOST, 96,
                                   x->e_mbd.allow_high_precision_mv);
          break;
        case ZEROMV:
          x->e_mbd.mode_info_context->mbmi.mv.as_int        = 0;
          x->e_mbd.mode_info_context->mbmi.second_mv.as_int = 0;
          break;
        case NEARMV:
          if (frame_near_mv[ref1].as_int == 0 || frame_near_mv[ref2].as_int == 0)
            continue;
          x->e_mbd.mode_info_context->mbmi.mv.as_int        = frame_near_mv[ref1].as_int;
          x->e_mbd.mode_info_context->mbmi.second_mv.as_int = frame_near_mv[ref2].as_int;
          break;
        case NEARESTMV:
          if (frame_nearest_mv[ref1].as_int == 0 || frame_nearest_mv[ref2].as_int == 0)
            continue;
          x->e_mbd.mode_info_context->mbmi.mv.as_int        = frame_nearest_mv[ref1].as_int;
          x->e_mbd.mode_info_context->mbmi.second_mv.as_int = frame_nearest_mv[ref2].as_int;
          break;
        default:
          break;
      }

      /* We don't include the cost of the second reference here, because there are only
       * three options: Last/Golden, ARF/Last or Golden/ARF, or in other words if you
       * present them in that order, the second one is always known if the first is known */
      compmode_cost =
        vp8_cost_bit(get_pred_prob(cm, xd, PRED_COMP), 1);

      /* Add in the Mv/mode cost */
      rate2 += vp8_cost_mv_ref(cpi, this_mode, mdcounts);

      vp8_clamp_mv2(&x->e_mbd.mode_info_context->mbmi.mv, xd);
      vp8_clamp_mv2(&x->e_mbd.mode_info_context->mbmi.second_mv, xd);
      if (mv_check_bounds(x, &x->e_mbd.mode_info_context->mbmi.mv))
        continue;
      if (mv_check_bounds(x, &x->e_mbd.mode_info_context->mbmi.second_mv))
        continue;

      /* build first and second prediction */
      vp8_build_1st_inter16x16_predictors_mby(&x->e_mbd, x->e_mbd.predictor,
                                              16);
      vp8_build_1st_inter16x16_predictors_mbuv(&x->e_mbd, &xd->predictor[256],
                                               &xd->predictor[320], 8);
      /* do second round and average the results */
      vp8_build_2nd_inter16x16_predictors_mb(&x->e_mbd, x->e_mbd.predictor,
                                             &x->e_mbd.predictor[256],
                                             &x->e_mbd.predictor[320], 16, 8);

      inter_mode_cost(cpi, x, this_mode, &rate2, &distortion2,
                      &rate_y, &distortion, &rate_uv, &distortion_uv);

      /* don't bother w/ skip, we would never have come here if skip were enabled */
      x->e_mbd.mode_info_context->mbmi.mode = this_mode;
    }

    if (cpi->common.comp_pred_mode == HYBRID_PREDICTION)
      rate2 += compmode_cost;

    // Estimate the reference frame signaling cost and add it
    // to the rolling cost variable.
    rate2 += ref_costs[x->e_mbd.mode_info_context->mbmi.ref_frame];

    if (!disable_skip) {
      // Test for the condition where skip block will be activated
      // because there are no non zero coefficients and make any
      // necessary adjustment for rate. Ignore if skip is coded at
      // segment level as the cost wont have been added in.
      if (cpi->common.mb_no_coeff_skip) {
        int mb_skippable;
        int mb_skip_allowed;
        int has_y2 = (this_mode != SPLITMV
                      && this_mode != B_PRED
                      && this_mode != I8X8_PRED);

#if CONFIGURE_TX16X16
        if (this_mode <= TM_PRED ||
            this_mode == NEWMV ||
            this_mode == ZEROMV ||
            this_mode == NEARESTMV ||
            this_mode == NEARMV)
          mb_skippable = mb_is_skippable_16x16(&x->e_mbd);
        else
#endif
        if ((cpi->common.txfm_mode == ALLOW_8X8) && has_y2) {
          if (x->e_mbd.mode_info_context->mbmi.ref_frame != INTRA_FRAME)
            mb_skippable = mb_is_skippable_8x8(&x->e_mbd);
          else
            mb_skippable = uv_intra_skippable_8x8
                           & mby_is_skippable_8x8(&x->e_mbd);
        } else {
          if (x->e_mbd.mode_info_context->mbmi.ref_frame != INTRA_FRAME)
            mb_skippable = mb_is_skippable(&x->e_mbd, has_y2);
          else
            mb_skippable = uv_intra_skippable
                           & mby_is_skippable(&x->e_mbd, has_y2);
        }

        // Is Mb level skip allowed for this mb.
        mb_skip_allowed =
          !segfeature_active(xd, segment_id, SEG_LVL_EOB) ||
          get_segdata(xd, segment_id, SEG_LVL_EOB);

        if (mb_skippable) {
          // Back out the coefficient coding costs
          rate2 -= (rate_y + rate_uv);
          // for best_yrd calculation
          rate_uv = 0;

          if (mb_skip_allowed) {
            int prob_skip_cost;

            // Cost the skip mb case
            vp8_prob skip_prob =
              get_pred_prob(cm, &x->e_mbd, PRED_MBSKIP);

            if (skip_prob) {
              prob_skip_cost = vp8_cost_bit(skip_prob, 1);
              rate2 += prob_skip_cost;
              other_cost += prob_skip_cost;
            }
          }
        }
        // Add in the cost of the no skip flag.
        else if (mb_skip_allowed) {
          int prob_skip_cost = vp8_cost_bit(
                                 get_pred_prob(cm, &x->e_mbd, PRED_MBSKIP), 0);
          rate2 += prob_skip_cost;
          other_cost += prob_skip_cost;
        }
      }

      // Calculate the final RD estimate for this mode.
      this_rd = RDCOST(x->rdmult, x->rddiv, rate2, distortion2);
    }

    // Keep record of best intra distortion
    if ((x->e_mbd.mode_info_context->mbmi.ref_frame == INTRA_FRAME) &&
        (this_rd < best_intra_rd)) {
      best_intra_rd = this_rd;
      *returnintra = distortion2;
    }

    if (!disable_skip && x->e_mbd.mode_info_context->mbmi.ref_frame == INTRA_FRAME) {
      best_comp_rd = MIN(best_comp_rd, this_rd);
      best_single_rd = MIN(best_single_rd, this_rd);
      best_hybrid_rd = MIN(best_hybrid_rd, this_rd);
    }

#if CONFIG_PRED_FILTER
    // Keep track of the best mode irrespective of prediction filter state
    if (this_rd < best_overall_rd) {
      best_overall_rd = this_rd;
      best_filter_state = xd->mode_info_context->mbmi.pred_filter_enabled;
    }

    // Ignore modes where the prediction filter state doesn't
    // match the state signaled at the frame level
    if ((cm->pred_filter_mode == 2) ||
        (cm->pred_filter_mode ==
         xd->mode_info_context->mbmi.pred_filter_enabled)) {
#endif
      // Did this mode help.. i.e. is it the new best mode
      if (this_rd < best_rd || x->skip) {
        if (!mode_excluded) {
          // Note index of best mode so far
          best_mode_index = mode_index;

          if (this_mode <= B_PRED) {
            if (cpi->common.txfm_mode == ALLOW_8X8
                && this_mode != B_PRED
                && this_mode != I8X8_PRED)
              x->e_mbd.mode_info_context->mbmi.uv_mode = uv_intra_mode_8x8;
            else
              x->e_mbd.mode_info_context->mbmi.uv_mode = uv_intra_mode;
            /* required for left and above block mv */
            x->e_mbd.mode_info_context->mbmi.mv.as_int = 0;
          }

          other_cost +=
            ref_costs[x->e_mbd.mode_info_context->mbmi.ref_frame];

          /* Calculate the final y RD estimate for this mode */
          best_yrd = RDCOST(x->rdmult, x->rddiv, (rate2 - rate_uv - other_cost),
                            (distortion2 - distortion_uv));

          *returnrate = rate2;
          *returndistortion = distortion2;
          best_rd = this_rd;
          vpx_memcpy(&best_mbmode, &x->e_mbd.mode_info_context->mbmi, sizeof(MB_MODE_INFO));
          vpx_memcpy(&best_partition, x->partition_info, sizeof(PARTITION_INFO));

          if ((this_mode == B_PRED)
              || (this_mode == I8X8_PRED)
              || (this_mode == SPLITMV))
            for (i = 0; i < 16; i++) {
              best_bmodes[i] = x->e_mbd.block[i].bmi;
            }
        }

        // Testing this mode gave rise to an improvement in best error score.
        // Lower threshold a bit for next time
        cpi->rd_thresh_mult[mode_index] =
            (cpi->rd_thresh_mult[mode_index] >= (MIN_THRESHMULT + 2)) ?
            cpi->rd_thresh_mult[mode_index] - 2 : MIN_THRESHMULT;
        cpi->rd_threshes[mode_index] =
            (cpi->rd_baseline_thresh[mode_index] >> 7) *
            cpi->rd_thresh_mult[mode_index];
      }
      // If the mode did not help improve the best error case then raise the
      // threshold for testing that mode next time around.
      else {
        cpi->rd_thresh_mult[mode_index] += 4;

        if (cpi->rd_thresh_mult[mode_index] > MAX_THRESHMULT)
          cpi->rd_thresh_mult[mode_index] = MAX_THRESHMULT;

        cpi->rd_threshes[mode_index] = (cpi->rd_baseline_thresh[mode_index] >> 7) * cpi->rd_thresh_mult[mode_index];
      }

      /* keep record of best compound/single-only prediction */
      if (!disable_skip &&
          x->e_mbd.mode_info_context->mbmi.ref_frame != INTRA_FRAME) {
        int64_t single_rd, hybrid_rd;
        int single_rate, hybrid_rate;

        if (cpi->common.comp_pred_mode == HYBRID_PREDICTION) {
          single_rate = rate2 - compmode_cost;
          hybrid_rate = rate2;
        } else {
          single_rate = rate2;
          hybrid_rate = rate2 + compmode_cost;
        }

        single_rd = RDCOST(x->rdmult, x->rddiv, single_rate, distortion2);
        hybrid_rd = RDCOST(x->rdmult, x->rddiv, hybrid_rate, distortion2);

        if (x->e_mbd.mode_info_context->mbmi.second_ref_frame == INTRA_FRAME &&
            single_rd < best_single_rd) {
          best_single_rd = single_rd;
        } else if (x->e_mbd.mode_info_context->mbmi.second_ref_frame != INTRA_FRAME &&
                   single_rd < best_comp_rd) {
          best_comp_rd = single_rd;
        }
        if (hybrid_rd < best_hybrid_rd)
          best_hybrid_rd = hybrid_rd;
      }
#if CONFIG_PRED_FILTER
    }
#endif

    if (x->skip)
      break;
  }

#if CONFIG_PRED_FILTER
  // Update counts for prediction filter usage
  if (best_filter_state != 0)
    ++cpi->pred_filter_on_count;
  else
    ++cpi->pred_filter_off_count;
#endif
#if CONFIG_SWITCHABLE_INTERP
  if (cpi->common.mcomp_filter_type == SWITCHABLE &&
      best_mbmode.mode >= NEARESTMV &&
      best_mbmode.mode <= SPLITMV) {
    ++cpi->switchable_interp_count
        [get_pred_context(&cpi->common, xd, PRED_SWITCHABLE_INTERP)]
        [vp8_switchable_interp_map[best_mbmode.interp_filter]];
  }
#endif

  // Reduce the activation RD thresholds for the best choice mode
  if ((cpi->rd_baseline_thresh[best_mode_index] > 0) &&
      (cpi->rd_baseline_thresh[best_mode_index] < (INT_MAX >> 2))) {
    int best_adjustment = (cpi->rd_thresh_mult[best_mode_index] >> 2);

    cpi->rd_thresh_mult[best_mode_index] =
        (cpi->rd_thresh_mult[best_mode_index] >=
         (MIN_THRESHMULT + best_adjustment)) ?
        cpi->rd_thresh_mult[best_mode_index] - best_adjustment : MIN_THRESHMULT;
    cpi->rd_threshes[best_mode_index] =
        (cpi->rd_baseline_thresh[best_mode_index] >> 7) *
        cpi->rd_thresh_mult[best_mode_index];
  }

  // This code force Altref,0,0 and skip for the frame that overlays a
  // an alrtef unless Altref is filtered. However, this is unsafe if
  // segment level coding of ref frame or mode is enabled for this
  // segment.
  if (!segfeature_active(xd, segment_id, SEG_LVL_REF_FRAME) &&
      !segfeature_active(xd, segment_id, SEG_LVL_MODE) &&
      cpi->is_src_frame_alt_ref &&
      (cpi->oxcf.arnr_max_frames == 0) &&
      (best_mbmode.mode != ZEROMV || best_mbmode.ref_frame != ALTREF_FRAME)) {
    x->e_mbd.mode_info_context->mbmi.mode = ZEROMV;
    x->e_mbd.mode_info_context->mbmi.ref_frame = ALTREF_FRAME;
    x->e_mbd.mode_info_context->mbmi.mv.as_int = 0;
    x->e_mbd.mode_info_context->mbmi.uv_mode = DC_PRED;
    x->e_mbd.mode_info_context->mbmi.mb_skip_coeff =
      (cpi->common.mb_no_coeff_skip) ? 1 : 0;
    x->e_mbd.mode_info_context->mbmi.partitioning = 0;

    *best_single_rd_diff = *best_comp_rd_diff = *best_hybrid_rd_diff = 0;

    store_coding_context(x, xd->mb_index, best_mode_index, &best_partition,
                         &frame_best_ref_mv[xd->mode_info_context->mbmi.ref_frame],
                         &frame_best_ref_mv[xd->mode_info_context->mbmi.second_ref_frame]);
    return;
  }

  // macroblock modes
  vpx_memcpy(&x->e_mbd.mode_info_context->mbmi,
             &best_mbmode, sizeof(MB_MODE_INFO));
#if CONFIG_NEWBESTREFMV
  x->e_mbd.mode_info_context->mbmi.ref_mv =
    ref_mv[best_mbmode.ref_frame];
  x->e_mbd.mode_info_context->mbmi.second_ref_mv =
    ref_mv[best_mbmode.second_ref_frame];
#endif
  if (best_mbmode.mode == B_PRED) {
    for (i = 0; i < 16; i++) {
      xd->mode_info_context->bmi[i].as_mode = best_bmodes[i].as_mode;
      xd->block[i].bmi.as_mode = xd->mode_info_context->bmi[i].as_mode;
    }
  }

  if (best_mbmode.mode == I8X8_PRED)
    set_i8x8_block_modes(x, mode8x8);

  if (best_mbmode.mode == SPLITMV) {
    for (i = 0; i < 16; i++)
      xd->mode_info_context->bmi[i].as_mv.first.as_int = best_bmodes[i].as_mv.first.as_int;
    if (xd->mode_info_context->mbmi.second_ref_frame)
      for (i = 0; i < 16; i++)
        xd->mode_info_context->bmi[i].as_mv.second.as_int = best_bmodes[i].as_mv.second.as_int;

    vpx_memcpy(x->partition_info, &best_partition, sizeof(PARTITION_INFO));

    x->e_mbd.mode_info_context->mbmi.mv.as_int =
      x->partition_info->bmi[15].mv.as_int;
    x->e_mbd.mode_info_context->mbmi.second_mv.as_int =
      x->partition_info->bmi[15].second_mv.as_int;
  }

  if (best_single_rd == INT64_MAX)
    *best_single_rd_diff = INT_MIN;
  else
    *best_single_rd_diff = best_rd - best_single_rd;
  if (best_comp_rd == INT64_MAX)
    *best_comp_rd_diff = INT_MIN;
  else
    *best_comp_rd_diff   = best_rd - best_comp_rd;
  if (best_hybrid_rd == INT64_MAX)
    *best_hybrid_rd_diff = INT_MIN;
  else
    *best_hybrid_rd_diff = best_rd - best_hybrid_rd;

  store_coding_context(x, xd->mb_index, best_mode_index, &best_partition,
                       &frame_best_ref_mv[xd->mode_info_context->mbmi.ref_frame],
                       &frame_best_ref_mv[xd->mode_info_context->mbmi.second_ref_frame]);
}

int vp8_rd_pick_intra_mode(VP8_COMP *cpi, MACROBLOCK *x) {
  MACROBLOCKD *xd = &x->e_mbd;
  int64_t error4x4, error16x16;
#if CONFIG_COMP_INTRA_PRED
  int64_t error4x4d;
  int rate4x4d, dist4x4d;
#endif
  int rate4x4, rate16x16 = 0, rateuv;
  int dist4x4, dist16x16, distuv;
  int rate;
  int rate4x4_tokenonly = 0;
  int rate16x16_tokenonly = 0;
  int rateuv_tokenonly = 0;
  int64_t error8x8;
  int rate8x8_tokenonly=0;
  int rate8x8, dist8x8;
  int mode16x16;
  int mode8x8[2][4];

  xd->mode_info_context->mbmi.ref_frame = INTRA_FRAME;

  rd_pick_intra_mbuv_mode(cpi, x, &rateuv, &rateuv_tokenonly, &distuv);
  rate = rateuv;

  // current macroblock under rate-distortion optimization test loop
#if CONFIG_HYBRIDTRANSFORM
  xd->mode_info_context->mbmi.mode_rdopt = DC_PRED;
#endif

  error16x16 = rd_pick_intra16x16mby_mode(cpi, x, &rate16x16,
                                          &rate16x16_tokenonly, &dist16x16);
  mode16x16 = xd->mode_info_context->mbmi.mode;

#if CONFIG_HYBRIDTRANSFORM
  xd->mode_info_context->mbmi.mode_rdopt = I8X8_PRED;
#endif

  error8x8 = rd_pick_intra8x8mby_modes(cpi, x, &rate8x8, &rate8x8_tokenonly,
                                       &dist8x8, error16x16);
  mode8x8[0][0]= xd->mode_info_context->bmi[0].as_mode.first;
  mode8x8[0][1]= xd->mode_info_context->bmi[2].as_mode.first;
  mode8x8[0][2]= xd->mode_info_context->bmi[8].as_mode.first;
  mode8x8[0][3]= xd->mode_info_context->bmi[10].as_mode.first;
#if CONFIG_COMP_INTRA_PRED
  mode8x8[1][0] = xd->mode_info_context->bmi[0].as_mode.second;
  mode8x8[1][1] = xd->mode_info_context->bmi[2].as_mode.second;
  mode8x8[1][2] = xd->mode_info_context->bmi[8].as_mode.second;
  mode8x8[1][3] = xd->mode_info_context->bmi[10].as_mode.second;
#endif

#if CONFIG_HYBRIDTRANSFORM
  xd->mode_info_context->mbmi.mode_rdopt = B_PRED;
#endif

  error4x4 = rd_pick_intra4x4mby_modes(cpi, x,
                                       &rate4x4, &rate4x4_tokenonly,
                                       &dist4x4, error16x16,
#if CONFIG_COMP_INTRA_PRED
                                       0,
#endif
                                       0);
#if CONFIG_COMP_INTRA_PRED
  error4x4d = rd_pick_intra4x4mby_modes(cpi, x,
                                        &rate4x4d, &rate4x4_tokenonly,
                                        &dist4x4d, error16x16, 1, 0);
#endif

  if (error8x8 > error16x16) {
    if (error4x4 < error16x16) {
#if CONFIG_COMP_INTRA_PRED
      rate += (error4x4d < error4x4) ? rate4x4d : rate4x4;
      if (error4x4d >= error4x4) // FIXME save original modes etc.
        error4x4 = rd_pick_intra4x4mby_modes(cpi, x, &rate4x4,
                                             &rate4x4_tokenonly,
                                             &dist4x4, error16x16, 0,
                                             cpi->update_context);
#else
      rate += rate4x4;
#endif
      xd->mode_info_context->mbmi.mode = B_PRED;
    } else {
      xd->mode_info_context->mbmi.mode = mode16x16;
      rate += rate16x16;
    }
  } else {
    if (error4x4 < error8x8) {
#if CONFIG_COMP_INTRA_PRED
      rate += (error4x4d < error4x4) ? rate4x4d : rate4x4;
      if (error4x4d >= error4x4) // FIXME save original modes etc.
        error4x4 = rd_pick_intra4x4mby_modes(cpi, x, &rate4x4,
                                             &rate4x4_tokenonly,
                                             &dist4x4, error16x16, 0,
                                             cpi->update_context);
#else
      rate += rate4x4;
#endif
      xd->mode_info_context->mbmi.mode = B_PRED;
    } else {
      xd->mode_info_context->mbmi.mode = I8X8_PRED;
      set_i8x8_block_modes(x, mode8x8);
      rate += rate8x8;
    }
  }
  return rate;
}

int vp8cx_pick_mode_inter_macroblock(VP8_COMP *cpi, MACROBLOCK *x,
                                     int recon_yoffset, int recon_uvoffset) {
  VP8_COMMON *cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  int rate;
  int distortion;
  int64_t intra_error = 0;
  unsigned char *segment_id = &xd->mode_info_context->mbmi.segment_id;

  if (xd->segmentation_enabled)
    x->encode_breakout = cpi->segment_encode_breakout[*segment_id];
  else
    x->encode_breakout = cpi->oxcf.encode_breakout;

  // if (cpi->sf.RD)
  // For now this codebase is limited to a single rd encode path
  {
    int zbin_mode_boost_enabled = cpi->zbin_mode_boost_enabled;
    int64_t single, compound, hybrid;

    vp8_rd_pick_inter_mode(cpi, x, recon_yoffset, recon_uvoffset, &rate,
                           &distortion, &intra_error, &single, &compound,
                           &hybrid);

    // TODO Save these to add in only if MB coding mode is selected?
    cpi->rd_single_diff += single;
    cpi->rd_comp_diff   += compound;
    cpi->rd_hybrid_diff += hybrid;
    if (xd->mode_info_context->mbmi.ref_frame) {
      unsigned char pred_context;

      pred_context = get_pred_context(cm, xd, PRED_COMP);

      if (xd->mode_info_context->mbmi.second_ref_frame == INTRA_FRAME)
        cpi->single_pred_count[pred_context]++;
      else
        cpi->comp_pred_count[pred_context]++;
    }

    /* restore cpi->zbin_mode_boost_enabled */
    cpi->zbin_mode_boost_enabled = zbin_mode_boost_enabled;
  }
  // else
  // The non rd encode path has been deleted from this code base
  // to simplify development
  //    vp8_pick_inter_mode

  // Store metrics so they can be added in to totals if this mode is picked
  x->mb_context[xd->mb_index].distortion  = distortion;
  x->mb_context[xd->mb_index].intra_error = intra_error;

  return rate;
}
