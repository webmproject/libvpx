
/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <limits.h>

#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_treecoder.h"

// TBD prediction functions for various bitstream signals

// Returns a context number for the given MB prediction signal
unsigned char vp9_get_pred_context(const VP9_COMMON *const cm,
                                   const MACROBLOCKD *const xd,
                                   PRED_ID pred_id) {
  int pred_context;
  const MODE_INFO *const mi = xd->mode_info_context;
  const MODE_INFO *const above_mi = mi - cm->mode_info_stride;
  const MODE_INFO *const left_mi = mi - 1;
  const int left_in_image = xd->left_available && left_mi->mbmi.mb_in_image;
  const int above_in_image = xd->up_available && above_mi->mbmi.mb_in_image;
  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries correpsonding to real macroblocks.
  // The prediction flags in these dummy entries are initialised to 0.
  switch (pred_id) {
    case PRED_SEG_ID:
      pred_context = above_mi->mbmi.seg_id_predicted;
      if (xd->left_available)
        pred_context += left_mi->mbmi.seg_id_predicted;
      break;

    case PRED_MBSKIP:
      pred_context = above_mi->mbmi.mb_skip_coeff;
      if (xd->left_available)
        pred_context += left_mi->mbmi.mb_skip_coeff;
      break;

    case PRED_SWITCHABLE_INTERP: {
      // left
      const int left_mv_pred = is_inter_mode(left_mi->mbmi.mode);
      const int left_interp = left_in_image && left_mv_pred ?
                    vp9_switchable_interp_map[left_mi->mbmi.interp_filter] :
                    VP9_SWITCHABLE_FILTERS;

      // above
      const int above_mv_pred = is_inter_mode(above_mi->mbmi.mode);
      const int above_interp = above_in_image && above_mv_pred ?
                    vp9_switchable_interp_map[above_mi->mbmi.interp_filter] :
                    VP9_SWITCHABLE_FILTERS;

      assert(left_interp != -1);
      assert(above_interp != -1);

      if (left_interp == above_interp)
        pred_context = left_interp;
      else if (left_interp == VP9_SWITCHABLE_FILTERS &&
               above_interp != VP9_SWITCHABLE_FILTERS)
         pred_context = above_interp;
      else if (left_interp != VP9_SWITCHABLE_FILTERS &&
               above_interp == VP9_SWITCHABLE_FILTERS)
        pred_context = left_interp;
      else
        pred_context = VP9_SWITCHABLE_FILTERS;

      break;
    }

    case PRED_INTRA_INTER: {
      if (above_in_image && left_in_image) {  // both edges available
        if (left_mi->mbmi.ref_frame[0] == INTRA_FRAME &&
            above_mi->mbmi.ref_frame[0] == INTRA_FRAME) {  // intra/intra (3)
          pred_context = 3;
        } else {  // intra/inter (1) or inter/inter (0)
          pred_context = left_mi->mbmi.ref_frame[0] == INTRA_FRAME ||
                         above_mi->mbmi.ref_frame[0] == INTRA_FRAME;
        }
      } else if (above_in_image || left_in_image) {  // one edge available
        const MODE_INFO *edge = above_in_image ? above_mi : left_mi;

        // inter: 0, intra: 2
        pred_context = 2 * (edge->mbmi.ref_frame[0] == INTRA_FRAME);
      } else {
        pred_context = 0;
      }
      assert(pred_context >= 0 && pred_context < INTRA_INTER_CONTEXTS);
      break;
    }

    case PRED_COMP_INTER_INTER: {
      if (above_in_image && left_in_image) {  // both edges available
        if (above_mi->mbmi.ref_frame[1] <= INTRA_FRAME &&
            left_mi->mbmi.ref_frame[1] <= INTRA_FRAME) {
          // neither edge uses comp pred (0/1)
          pred_context = ((above_mi->mbmi.ref_frame[0] == cm->comp_fixed_ref) ^
                          (left_mi->mbmi.ref_frame[0] == cm->comp_fixed_ref));
        } else if (above_mi->mbmi.ref_frame[1] <= INTRA_FRAME) {
          // one of two edges uses comp pred (2/3)
          pred_context = 2 +
              (above_mi->mbmi.ref_frame[0] == cm->comp_fixed_ref ||
               above_mi->mbmi.ref_frame[0] == INTRA_FRAME);
        } else if (left_mi->mbmi.ref_frame[1] <= INTRA_FRAME) {
          // one of two edges uses comp pred (2/3)
          pred_context = 2 +
              (left_mi->mbmi.ref_frame[0] == cm->comp_fixed_ref ||
               left_mi->mbmi.ref_frame[0] == INTRA_FRAME);
        } else {  // both edges use comp pred (4)
          pred_context = 4;
        }
      } else if (above_in_image || left_in_image) {  // one edge available
        const MODE_INFO *edge = above_in_image ? above_mi : left_mi;

        if (edge->mbmi.ref_frame[1] <= INTRA_FRAME) {
          // edge does not use comp pred (0/1)
          pred_context = edge->mbmi.ref_frame[0] == cm->comp_fixed_ref;
        } else {  // edge uses comp pred (3)
          pred_context = 3;
        }
      } else {  // no edges available (1)
        pred_context = 1;
      }
      assert(pred_context >= 0 && pred_context < COMP_INTER_CONTEXTS);
      break;
    }

    case PRED_COMP_REF_P: {
      const int fix_ref_idx = cm->ref_frame_sign_bias[cm->comp_fixed_ref];
      const int var_ref_idx = !fix_ref_idx;

      if (above_in_image && left_in_image) {  // both edges available
        if (above_mi->mbmi.ref_frame[0] == INTRA_FRAME &&
            left_mi->mbmi.ref_frame[0] == INTRA_FRAME) {  // intra/intra (2)
          pred_context = 2;
        } else if (above_mi->mbmi.ref_frame[0] == INTRA_FRAME ||
                   left_mi->mbmi.ref_frame[0] == INTRA_FRAME) {  // intra/inter
          const MODE_INFO *edge = above_mi->mbmi.ref_frame[0] == INTRA_FRAME ?
                                  left_mi : above_mi;

          if (edge->mbmi.ref_frame[1] <= INTRA_FRAME) {  // single pred (1/3)
            pred_context = 1 +
                2 * (edge->mbmi.ref_frame[0] != cm->comp_var_ref[1]);
          } else {  // comp pred (1/3)
            pred_context = 1 +
                2 * (edge->mbmi.ref_frame[var_ref_idx] != cm->comp_var_ref[1]);
          }
        } else {  // inter/inter
          int l_sg = left_mi->mbmi.ref_frame[1] <= INTRA_FRAME;
          int a_sg = above_mi->mbmi.ref_frame[1] <= INTRA_FRAME;
          MV_REFERENCE_FRAME vrfa = a_sg ? above_mi->mbmi.ref_frame[0] :
              above_mi->mbmi.ref_frame[var_ref_idx];
          MV_REFERENCE_FRAME vrfl = l_sg ? left_mi->mbmi.ref_frame[0] :
              left_mi->mbmi.ref_frame[var_ref_idx];

          if (vrfa == vrfl && cm->comp_var_ref[1] == vrfa) {
            pred_context = 0;
          } else if (l_sg && a_sg) {  // single/single
            if ((vrfa == cm->comp_fixed_ref && vrfl == cm->comp_var_ref[0]) ||
                (vrfl == cm->comp_fixed_ref && vrfa == cm->comp_var_ref[0])) {
              pred_context = 4;
            } else if (vrfa == vrfl) {
              pred_context = 3;
            } else {
              pred_context = 1;
            }
          } else if (l_sg || a_sg) {  // single/comp
            MV_REFERENCE_FRAME vrfc = l_sg ? vrfa : vrfl;
            MV_REFERENCE_FRAME rfs = a_sg ? vrfa : vrfl;

            if (vrfc == cm->comp_var_ref[1] && rfs != cm->comp_var_ref[1]) {
              pred_context = 1;
            } else if (rfs == cm->comp_var_ref[1] &&
                       vrfc != cm->comp_var_ref[1]) {
              pred_context = 2;
            } else {
              pred_context = 4;
            }
          } else if (vrfa == vrfl) {  // comp/comp
            pred_context = 4;
          } else {
            pred_context = 2;
          }
        }
      } else if (above_in_image || left_in_image) {  // one edge available
        const MODE_INFO *edge = above_in_image ? above_mi : left_mi;

        if (edge->mbmi.ref_frame[0] == INTRA_FRAME) {
          pred_context = 2;
        } else if (edge->mbmi.ref_frame[1] > INTRA_FRAME) {
          pred_context =
              4 * (edge->mbmi.ref_frame[var_ref_idx] != cm->comp_var_ref[1]);
        } else {
          pred_context = 3 * (edge->mbmi.ref_frame[0] != cm->comp_var_ref[1]);
        }
      } else {  // no edges available (2)
        pred_context = 2;
      }
      assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
      break;
    }

    case PRED_SINGLE_REF_P1: {
      if (above_in_image && left_in_image) {  // both edges available
        if (above_mi->mbmi.ref_frame[0] == INTRA_FRAME &&
            left_mi->mbmi.ref_frame[0] == INTRA_FRAME) {
          pred_context = 2;
        } else if (above_mi->mbmi.ref_frame[0] == INTRA_FRAME ||
                   left_mi->mbmi.ref_frame[0] == INTRA_FRAME) {
          const MODE_INFO *edge = above_mi->mbmi.ref_frame[0] == INTRA_FRAME ?
                                  left_mi : above_mi;

          if (edge->mbmi.ref_frame[1] <= INTRA_FRAME) {
            pred_context = 4 * (edge->mbmi.ref_frame[0] == LAST_FRAME);
          } else {
            pred_context = 1 + (edge->mbmi.ref_frame[0] == LAST_FRAME ||
                                edge->mbmi.ref_frame[1] == LAST_FRAME);
          }
        } else if (above_mi->mbmi.ref_frame[1] <= INTRA_FRAME &&
                   left_mi->mbmi.ref_frame[1] <= INTRA_FRAME) {
          pred_context = 2 * (above_mi->mbmi.ref_frame[0] == LAST_FRAME) +
                         2 * (left_mi->mbmi.ref_frame[0] == LAST_FRAME);
        } else if (above_mi->mbmi.ref_frame[1] > INTRA_FRAME &&
                   left_mi->mbmi.ref_frame[1] > INTRA_FRAME) {
          pred_context = 1 + (above_mi->mbmi.ref_frame[0] == LAST_FRAME ||
                              above_mi->mbmi.ref_frame[1] == LAST_FRAME ||
                              left_mi->mbmi.ref_frame[0] == LAST_FRAME ||
                              left_mi->mbmi.ref_frame[1] == LAST_FRAME);
        } else {
          MV_REFERENCE_FRAME rfs = above_mi->mbmi.ref_frame[1] <= INTRA_FRAME ?
              above_mi->mbmi.ref_frame[0] : left_mi->mbmi.ref_frame[0];
          MV_REFERENCE_FRAME crf1 = above_mi->mbmi.ref_frame[1] > INTRA_FRAME ?
              above_mi->mbmi.ref_frame[0] : left_mi->mbmi.ref_frame[0];
          MV_REFERENCE_FRAME crf2 = above_mi->mbmi.ref_frame[1] > INTRA_FRAME ?
              above_mi->mbmi.ref_frame[1] : left_mi->mbmi.ref_frame[1];

          if (rfs == LAST_FRAME) {
            pred_context = 3 + (crf1 == LAST_FRAME || crf2 == LAST_FRAME);
          } else {
            pred_context = crf1 == LAST_FRAME || crf2 == LAST_FRAME;
          }
        }
      } else if (above_in_image || left_in_image) {  // one edge available
        const MODE_INFO *edge = above_in_image ? above_mi : left_mi;

        if (edge->mbmi.ref_frame[0] == INTRA_FRAME) {
          pred_context = 2;
        } else if (edge->mbmi.ref_frame[1] <= INTRA_FRAME) {
          pred_context = 4 * (edge->mbmi.ref_frame[0] == LAST_FRAME);
        } else {
          pred_context = 1 + (edge->mbmi.ref_frame[0] == LAST_FRAME ||
                              edge->mbmi.ref_frame[1] == LAST_FRAME);
        }
      } else {  // no edges available (2)
        pred_context = 2;
      }
      assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
      break;
    }

    case PRED_SINGLE_REF_P2: {
      if (above_in_image && left_in_image) {  // both edges available
        if (above_mi->mbmi.ref_frame[0] == INTRA_FRAME &&
            left_mi->mbmi.ref_frame[0] == INTRA_FRAME) {
          pred_context = 2;
        } else if (above_mi->mbmi.ref_frame[0] == INTRA_FRAME ||
                   left_mi->mbmi.ref_frame[0] == INTRA_FRAME) {
          const MODE_INFO *edge = above_mi->mbmi.ref_frame[0] == INTRA_FRAME ?
                                  left_mi : above_mi;

          if (edge->mbmi.ref_frame[1] <= INTRA_FRAME) {
            if (edge->mbmi.ref_frame[0] == LAST_FRAME) {
              pred_context = 3;
            } else {
              pred_context = 4 * (edge->mbmi.ref_frame[0] == GOLDEN_FRAME);
            }
          } else {
            pred_context = 1 + 2 * (edge->mbmi.ref_frame[0] == GOLDEN_FRAME ||
                                    edge->mbmi.ref_frame[1] == GOLDEN_FRAME);
          }
        } else if (above_mi->mbmi.ref_frame[1] <= INTRA_FRAME &&
                   left_mi->mbmi.ref_frame[1] <= INTRA_FRAME) {
          if (above_mi->mbmi.ref_frame[0] == LAST_FRAME &&
              left_mi->mbmi.ref_frame[0] == LAST_FRAME) {
            pred_context = 3;
          } else if (above_mi->mbmi.ref_frame[0] == LAST_FRAME ||
                     left_mi->mbmi.ref_frame[0] == LAST_FRAME) {
            const MODE_INFO *edge = above_mi->mbmi.ref_frame[0] == LAST_FRAME ?
                                    left_mi : above_mi;

            pred_context = 4 * (edge->mbmi.ref_frame[0] == GOLDEN_FRAME);
          } else {
            pred_context = 2 * (above_mi->mbmi.ref_frame[0] == GOLDEN_FRAME) +
                           2 * (left_mi->mbmi.ref_frame[0] == GOLDEN_FRAME);
          }
        } else if (above_mi->mbmi.ref_frame[1] > INTRA_FRAME &&
                   left_mi->mbmi.ref_frame[1] > INTRA_FRAME) {
          if (above_mi->mbmi.ref_frame[0] == left_mi->mbmi.ref_frame[0] &&
              above_mi->mbmi.ref_frame[1] == left_mi->mbmi.ref_frame[1]) {
            pred_context = 3 * (above_mi->mbmi.ref_frame[0] == GOLDEN_FRAME ||
                                above_mi->mbmi.ref_frame[1] == GOLDEN_FRAME ||
                                left_mi->mbmi.ref_frame[0] == GOLDEN_FRAME ||
                                left_mi->mbmi.ref_frame[1] == GOLDEN_FRAME);
          } else {
            pred_context = 2;
          }
        } else {
          MV_REFERENCE_FRAME rfs = above_mi->mbmi.ref_frame[1] <= INTRA_FRAME ?
              above_mi->mbmi.ref_frame[0] : left_mi->mbmi.ref_frame[0];
          MV_REFERENCE_FRAME crf1 = above_mi->mbmi.ref_frame[1] > INTRA_FRAME ?
              above_mi->mbmi.ref_frame[0] : left_mi->mbmi.ref_frame[0];
          MV_REFERENCE_FRAME crf2 = above_mi->mbmi.ref_frame[1] > INTRA_FRAME ?
              above_mi->mbmi.ref_frame[1] : left_mi->mbmi.ref_frame[1];

          if (rfs == GOLDEN_FRAME) {
            pred_context = 3 + (crf1 == GOLDEN_FRAME || crf2 == GOLDEN_FRAME);
          } else if (rfs == ALTREF_FRAME) {
            pred_context = crf1 == GOLDEN_FRAME || crf2 == GOLDEN_FRAME;
          } else {
            pred_context =
                1 + 2 * (crf1 == GOLDEN_FRAME || crf2 == GOLDEN_FRAME);
          }
        }
      } else if (above_in_image || left_in_image) {  // one edge available
        const MODE_INFO *edge = above_in_image ? above_mi : left_mi;

        if (edge->mbmi.ref_frame[0] == INTRA_FRAME ||
            (edge->mbmi.ref_frame[0] == LAST_FRAME &&
             edge->mbmi.ref_frame[1] <= INTRA_FRAME)) {
          pred_context = 2;
        } else if (edge->mbmi.ref_frame[1] <= INTRA_FRAME) {
          pred_context = 4 * (edge->mbmi.ref_frame[0] == GOLDEN_FRAME);
        } else {
          pred_context = 3 * (edge->mbmi.ref_frame[0] == GOLDEN_FRAME ||
                              edge->mbmi.ref_frame[1] == GOLDEN_FRAME);
        }
      } else {  // no edges available (2)
        pred_context = 2;
      }
      assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
      break;
    }

    case PRED_TX_SIZE: {
      int above_context, left_context;
      int max_tx_size;
      if (mi->mbmi.sb_type < BLOCK_SIZE_SB8X8)
        max_tx_size = TX_4X4;
      else if (mi->mbmi.sb_type < BLOCK_SIZE_MB16X16)
        max_tx_size = TX_8X8;
      else if (mi->mbmi.sb_type < BLOCK_SIZE_SB32X32)
        max_tx_size = TX_16X16;
      else
        max_tx_size = TX_32X32;
      above_context = left_context = max_tx_size;
      if (above_in_image) {
        above_context = (above_mi->mbmi.mb_skip_coeff ?
                         max_tx_size : above_mi->mbmi.txfm_size);
      }
      if (left_in_image) {
        left_context = (left_mi->mbmi.mb_skip_coeff ?
                        max_tx_size : left_mi->mbmi.txfm_size);
      }
      if (!left_in_image) {
        left_context = above_context;
      }
      if (!above_in_image) {
        above_context = left_context;
      }
      pred_context = (above_context + left_context > max_tx_size);
      break;
    }

    default:
      assert(0);
      pred_context = 0;  // *** add error trap code.
      break;
  }

  return pred_context;
}

// This function returns a context probability for coding a given
// prediction signal
vp9_prob vp9_get_pred_prob(const VP9_COMMON *const cm,
                          const MACROBLOCKD *const xd,
                          PRED_ID pred_id) {
  const int pred_context = vp9_get_pred_context(cm, xd, pred_id);

  switch (pred_id) {
    case PRED_SEG_ID:
      return cm->segment_pred_probs[pred_context];
    case PRED_MBSKIP:
      return cm->fc.mbskip_probs[pred_context];
    case PRED_INTRA_INTER:
      return cm->fc.intra_inter_prob[pred_context];
    case PRED_COMP_INTER_INTER:
      return cm->fc.comp_inter_prob[pred_context];
    case PRED_COMP_REF_P:
      return cm->fc.comp_ref_prob[pred_context];
    case PRED_SINGLE_REF_P1:
      return cm->fc.single_ref_prob[pred_context][0];
    case PRED_SINGLE_REF_P2:
      return cm->fc.single_ref_prob[pred_context][1];
    default:
      assert(0);
      return 128;  // *** add error trap code.
  }
}

// This function returns a context probability ptr for coding a given
// prediction signal
const vp9_prob *vp9_get_pred_probs(const VP9_COMMON *const cm,
                                   const MACROBLOCKD *const xd,
                                   PRED_ID pred_id) {
  const MODE_INFO *const mi = xd->mode_info_context;
  const int pred_context = vp9_get_pred_context(cm, xd, pred_id);

  switch (pred_id) {
    case PRED_SWITCHABLE_INTERP:
      return &cm->fc.switchable_interp_prob[pred_context][0];

    case PRED_TX_SIZE:
      if (mi->mbmi.sb_type < BLOCK_SIZE_MB16X16)
        return cm->fc.tx_probs_8x8p[pred_context];
      else if (mi->mbmi.sb_type < BLOCK_SIZE_SB32X32)
        return cm->fc.tx_probs_16x16p[pred_context];
      else
        return cm->fc.tx_probs_32x32p[pred_context];

    default:
      assert(0);
      return NULL;  // *** add error trap code.
  }
}

// This function returns the status of the given prediction signal.
// I.e. is the predicted value for the given signal correct.
unsigned char vp9_get_pred_flag(const MACROBLOCKD *const xd,
                                PRED_ID pred_id) {
  switch (pred_id) {
    case PRED_SEG_ID:
      return xd->mode_info_context->mbmi.seg_id_predicted;
    case PRED_MBSKIP:
      return xd->mode_info_context->mbmi.mb_skip_coeff;
    default:
      assert(0);
      return 0;  // *** add error trap code.
  }
}

// This function sets the status of the given prediction signal.
// I.e. is the predicted value for the given signal correct.
void vp9_set_pred_flag(MACROBLOCKD *const xd,
                       PRED_ID pred_id,
                       unsigned char pred_flag) {
  const int mis = xd->mode_info_stride;
  BLOCK_SIZE_TYPE bsize = xd->mode_info_context->mbmi.sb_type;
  const int bh = 1 << mi_height_log2(bsize);
  const int bw = 1 << mi_width_log2(bsize);
#define sub(a, b) (b) < 0 ? (a) + (b) : (a)
  const int x_mis = sub(bw, xd->mb_to_right_edge >> (3 + LOG2_MI_SIZE));
  const int y_mis = sub(bh, xd->mb_to_bottom_edge >> (3 + LOG2_MI_SIZE));
#undef sub
  int x, y;

  switch (pred_id) {
    case PRED_SEG_ID:
      for (y = 0; y < y_mis; y++) {
        for (x = 0; x < x_mis; x++) {
          xd->mode_info_context[y * mis + x].mbmi.seg_id_predicted = pred_flag;
        }
      }
      break;

    case PRED_MBSKIP:
      for (y = 0; y < y_mis; y++) {
        for (x = 0; x < x_mis; x++) {
          xd->mode_info_context[y * mis + x].mbmi.mb_skip_coeff = pred_flag;
        }
      }
      break;

    default:
      assert(0);
      // *** add error trap code.
      break;
  }
}


// The following contain the guts of the prediction code used to
// peredict various bitstream signals.

// Macroblock segment id prediction function
int vp9_get_pred_mi_segid(VP9_COMMON *cm, BLOCK_SIZE_TYPE sb_type,
                          int mi_row, int mi_col) {
  const int mi_index = mi_row * cm->mi_cols + mi_col;
  const int bw = 1 << mi_width_log2(sb_type);
  const int bh = 1 << mi_height_log2(sb_type);
  const int ymis = MIN(cm->mi_rows - mi_row, bh);
  const int xmis = MIN(cm->mi_cols - mi_col, bw);
  int segment_id = INT_MAX;
  int x, y;

  for (y = 0; y < ymis; y++) {
    for (x = 0; x < xmis; x++) {
      const int index = mi_index + (y * cm->mi_cols + x);
      segment_id = MIN(segment_id, cm->last_frame_seg_map[index]);
    }
  }
  return segment_id;
}
