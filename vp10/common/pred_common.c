
/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp10/common/common.h"
#include "vp10/common/pred_common.h"
#include "vp10/common/seg_common.h"

// Returns a context number for the given MB prediction signal
int vp10_get_pred_context_switchable_interp(const MACROBLOCKD *xd) {
  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries corresponding to real macroblocks.
  // The prediction flags in these dummy entries are initialized to 0.
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int left_type = xd->left_available && is_inter_block(left_mbmi) ?
      left_mbmi->interp_filter : SWITCHABLE_FILTERS;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const int above_type = xd->up_available && is_inter_block(above_mbmi) ?
      above_mbmi->interp_filter : SWITCHABLE_FILTERS;

  if (left_type == above_type)
    return left_type;
  else if (left_type == SWITCHABLE_FILTERS && above_type != SWITCHABLE_FILTERS)
    return above_type;
  else if (left_type != SWITCHABLE_FILTERS && above_type == SWITCHABLE_FILTERS)
    return left_type;
  else
    return SWITCHABLE_FILTERS;
}

#if CONFIG_EXT_INTRA
int vp10_get_pred_context_intra_interp(const MACROBLOCKD *xd) {
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  int left_type = INTRA_FILTERS, above_type = INTRA_FILTERS;

  if (xd->left_available && left_mbmi->sb_type >= BLOCK_8X8) {
    PREDICTION_MODE mode = left_mbmi->mode;
    if (is_inter_block(left_mbmi)) {
      switch (left_mbmi->interp_filter) {
        case EIGHTTAP_REGULAR:
          left_type = INTRA_FILTER_8TAP;
          break;
        case EIGHTTAP_SMOOTH:
          left_type = INTRA_FILTER_8TAP_SMOOTH;
          break;
        case MULTITAP_SHARP:
          left_type = INTRA_FILTER_8TAP_SHARP;
          break;
        case BILINEAR:
          left_type = INTRA_FILTERS;
          break;
        default:
          break;
      }
    } else {
      if (mode != DC_PRED && mode != TM_PRED) {
        int p_angle;
        p_angle = mode_to_angle_map[mode] +
            left_mbmi->angle_delta[0] * ANGLE_STEP;
        if (pick_intra_filter(p_angle)) {
          left_type = left_mbmi->intra_filter;
        }
      }
    }
  }

  if (xd->up_available && above_mbmi->sb_type >= BLOCK_8X8) {
    if (is_inter_block(above_mbmi)) {
      switch (above_mbmi->interp_filter) {
        case EIGHTTAP_REGULAR:
          above_type = INTRA_FILTER_8TAP;
          break;
        case EIGHTTAP_SMOOTH:
          above_type = INTRA_FILTER_8TAP_SMOOTH;
          break;
        case MULTITAP_SHARP:
          above_type = INTRA_FILTER_8TAP_SHARP;
          break;
        case BILINEAR:
          above_type = INTRA_FILTERS;
          break;
        default:
          break;
      }
    } else {
      PREDICTION_MODE mode = above_mbmi->mode;
      if (mode != DC_PRED && mode != TM_PRED) {
        int p_angle;
        p_angle = mode_to_angle_map[mode] +
            above_mbmi->angle_delta[0] * ANGLE_STEP;
        if (pick_intra_filter(p_angle)) {
          above_type = above_mbmi->intra_filter;
        }
      }
    }
  }

  if (left_type == above_type)
    return left_type;
  else if (left_type == INTRA_FILTERS && above_type != INTRA_FILTERS)
    return above_type;
  else if (left_type != INTRA_FILTERS && above_type == INTRA_FILTERS)
    return left_type;
  else
    return INTRA_FILTERS;
}
#endif  // CONFIG_EXT_INTRA

// The mode info data structure has a one element border above and to the
// left of the entries corresponding to real macroblocks.
// The prediction flags in these dummy entries are initialized to 0.
// 0 - inter/inter, inter/--, --/inter, --/--
// 1 - intra/inter, inter/intra
// 2 - intra/--, --/intra
// 3 - intra/intra
int vp10_get_intra_inter_context(const MACROBLOCKD *xd) {
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int has_above = xd->up_available;
  const int has_left = xd->left_available;

  if (has_above && has_left) {  // both edges available
    const int above_intra = !is_inter_block(above_mbmi);
    const int left_intra = !is_inter_block(left_mbmi);
    return left_intra && above_intra ? 3
                                     : left_intra || above_intra;
  } else if (has_above || has_left) {  // one edge available
    return 2 * !is_inter_block(has_above ? above_mbmi : left_mbmi);
  } else {
    return 0;
  }
}

int vp10_get_reference_mode_context(const VP10_COMMON *cm,
                                   const MACROBLOCKD *xd) {
  int ctx;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int has_above = xd->up_available;
  const int has_left = xd->left_available;
  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries corresponding to real macroblocks.
  // The prediction flags in these dummy entries are initialized to 0.
  if (has_above && has_left) {  // both edges available
    if (!has_second_ref(above_mbmi) && !has_second_ref(left_mbmi))
      // neither edge uses comp pred (0/1)
      ctx = (above_mbmi->ref_frame[0] == cm->comp_fixed_ref) ^
            (left_mbmi->ref_frame[0] == cm->comp_fixed_ref);
    else if (!has_second_ref(above_mbmi))
      // one of two edges uses comp pred (2/3)
      ctx = 2 + (above_mbmi->ref_frame[0] == cm->comp_fixed_ref ||
                 !is_inter_block(above_mbmi));
    else if (!has_second_ref(left_mbmi))
      // one of two edges uses comp pred (2/3)
      ctx = 2 + (left_mbmi->ref_frame[0] == cm->comp_fixed_ref ||
                 !is_inter_block(left_mbmi));
    else  // both edges use comp pred (4)
      ctx = 4;
  } else if (has_above || has_left) {  // one edge available
    const MB_MODE_INFO *edge_mbmi = has_above ? above_mbmi : left_mbmi;

    if (!has_second_ref(edge_mbmi))
      // edge does not use comp pred (0/1)
      ctx = edge_mbmi->ref_frame[0] == cm->comp_fixed_ref;
    else
      // edge uses comp pred (3)
      ctx = 3;
  } else {  // no edges available (1)
    ctx = 1;
  }
  assert(ctx >= 0 && ctx < COMP_INTER_CONTEXTS);
  return ctx;
}

#if CONFIG_EXT_REFS

// TODO(zoeliu): Future work will be conducted to optimize the context design
// for the coding of the reference frames.

#define CHECK_LAST_OR_LAST2(ref_frame) \
  ((ref_frame == LAST_FRAME) || (ref_frame == LAST2_FRAME))

#define CHECK_GOLDEN_LAST3_LAST4(ref_frame) \
  ((ref_frame == GOLDEN_FRAME) || (ref_frame == LAST3_FRAME) || \
  (ref_frame == LAST4_FRAME))

// Returns a context number for the given MB prediction signal
// Signal the first reference frame for a compound mode is either
// GOLDEN/LAST3/LAST4, or LAST/LAST2.
//
// NOTE(zoeliu): The probability of ref_frame[0] is either
//               GOLDEN_FRAME/LAST3_FRAME/LAST4_FRAME.
int vp10_get_pred_context_comp_ref_p(const VP10_COMMON *cm,
                                     const MACROBLOCKD *xd) {
  int pred_context;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int above_in_image = xd->up_available;
  const int left_in_image = xd->left_available;

  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries correpsonding to real macroblocks.
  // The prediction flags in these dummy entries are initialised to 0.
  const int fix_ref_idx = cm->ref_frame_sign_bias[cm->comp_fixed_ref];
  const int var_ref_idx = !fix_ref_idx;

  if (above_in_image && left_in_image) {  // both edges available
    const int above_intra = !is_inter_block(above_mbmi);
    const int left_intra = !is_inter_block(left_mbmi);

    if (above_intra && left_intra) {  // intra/intra (2)
      pred_context = 2;
    } else if (above_intra || left_intra) {  // intra/inter
      const MB_MODE_INFO *edge_mbmi = above_intra ? left_mbmi : above_mbmi;

      if (!has_second_ref(edge_mbmi))  // single pred (1/3)
        pred_context = 1 +
            2 * (!CHECK_GOLDEN_LAST3_LAST4(edge_mbmi->ref_frame[0]));
      else  // comp pred (1/3)
        pred_context = 1 +
            2 * (!CHECK_GOLDEN_LAST3_LAST4(edge_mbmi->ref_frame[var_ref_idx]));
    } else {  // inter/inter
      const int l_sg = !has_second_ref(left_mbmi);
      const int a_sg = !has_second_ref(above_mbmi);
      const MV_REFERENCE_FRAME vrfa = a_sg ? above_mbmi->ref_frame[0]
                                           : above_mbmi->ref_frame[var_ref_idx];
      const MV_REFERENCE_FRAME vrfl = l_sg ? left_mbmi->ref_frame[0]
                                           : left_mbmi->ref_frame[var_ref_idx];

      if (vrfa == vrfl && CHECK_GOLDEN_LAST3_LAST4(vrfa)) {
        pred_context = 0;
      } else if (l_sg && a_sg) {  // single/single
        if ((vrfa == ALTREF_FRAME && CHECK_LAST_OR_LAST2(vrfl)) ||
            (vrfl == ALTREF_FRAME && CHECK_LAST_OR_LAST2(vrfa))) {
          pred_context = 4;
        } else if (vrfa == vrfl || (CHECK_LAST_OR_LAST2(vrfa) &&
                                    CHECK_LAST_OR_LAST2(vrfl))) {
          pred_context = 3;
        } else {  // Either vrfa or vrfl is GOLDEN / LAST3 / LAST4
          // NOTE(zoeliu): Following assert may be removed once confirmed.
          assert(CHECK_GOLDEN_LAST3_LAST4(vrfa) ||
                 CHECK_GOLDEN_LAST3_LAST4(vrfl));
          pred_context = 1;
        }
      } else if (l_sg || a_sg) {  // single/comp
        const MV_REFERENCE_FRAME vrfc = l_sg ? vrfa : vrfl;
        const MV_REFERENCE_FRAME rfs = a_sg ? vrfa : vrfl;

        if (CHECK_GOLDEN_LAST3_LAST4(vrfc) && !CHECK_GOLDEN_LAST3_LAST4(rfs))
          pred_context = 1;
        else if (CHECK_GOLDEN_LAST3_LAST4(rfs) &&
                 !CHECK_GOLDEN_LAST3_LAST4(vrfc))
          pred_context = 2;
        else
          pred_context = 4;
      } else {  // comp/comp
        if ((CHECK_LAST_OR_LAST2(vrfa) && CHECK_LAST_OR_LAST2(vrfl))) {
          pred_context = 4;
        } else {
          // NOTE(zoeliu): Following assert may be removed once confirmed.
          assert(CHECK_GOLDEN_LAST3_LAST4(vrfa) ||
                 CHECK_GOLDEN_LAST3_LAST4(vrfl));
          pred_context = 2;
        }
      }
    }
  } else if (above_in_image || left_in_image) {  // one edge available
    const MB_MODE_INFO *edge_mbmi = above_in_image ? above_mbmi : left_mbmi;

    if (!is_inter_block(edge_mbmi)) {
      pred_context = 2;
    } else {
      if (has_second_ref(edge_mbmi))
        pred_context =
            4 * (!CHECK_GOLDEN_LAST3_LAST4(edge_mbmi->ref_frame[var_ref_idx]));
      else
        pred_context = 3 * (!CHECK_GOLDEN_LAST3_LAST4(edge_mbmi->ref_frame[0]));
    }
  } else {  // no edges available (2)
    pred_context = 2;
  }

  assert(pred_context >= 0 && pred_context < REF_CONTEXTS);

  return pred_context;
}

// Returns a context number for the given MB prediction signal
// Signal the first reference frame for a compound mode is LAST,
// conditioning on that it is known either LAST/LAST2.
//
// NOTE(zoeliu): The probability of ref_frame[0] is LAST_FRAME,
// conditioning on it is either LAST_FRAME or LAST2_FRAME.
int vp10_get_pred_context_comp_ref_p1(const VP10_COMMON *cm,
                                      const MACROBLOCKD *xd) {
  int pred_context;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int above_in_image = xd->up_available;
  const int left_in_image = xd->left_available;

  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries correpsonding to real macroblocks.
  // The prediction flags in these dummy entries are initialised to 0.
  const int fix_ref_idx = cm->ref_frame_sign_bias[cm->comp_fixed_ref];
  const int var_ref_idx = !fix_ref_idx;

  if (above_in_image && left_in_image) {  // both edges available
    const int above_intra = !is_inter_block(above_mbmi);
    const int left_intra = !is_inter_block(left_mbmi);

    if (above_intra && left_intra) {  // intra/intra (2)
      pred_context = 2;
    } else if (above_intra || left_intra) {  // intra/inter
      const MB_MODE_INFO *edge_mbmi = above_intra ? left_mbmi : above_mbmi;

      if (!has_second_ref(edge_mbmi))  // single pred (1/3)
        pred_context = 1 + 2 * (edge_mbmi->ref_frame[0] != LAST_FRAME);
      else  // comp pred (1/3)
        pred_context = 1 + 2 * (edge_mbmi->ref_frame[var_ref_idx]
                                != LAST_FRAME);
    } else {  // inter/inter
      const int l_sg = !has_second_ref(left_mbmi);
      const int a_sg = !has_second_ref(above_mbmi);
      const MV_REFERENCE_FRAME vrfa = a_sg ? above_mbmi->ref_frame[0]
                                           : above_mbmi->ref_frame[var_ref_idx];
      const MV_REFERENCE_FRAME vrfl = l_sg ? left_mbmi->ref_frame[0]
                                           : left_mbmi->ref_frame[var_ref_idx];

      if (vrfa == vrfl && vrfa == LAST_FRAME)
        pred_context = 0;
      else if (l_sg && a_sg) {  // single/single
        if (vrfa == LAST_FRAME || vrfl == LAST_FRAME)
          pred_context = 1;
        else if (CHECK_GOLDEN_LAST3_LAST4(vrfa) ||
                 CHECK_GOLDEN_LAST3_LAST4(vrfl))
          pred_context = 2 + (vrfa != vrfl);
        else if (vrfa == vrfl)
          pred_context = 3;
        else
          pred_context = 4;
      } else if (l_sg || a_sg) {  // single/comp
        const MV_REFERENCE_FRAME vrfc = l_sg ? vrfa : vrfl;
        const MV_REFERENCE_FRAME rfs = a_sg ? vrfa : vrfl;

        if (vrfc == LAST_FRAME && rfs != LAST_FRAME)
          pred_context = 1;
        else if (rfs == LAST_FRAME && vrfc != LAST_FRAME)
          pred_context = 2;
        else
          pred_context = 3 +
              (vrfc == LAST2_FRAME || CHECK_GOLDEN_LAST3_LAST4(rfs));
      } else {  // comp/comp
        if (vrfa == LAST_FRAME || vrfl == LAST_FRAME)
          pred_context = 2;
        else
          pred_context = 3 + (CHECK_GOLDEN_LAST3_LAST4(vrfa) ||
                              CHECK_GOLDEN_LAST3_LAST4(vrfl));
      }
    }
  } else if (above_in_image || left_in_image) {  // one edge available
    const MB_MODE_INFO *edge_mbmi = above_in_image ? above_mbmi : left_mbmi;

    if (!is_inter_block(edge_mbmi)) {
      pred_context = 2;
    } else {
      if (has_second_ref(edge_mbmi)) {
        pred_context = 4 * (edge_mbmi->ref_frame[var_ref_idx] != LAST_FRAME);
      } else {
        if (edge_mbmi->ref_frame[0] == LAST_FRAME)
          pred_context = 0;
        else
          pred_context = 2 + CHECK_GOLDEN_LAST3_LAST4(edge_mbmi->ref_frame[0]);
      }
    }
  } else {  // no edges available (2)
    pred_context = 2;
  }

  assert(pred_context >= 0 && pred_context < REF_CONTEXTS);

  return pred_context;
}

#define CHECK_LAST3_OR_LAST4(ref_frame) \
  ((ref_frame == LAST3_FRAME) || (ref_frame == LAST4_FRAME))

// Returns a context number for the given MB prediction signal
// Signal the first reference frame for a compound mode is GOLDEN,
// conditioning on that it is known either GOLDEN/LAST3/LAST4.
//
// NOTE(zoeliu): The probability of ref_frame[0] is GOLDEN_FRAME,
// conditioning on it is either GOLDEN / LAST3 / LAST4.
int vp10_get_pred_context_comp_ref_p2(const VP10_COMMON *cm,
                                      const MACROBLOCKD *xd) {
  int pred_context;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int above_in_image = xd->up_available;
  const int left_in_image = xd->left_available;

  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries correpsonding to real macroblocks.
  // The prediction flags in these dummy entries are initialised to 0.
  const int fix_ref_idx = cm->ref_frame_sign_bias[cm->comp_fixed_ref];
  const int var_ref_idx = !fix_ref_idx;

  if (above_in_image && left_in_image) {  // both edges available
    const int above_intra = !is_inter_block(above_mbmi);
    const int left_intra = !is_inter_block(left_mbmi);

    if (above_intra && left_intra) {  // intra/intra (2)
      pred_context = 2;
    } else if (above_intra || left_intra) {  // intra/inter
      const MB_MODE_INFO *edge_mbmi = above_intra ? left_mbmi : above_mbmi;

      if (!has_second_ref(edge_mbmi))  // single pred (1/3)
        pred_context = 1 + 2 * (edge_mbmi->ref_frame[0] != GOLDEN_FRAME);
      else  // comp pred (1/3)
        pred_context = 1 +
            2 * (edge_mbmi->ref_frame[var_ref_idx] != GOLDEN_FRAME);
    } else {  // inter/inter
      const int l_sg = !has_second_ref(left_mbmi);
      const int a_sg = !has_second_ref(above_mbmi);
      const MV_REFERENCE_FRAME vrfa = a_sg ? above_mbmi->ref_frame[0]
                                           : above_mbmi->ref_frame[var_ref_idx];
      const MV_REFERENCE_FRAME vrfl = l_sg ? left_mbmi->ref_frame[0]
                                           : left_mbmi->ref_frame[var_ref_idx];

      if (vrfa == vrfl && vrfa == GOLDEN_FRAME)
        pred_context = 0;
      else if (l_sg && a_sg) {  // single/single
        if (vrfa == GOLDEN_FRAME || vrfl == GOLDEN_FRAME)
          pred_context = 1;
        else if (CHECK_LAST_OR_LAST2(vrfa) || CHECK_LAST_OR_LAST2(vrfl))
          pred_context = 2 + (vrfa != vrfl);
        else if (vrfa == vrfl)
          pred_context = 3;
        else
          pred_context = 4;
      } else if (l_sg || a_sg) {  // single/comp
        const MV_REFERENCE_FRAME vrfc = l_sg ? vrfa : vrfl;
        const MV_REFERENCE_FRAME rfs = a_sg ? vrfa : vrfl;

        if (vrfc == GOLDEN_FRAME && rfs != GOLDEN_FRAME)
          pred_context = 1;
        else if (rfs == GOLDEN_FRAME && vrfc != GOLDEN_FRAME)
          pred_context = 2;
        else
          pred_context = 3 +
              (CHECK_LAST3_OR_LAST4(vrfc) || CHECK_LAST_OR_LAST2(rfs));
      } else {  // comp/comp
        if (vrfa == GOLDEN_FRAME || vrfl == GOLDEN_FRAME)
          pred_context = 2;
        else
          pred_context = 3 +
              (CHECK_LAST_OR_LAST2(vrfa) || CHECK_LAST_OR_LAST2(vrfl));
      }
    }
  } else if (above_in_image || left_in_image) {  // one edge available
    const MB_MODE_INFO *edge_mbmi = above_in_image ? above_mbmi : left_mbmi;

    if (!is_inter_block(edge_mbmi)) {
      pred_context = 2;
    } else {
      if (has_second_ref(edge_mbmi)) {
        pred_context = 4 * (edge_mbmi->ref_frame[var_ref_idx] != GOLDEN_FRAME);
      } else {
        if (edge_mbmi->ref_frame[0] == GOLDEN_FRAME)
          pred_context = 0;
        else
          pred_context = 2 + CHECK_LAST_OR_LAST2(edge_mbmi->ref_frame[0]);
      }
    }
  } else {  // no edges available (2)
    pred_context = 2;
  }

  assert(pred_context >= 0 && pred_context < REF_CONTEXTS);

  return pred_context;
}

#define CHECK_LAST_LAST2_GOLDEN(ref_frame) \
  ((ref_frame == LAST_FRAME) || (ref_frame == LAST2_FRAME) || \
  (ref_frame == GOLDEN_FRAME))

// Returns a context number for the given MB prediction signal
// Signal the first reference frame for a compound mode is LAST3,
// conditioning on that it is known either LAST3/LAST4.
//
// NOTE(zoeliu): The probability of ref_frame[0] is LAST3_FRAME,
// conditioning on it is either LAST3 / LAST4.
int vp10_get_pred_context_comp_ref_p3(const VP10_COMMON *cm,
                                      const MACROBLOCKD *xd) {
  int pred_context;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int above_in_image = xd->up_available;
  const int left_in_image = xd->left_available;

  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries correpsonding to real macroblocks.
  // The prediction flags in these dummy entries are initialised to 0.
  const int fix_ref_idx = cm->ref_frame_sign_bias[cm->comp_fixed_ref];
  const int var_ref_idx = !fix_ref_idx;

  if (above_in_image && left_in_image) {  // both edges available
    const int above_intra = !is_inter_block(above_mbmi);
    const int left_intra = !is_inter_block(left_mbmi);

    if (above_intra && left_intra) {  // intra/intra (2)
      pred_context = 2;
    } else if (above_intra || left_intra) {  // intra/inter
      const MB_MODE_INFO *edge_mbmi = above_intra ? left_mbmi : above_mbmi;

      if (!has_second_ref(edge_mbmi))  // single pred (1/3)
        pred_context = 1 + 2 * (edge_mbmi->ref_frame[0] != LAST3_FRAME);
      else  // comp pred (1/3)
        pred_context = 1 +
            2 * (edge_mbmi->ref_frame[var_ref_idx] != LAST3_FRAME);
    } else {  // inter/inter
      const int l_sg = !has_second_ref(left_mbmi);
      const int a_sg = !has_second_ref(above_mbmi);
      const MV_REFERENCE_FRAME vrfa = a_sg ? above_mbmi->ref_frame[0]
                                           : above_mbmi->ref_frame[var_ref_idx];
      const MV_REFERENCE_FRAME vrfl = l_sg ? left_mbmi->ref_frame[0]
                                           : left_mbmi->ref_frame[var_ref_idx];

      if (vrfa == vrfl && vrfa == LAST3_FRAME)
        pred_context = 0;
      else if (l_sg && a_sg) {  // single/single
        if (vrfa == LAST3_FRAME || vrfl == LAST3_FRAME)
          pred_context = 1;
        else if (CHECK_LAST_LAST2_GOLDEN(vrfa) || CHECK_LAST_LAST2_GOLDEN(vrfl))
          pred_context = 2 + (vrfa != vrfl);
        else if (vrfa == vrfl)
          pred_context = 3;
        else
          pred_context = 4;
      } else if (l_sg || a_sg) {  // single/comp
        const MV_REFERENCE_FRAME vrfc = l_sg ? vrfa : vrfl;
        const MV_REFERENCE_FRAME rfs = a_sg ? vrfa : vrfl;

        if (vrfc == LAST3_FRAME && rfs != LAST3_FRAME)
          pred_context = 1;
        else if (rfs == LAST3_FRAME && vrfc != LAST3_FRAME)
          pred_context = 2;
        else
          pred_context = 3 +
              (vrfc == LAST4_FRAME || CHECK_LAST_LAST2_GOLDEN(rfs));
      } else {  // comp/comp
        if (vrfa == LAST3_FRAME || vrfl == LAST3_FRAME)
          pred_context = 2;
        else
          pred_context = 3 +
              (CHECK_LAST_LAST2_GOLDEN(vrfa) || CHECK_LAST_LAST2_GOLDEN(vrfl));
      }
    }
  } else if (above_in_image || left_in_image) {  // one edge available
    const MB_MODE_INFO *edge_mbmi = above_in_image ? above_mbmi : left_mbmi;

    if (!is_inter_block(edge_mbmi)) {
      pred_context = 2;
    } else {
      if (has_second_ref(edge_mbmi)) {
        pred_context = 4 * (edge_mbmi->ref_frame[var_ref_idx] != LAST3_FRAME);
      } else {
        if (edge_mbmi->ref_frame[0] == LAST3_FRAME)
          pred_context = 0;
        else
          pred_context = 2 + CHECK_LAST_LAST2_GOLDEN(edge_mbmi->ref_frame[0]);
      }
    }
  } else {  // no edges available (2)
    pred_context = 2;
  }

  assert(pred_context >= 0 && pred_context < REF_CONTEXTS);

  return pred_context;
}

#else  // CONFIG_EXT_REFS

// Returns a context number for the given MB prediction signal
int vp10_get_pred_context_comp_ref_p(const VP10_COMMON *cm,
                                     const MACROBLOCKD *xd) {
  int pred_context;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int above_in_image = xd->up_available;
  const int left_in_image = xd->left_available;

  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries corresponding to real macroblocks.
  // The prediction flags in these dummy entries are initialized to 0.
  const int fix_ref_idx = cm->ref_frame_sign_bias[cm->comp_fixed_ref];
  const int var_ref_idx = !fix_ref_idx;

  if (above_in_image && left_in_image) {  // both edges available
    const int above_intra = !is_inter_block(above_mbmi);
    const int left_intra = !is_inter_block(left_mbmi);

    if (above_intra && left_intra) {  // intra/intra (2)
      pred_context = 2;
    } else if (above_intra || left_intra) {  // intra/inter
      const MB_MODE_INFO *edge_mbmi = above_intra ? left_mbmi : above_mbmi;

      if (!has_second_ref(edge_mbmi))  // single pred (1/3)
        pred_context = 1 + 2 * (edge_mbmi->ref_frame[0] != cm->comp_var_ref[1]);
      else  // comp pred (1/3)
        pred_context = 1 + 2 * (edge_mbmi->ref_frame[var_ref_idx]
                                    != cm->comp_var_ref[1]);
    } else {  // inter/inter
      const int l_sg = !has_second_ref(left_mbmi);
      const int a_sg = !has_second_ref(above_mbmi);
      const MV_REFERENCE_FRAME vrfa = a_sg ? above_mbmi->ref_frame[0]
                                           : above_mbmi->ref_frame[var_ref_idx];
      const MV_REFERENCE_FRAME vrfl = l_sg ? left_mbmi->ref_frame[0]
                                           : left_mbmi->ref_frame[var_ref_idx];

      if (vrfa == vrfl && cm->comp_var_ref[1] == vrfa) {
        pred_context = 0;
      } else if (l_sg && a_sg) {  // single/single
        if ((vrfa == cm->comp_fixed_ref && vrfl == cm->comp_var_ref[0]) ||
            (vrfl == cm->comp_fixed_ref && vrfa == cm->comp_var_ref[0]))
          pred_context = 4;
        else if (vrfa == vrfl)
          pred_context = 3;
        else
          pred_context = 1;
      } else if (l_sg || a_sg) {  // single/comp
        const MV_REFERENCE_FRAME vrfc = l_sg ? vrfa : vrfl;
        const MV_REFERENCE_FRAME rfs = a_sg ? vrfa : vrfl;
        if (vrfc == cm->comp_var_ref[1] && rfs != cm->comp_var_ref[1])
          pred_context = 1;
        else if (rfs == cm->comp_var_ref[1] && vrfc != cm->comp_var_ref[1])
          pred_context = 2;
        else
          pred_context = 4;
      } else if (vrfa == vrfl) {  // comp/comp
        pred_context = 4;
      } else {
        pred_context = 2;
      }
    }
  } else if (above_in_image || left_in_image) {  // one edge available
    const MB_MODE_INFO *edge_mbmi = above_in_image ? above_mbmi : left_mbmi;

    if (!is_inter_block(edge_mbmi)) {
      pred_context = 2;
    } else {
      if (has_second_ref(edge_mbmi))
        pred_context = 4 * (edge_mbmi->ref_frame[var_ref_idx]
                              != cm->comp_var_ref[1]);
      else
        pred_context = 3 * (edge_mbmi->ref_frame[0] != cm->comp_var_ref[1]);
    }
  } else {  // no edges available (2)
    pred_context = 2;
  }
  assert(pred_context >= 0 && pred_context < REF_CONTEXTS);

  return pred_context;
}

#endif  // CONFIG_EXT_REFS

#if CONFIG_EXT_REFS

#define CHECK_GOLDEN_OR_ALTREF(ref_frame) \
  ((ref_frame == GOLDEN_FRAME) || (ref_frame == ALTREF_FRAME))

// For the bit to signal whether the single reference is a ALTREF_FRAME
// or a GOLDEN_FRAME.
//
// NOTE(zoeliu): The probability of ref_frame[0] is ALTREF/GOLDEN.
int vp10_get_pred_context_single_ref_p1(const MACROBLOCKD *xd) {
  int pred_context;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int has_above = xd->up_available;
  const int has_left = xd->left_available;

  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries correpsonding to real macroblocks.
  // The prediction flags in these dummy entries are initialised to 0.
  if (has_above && has_left) {  // both edges available
    const int above_intra = !is_inter_block(above_mbmi);
    const int left_intra = !is_inter_block(left_mbmi);

    if (above_intra && left_intra) {  // intra/intra
      pred_context = 2;
    } else if (above_intra || left_intra) {  // intra/inter or inter/intra
      const MB_MODE_INFO *edge_mbmi = above_intra ? left_mbmi : above_mbmi;

      if (!has_second_ref(edge_mbmi))
        pred_context = 4 * (!CHECK_GOLDEN_OR_ALTREF(edge_mbmi->ref_frame[0]));
      else
        pred_context = 1 + (!CHECK_GOLDEN_OR_ALTREF(edge_mbmi->ref_frame[0]) ||
                            !CHECK_GOLDEN_OR_ALTREF(edge_mbmi->ref_frame[1]));
    } else {  // inter/inter
      const int above_has_second = has_second_ref(above_mbmi);
      const int left_has_second  = has_second_ref(left_mbmi);

      const MV_REFERENCE_FRAME above0 = above_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME above1 = above_mbmi->ref_frame[1];
      const MV_REFERENCE_FRAME left0 = left_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME left1 = left_mbmi->ref_frame[1];

      if (above_has_second && left_has_second) {
        pred_context = 1 + (!CHECK_GOLDEN_OR_ALTREF(above0) ||
                            !CHECK_GOLDEN_OR_ALTREF(above1) ||
                            !CHECK_GOLDEN_OR_ALTREF(left0) ||
                            !CHECK_GOLDEN_OR_ALTREF(left1));
      } else if (above_has_second || left_has_second) {
        const MV_REFERENCE_FRAME rfs = !above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf1 = above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf2 = above_has_second ? above1 : left1;

        if (!CHECK_GOLDEN_OR_ALTREF(rfs))
          pred_context = 3 + (!CHECK_GOLDEN_OR_ALTREF(crf1) ||
                              !CHECK_GOLDEN_OR_ALTREF(crf2));
        else
          pred_context = !CHECK_GOLDEN_OR_ALTREF(crf1) ||
                         !CHECK_GOLDEN_OR_ALTREF(crf2);
      } else {
        pred_context = 2 * (!CHECK_GOLDEN_OR_ALTREF(above0)) +
                       2 * (!CHECK_GOLDEN_OR_ALTREF(left0));
      }
    }
  } else if (has_above || has_left) {  // one edge available
    const MB_MODE_INFO *edge_mbmi = has_above ? above_mbmi : left_mbmi;
    if (!is_inter_block(edge_mbmi)) {  // intra
      pred_context = 2;
    } else {  // inter
      if (!has_second_ref(edge_mbmi))
        pred_context = 4 * (!CHECK_GOLDEN_OR_ALTREF(edge_mbmi->ref_frame[0]));
      else
        pred_context = 1 + (!CHECK_GOLDEN_OR_ALTREF(edge_mbmi->ref_frame[0]) ||
                            !CHECK_GOLDEN_OR_ALTREF(edge_mbmi->ref_frame[1]));
    }
  } else {  // no edges available
    pred_context = 2;
  }

  assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
  return pred_context;
}

// For the bit to signal whether the single reference is ALTREF_FRAME or
// GOLDEN_FRAME, knowing that it shall be either of these 2 choices.
//
// NOTE(zoeliu): The probability of ref_frame[0] is ALTREF_FRAME, conditioning
// on it is either ALTREF_FRAME/GOLDEN_FRAME.
int vp10_get_pred_context_single_ref_p2(const MACROBLOCKD *xd) {
  int pred_context;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int has_above = xd->up_available;
  const int has_left = xd->left_available;

  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries correpsonding to real macroblocks.
  // The prediction flags in these dummy entries are initialised to 0.
  if (has_above && has_left) {  // both edges available
    const int above_intra = !is_inter_block(above_mbmi);
    const int left_intra = !is_inter_block(left_mbmi);

    if (above_intra && left_intra) {  // intra/intra
      pred_context = 2;
    } else if (above_intra || left_intra) {  // intra/inter or inter/intra
      const MB_MODE_INFO *edge_mbmi = above_intra ? left_mbmi : above_mbmi;
      if (!has_second_ref(edge_mbmi)) {
        if (!CHECK_GOLDEN_OR_ALTREF(edge_mbmi->ref_frame[0]))
          pred_context = 3;
        else
          pred_context = 4 * (edge_mbmi->ref_frame[0] == GOLDEN_FRAME);
      } else {
        pred_context = 1 + 2 * (edge_mbmi->ref_frame[0] == GOLDEN_FRAME ||
                                edge_mbmi->ref_frame[1] == GOLDEN_FRAME);
      }
    } else {  // inter/inter
      const int above_has_second = has_second_ref(above_mbmi);
      const int left_has_second  = has_second_ref(left_mbmi);
      const MV_REFERENCE_FRAME above0 = above_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME above1 = above_mbmi->ref_frame[1];
      const MV_REFERENCE_FRAME left0 = left_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME left1 = left_mbmi->ref_frame[1];

      if (above_has_second && left_has_second) {
        if (above0 == left0 && above1 == left1)
          pred_context = 3 * (above0 == GOLDEN_FRAME ||
                              above1 == GOLDEN_FRAME ||
                              left0 == GOLDEN_FRAME ||
                              left1 == GOLDEN_FRAME);
        else
          pred_context = 2;
      } else if (above_has_second || left_has_second) {
        const MV_REFERENCE_FRAME rfs = !above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf1 = above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf2 = above_has_second ? above1 : left1;

        if (rfs == GOLDEN_FRAME)
          pred_context = 3 + (crf1 == GOLDEN_FRAME || crf2 == GOLDEN_FRAME);
        else if (rfs == ALTREF_FRAME)
          pred_context = (crf1 == GOLDEN_FRAME || crf2 == GOLDEN_FRAME);
        else
          pred_context = 1 + 2 * (crf1 == GOLDEN_FRAME || crf2 == GOLDEN_FRAME);
      } else {
        if (!CHECK_GOLDEN_OR_ALTREF(above0) && !CHECK_GOLDEN_OR_ALTREF(left0)) {
          pred_context = 2 + (above0 == left0);
        } else if (!CHECK_GOLDEN_OR_ALTREF(above0) ||
                   !CHECK_GOLDEN_OR_ALTREF(left0)) {
          const MV_REFERENCE_FRAME edge0 =
              !CHECK_GOLDEN_OR_ALTREF(above0) ? left0 : above0;
          pred_context = 4 * (edge0 == GOLDEN_FRAME);
        } else {
          pred_context = 2 * (above0 == GOLDEN_FRAME) +
                         2 * (left0  == GOLDEN_FRAME);
        }
      }
    }
  } else if (has_above || has_left) {  // one edge available
    const MB_MODE_INFO *edge_mbmi = has_above ? above_mbmi : left_mbmi;

    if (!is_inter_block(edge_mbmi) ||
        (!CHECK_GOLDEN_OR_ALTREF(edge_mbmi->ref_frame[0]) &&
         !has_second_ref(edge_mbmi)))
      pred_context = 2;
    else if (!has_second_ref(edge_mbmi))
      pred_context = 4 * (edge_mbmi->ref_frame[0] == GOLDEN_FRAME);
    else
      pred_context = 3 * (edge_mbmi->ref_frame[0] == GOLDEN_FRAME ||
                          edge_mbmi->ref_frame[1] == GOLDEN_FRAME);
  } else {  // no edges available (2)
    pred_context = 2;
  }

  assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
  return pred_context;
}

// For the bit to signal whether the single reference is LAST3/LAST4 or
// LAST2/LAST, knowing that it shall be either of these 2 choices.
//
// NOTE(zoeliu): The probability of ref_frame[0] is LAST3/LAST4, conditioning
// on it is either LAST3/LAST4/LAST2/LAST.
int vp10_get_pred_context_single_ref_p3(const MACROBLOCKD *xd) {
  int pred_context;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int has_above = xd->up_available;
  const int has_left = xd->left_available;

  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries correpsonding to real macroblocks.
  // The prediction flags in these dummy entries are initialised to 0.
  if (has_above && has_left) {  // both edges available
    const int above_intra = !is_inter_block(above_mbmi);
    const int left_intra = !is_inter_block(left_mbmi);

    if (above_intra && left_intra) {  // intra/intra
      pred_context = 2;
    } else if (above_intra || left_intra) {  // intra/inter or inter/intra
      const MB_MODE_INFO *edge_mbmi = above_intra ? left_mbmi : above_mbmi;
      if (!has_second_ref(edge_mbmi)) {
        if (CHECK_GOLDEN_OR_ALTREF(edge_mbmi->ref_frame[0]))
          pred_context = 3;
        else
          pred_context = 4 * CHECK_LAST_OR_LAST2(edge_mbmi->ref_frame[0]);
      } else {
        pred_context = 1 +
            2 * (CHECK_LAST_OR_LAST2(edge_mbmi->ref_frame[0]) ||
                 CHECK_LAST_OR_LAST2(edge_mbmi->ref_frame[1]));
      }
    } else {  // inter/inter
      const int above_has_second = has_second_ref(above_mbmi);
      const int left_has_second  = has_second_ref(left_mbmi);
      const MV_REFERENCE_FRAME above0 = above_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME above1 = above_mbmi->ref_frame[1];
      const MV_REFERENCE_FRAME left0 = left_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME left1 = left_mbmi->ref_frame[1];

      if (above_has_second && left_has_second) {
        if (above0 == left0 && above1 == left1)
          pred_context = 3 * (CHECK_LAST_OR_LAST2(above0) ||
                              CHECK_LAST_OR_LAST2(above1) ||
                              CHECK_LAST_OR_LAST2(left0) ||
                              CHECK_LAST_OR_LAST2(left1));
        else
          pred_context = 2;
      } else if (above_has_second || left_has_second) {
        const MV_REFERENCE_FRAME rfs = !above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf1 = above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf2 = above_has_second ? above1 : left1;

        if (CHECK_LAST_OR_LAST2(rfs))
          pred_context = 3 + (CHECK_LAST_OR_LAST2(crf1) ||
                              CHECK_LAST_OR_LAST2(crf2));
        else if (rfs == LAST3_FRAME || rfs == LAST4_FRAME)
          pred_context = (CHECK_LAST_OR_LAST2(crf1) ||
                          CHECK_LAST_OR_LAST2(crf2));
        else
          pred_context = 1 + 2 * (CHECK_LAST_OR_LAST2(crf1) ||
                                  CHECK_LAST_OR_LAST2(crf2));
      } else {
        if (CHECK_GOLDEN_OR_ALTREF(above0) && CHECK_GOLDEN_OR_ALTREF(left0)) {
          pred_context = 2 + (above0 == left0);
        } else if (CHECK_GOLDEN_OR_ALTREF(above0) ||
                   CHECK_GOLDEN_OR_ALTREF(left0)) {
          const MV_REFERENCE_FRAME edge0 =
              CHECK_GOLDEN_OR_ALTREF(above0) ? left0 : above0;
          pred_context = 4 * CHECK_LAST_OR_LAST2(edge0);
        } else {
          pred_context = 2 * CHECK_LAST_OR_LAST2(above0) +
                         2 * CHECK_LAST_OR_LAST2(left0);
        }
      }
    }
  } else if (has_above || has_left) {  // one edge available
    const MB_MODE_INFO *edge_mbmi = has_above ? above_mbmi : left_mbmi;

    if (!is_inter_block(edge_mbmi) ||
        (CHECK_GOLDEN_OR_ALTREF(edge_mbmi->ref_frame[0]) &&
         !has_second_ref(edge_mbmi)))
      pred_context = 2;
    else if (!has_second_ref(edge_mbmi))
      pred_context = 4 * (CHECK_LAST_OR_LAST2(edge_mbmi->ref_frame[0]));
    else
      pred_context = 3 * (CHECK_LAST_OR_LAST2(edge_mbmi->ref_frame[0]) ||
                          CHECK_LAST_OR_LAST2(edge_mbmi->ref_frame[1]));
  } else {  // no edges available (2)
    pred_context = 2;
  }

  assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
  return pred_context;
}

// For the bit to signal whether the single reference is LAST2_FRAME or
// LAST_FRAME, knowing that it shall be either of these 2 choices.
//
// NOTE(zoeliu): The probability of ref_frame[0] is LAST2_FRAME, conditioning
// on it is either LAST2_FRAME/LAST_FRAME.
int vp10_get_pred_context_single_ref_p4(const MACROBLOCKD *xd) {
  int pred_context;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int has_above = xd->up_available;
  const int has_left = xd->left_available;

  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries correpsonding to real macroblocks.
  // The prediction flags in these dummy entries are initialised to 0.
  if (has_above && has_left) {  // both edges available
    const int above_intra = !is_inter_block(above_mbmi);
    const int left_intra = !is_inter_block(left_mbmi);

    if (above_intra && left_intra) {  // intra/intra
      pred_context = 2;
    } else if (above_intra || left_intra) {  // intra/inter or inter/intra
      const MB_MODE_INFO *edge_mbmi = above_intra ? left_mbmi : above_mbmi;
      if (!has_second_ref(edge_mbmi)) {
        if (!CHECK_LAST_OR_LAST2(edge_mbmi->ref_frame[0]))
          pred_context = 3;
        else
          pred_context = 4 * (edge_mbmi->ref_frame[0] == LAST_FRAME);
      } else {
        pred_context = 1 +
            2 * (edge_mbmi->ref_frame[0] == LAST_FRAME ||
                 edge_mbmi->ref_frame[1] == LAST_FRAME);
      }
    } else {  // inter/inter
      const int above_has_second = has_second_ref(above_mbmi);
      const int left_has_second  = has_second_ref(left_mbmi);
      const MV_REFERENCE_FRAME above0 = above_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME above1 = above_mbmi->ref_frame[1];
      const MV_REFERENCE_FRAME left0 = left_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME left1 = left_mbmi->ref_frame[1];

      if (above_has_second && left_has_second) {
        if (above0 == left0 && above1 == left1)
          pred_context = 3 * (above0 == LAST_FRAME || above1 == LAST_FRAME ||
                              left0 == LAST_FRAME || left1 == LAST_FRAME);
        else
          pred_context = 2;
      } else if (above_has_second || left_has_second) {
        const MV_REFERENCE_FRAME rfs = !above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf1 = above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf2 = above_has_second ? above1 : left1;

        if (rfs == LAST_FRAME)
          pred_context = 3 + (crf1 == LAST_FRAME || crf2 == LAST_FRAME);
        else if (rfs == LAST2_FRAME)
          pred_context = (crf1 == LAST_FRAME || crf2 == LAST_FRAME);
        else
          pred_context = 1 + 2 * (crf1 == LAST_FRAME || crf2 == LAST_FRAME);
      } else {
        if (!CHECK_LAST_OR_LAST2(above0) &&
            !CHECK_LAST_OR_LAST2(left0)) {
          pred_context = 2 + (above0 == left0);
        } else if (!CHECK_LAST_OR_LAST2(above0) ||
                   !CHECK_LAST_OR_LAST2(left0)) {
          const MV_REFERENCE_FRAME edge0 =
              !CHECK_LAST_OR_LAST2(above0) ? left0 : above0;
          pred_context = 4 * (edge0 == LAST_FRAME);
        } else {
          pred_context = 2 * (above0 == LAST_FRAME) + 2 * (left0 == LAST_FRAME);
        }
      }
    }
  } else if (has_above || has_left) {  // one edge available
    const MB_MODE_INFO *edge_mbmi = has_above ? above_mbmi : left_mbmi;

    if (!is_inter_block(edge_mbmi) ||
        (!CHECK_LAST_OR_LAST2(edge_mbmi->ref_frame[0]) &&
         !has_second_ref(edge_mbmi)))
      pred_context = 2;
    else if (!has_second_ref(edge_mbmi))
      pred_context = 4 * (edge_mbmi->ref_frame[0] == LAST_FRAME);
    else
      pred_context = 3 * (edge_mbmi->ref_frame[0] == LAST_FRAME ||
                          edge_mbmi->ref_frame[1] == LAST_FRAME);
  } else {  // no edges available (2)
    pred_context = 2;
  }

  assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
  return pred_context;
}

// For the bit to signal whether the single reference is LAST4_FRAME or
// LAST3_FRAME, knowing that it shall be either of these 2 choices.
//
// NOTE(zoeliu): The probability of ref_frame[0] is LAST4_FRAME, conditioning
// on it is either LAST4_FRAME/LAST3_FRAME.
int vp10_get_pred_context_single_ref_p5(const MACROBLOCKD *xd) {
  int pred_context;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int has_above = xd->up_available;
  const int has_left = xd->left_available;

  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries correpsonding to real macroblocks.
  // The prediction flags in these dummy entries are initialised to 0.
  if (has_above && has_left) {  // both edges available
    const int above_intra = !is_inter_block(above_mbmi);
    const int left_intra = !is_inter_block(left_mbmi);

    if (above_intra && left_intra) {  // intra/intra
      pred_context = 2;
    } else if (above_intra || left_intra) {  // intra/inter or inter/intra
      const MB_MODE_INFO *edge_mbmi = above_intra ? left_mbmi : above_mbmi;
      if (!has_second_ref(edge_mbmi)) {
        if (!CHECK_LAST3_OR_LAST4(edge_mbmi->ref_frame[0]))
          pred_context = 3;
        else
          pred_context = 4 * (edge_mbmi->ref_frame[0] == LAST3_FRAME);
      } else {
        pred_context = 1 +
            2 * (edge_mbmi->ref_frame[0] == LAST3_FRAME ||
                 edge_mbmi->ref_frame[1] == LAST3_FRAME);
      }
    } else {  // inter/inter
      const int above_has_second = has_second_ref(above_mbmi);
      const int left_has_second  = has_second_ref(left_mbmi);
      const MV_REFERENCE_FRAME above0 = above_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME above1 = above_mbmi->ref_frame[1];
      const MV_REFERENCE_FRAME left0 = left_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME left1 = left_mbmi->ref_frame[1];

      if (above_has_second && left_has_second) {
        if (above0 == left0 && above1 == left1)
          pred_context = 3 * (above0 == LAST3_FRAME || above1 == LAST3_FRAME ||
                              left0 == LAST3_FRAME || left1 == LAST3_FRAME);
        else
          pred_context = 2;
      } else if (above_has_second || left_has_second) {
        const MV_REFERENCE_FRAME rfs = !above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf1 = above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf2 = above_has_second ? above1 : left1;

        if (rfs == LAST3_FRAME)
          pred_context = 3 + (crf1 == LAST3_FRAME || crf2 == LAST3_FRAME);
        else if (rfs == LAST4_FRAME)
          pred_context = (crf1 == LAST3_FRAME || crf2 == LAST3_FRAME);
        else
          pred_context = 1 + 2 * (crf1 == LAST3_FRAME || crf2 == LAST3_FRAME);
      } else {
        if (!CHECK_LAST3_OR_LAST4(above0) &&
            !CHECK_LAST3_OR_LAST4(left0)) {
          pred_context = 2 + (above0 == left0);
        } else if (!CHECK_LAST3_OR_LAST4(above0) ||
                   !CHECK_LAST3_OR_LAST4(left0)) {
          const MV_REFERENCE_FRAME edge0 =
              !CHECK_LAST3_OR_LAST4(above0) ? left0 : above0;
          pred_context = 4 * (edge0 == LAST3_FRAME);
        } else {
          pred_context = 2 * (above0 == LAST3_FRAME) +
                         2 * (left0 == LAST3_FRAME);
        }
      }
    }
  } else if (has_above || has_left) {  // one edge available
    const MB_MODE_INFO *edge_mbmi = has_above ? above_mbmi : left_mbmi;

    if (!is_inter_block(edge_mbmi) ||
        (!CHECK_LAST3_OR_LAST4(edge_mbmi->ref_frame[0]) &&
         !has_second_ref(edge_mbmi)))
      pred_context = 2;
    else if (!has_second_ref(edge_mbmi))
      pred_context = 4 * (edge_mbmi->ref_frame[0] == LAST3_FRAME);
    else
      pred_context = 3 * (edge_mbmi->ref_frame[0] == LAST3_FRAME ||
                          edge_mbmi->ref_frame[1] == LAST3_FRAME);
  } else {  // no edges available (2)
    pred_context = 2;
  }

  assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
  return pred_context;
}

#else  // CONFIG_EXT_REFS

int vp10_get_pred_context_single_ref_p1(const MACROBLOCKD *xd) {
  int pred_context;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int has_above = xd->up_available;
  const int has_left = xd->left_available;
  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries corresponding to real macroblocks.
  // The prediction flags in these dummy entries are initialized to 0.
  if (has_above && has_left) {  // both edges available
    const int above_intra = !is_inter_block(above_mbmi);
    const int left_intra = !is_inter_block(left_mbmi);

    if (above_intra && left_intra) {  // intra/intra
      pred_context = 2;
    } else if (above_intra || left_intra) {  // intra/inter or inter/intra
      const MB_MODE_INFO *edge_mbmi = above_intra ? left_mbmi : above_mbmi;
      if (!has_second_ref(edge_mbmi))
        pred_context = 4 * (edge_mbmi->ref_frame[0] == LAST_FRAME);
      else
        pred_context = 1 + (edge_mbmi->ref_frame[0] == LAST_FRAME ||
                            edge_mbmi->ref_frame[1] == LAST_FRAME);
    } else {  // inter/inter
      const int above_has_second = has_second_ref(above_mbmi);
      const int left_has_second = has_second_ref(left_mbmi);
      const MV_REFERENCE_FRAME above0 = above_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME above1 = above_mbmi->ref_frame[1];
      const MV_REFERENCE_FRAME left0 = left_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME left1 = left_mbmi->ref_frame[1];

      if (above_has_second && left_has_second) {
        pred_context = 1 + (above0 == LAST_FRAME || above1 == LAST_FRAME ||
                            left0 == LAST_FRAME || left1 == LAST_FRAME);
      } else if (above_has_second || left_has_second) {
        const MV_REFERENCE_FRAME rfs = !above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf1 = above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf2 = above_has_second ? above1 : left1;

        if (rfs == LAST_FRAME)
          pred_context = 3 + (crf1 == LAST_FRAME || crf2 == LAST_FRAME);
        else
          pred_context = (crf1 == LAST_FRAME || crf2 == LAST_FRAME);
      } else {
        pred_context = 2 * (above0 == LAST_FRAME) + 2 * (left0 == LAST_FRAME);
      }
    }
  } else if (has_above || has_left) {  // one edge available
    const MB_MODE_INFO *edge_mbmi = has_above ? above_mbmi : left_mbmi;
    if (!is_inter_block(edge_mbmi)) {  // intra
      pred_context = 2;
    } else {  // inter
      if (!has_second_ref(edge_mbmi))
        pred_context = 4 * (edge_mbmi->ref_frame[0] == LAST_FRAME);
      else
        pred_context = 1 + (edge_mbmi->ref_frame[0] == LAST_FRAME ||
                            edge_mbmi->ref_frame[1] == LAST_FRAME);
    }
  } else {  // no edges available
    pred_context = 2;
  }

  assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
  return pred_context;
}

int vp10_get_pred_context_single_ref_p2(const MACROBLOCKD *xd) {
  int pred_context;
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int has_above = xd->up_available;
  const int has_left = xd->left_available;

  // Note:
  // The mode info data structure has a one element border above and to the
  // left of the entries corresponding to real macroblocks.
  // The prediction flags in these dummy entries are initialized to 0.
  if (has_above && has_left) {  // both edges available
    const int above_intra = !is_inter_block(above_mbmi);
    const int left_intra = !is_inter_block(left_mbmi);

    if (above_intra && left_intra) {  // intra/intra
      pred_context = 2;
    } else if (above_intra || left_intra) {  // intra/inter or inter/intra
      const MB_MODE_INFO *edge_mbmi = above_intra ? left_mbmi : above_mbmi;
      if (!has_second_ref(edge_mbmi)) {
        if (edge_mbmi->ref_frame[0] == LAST_FRAME)
          pred_context = 3;
        else
          pred_context = 4 * (edge_mbmi->ref_frame[0] == GOLDEN_FRAME);
      } else {
        pred_context = 1 + 2 * (edge_mbmi->ref_frame[0] == GOLDEN_FRAME ||
                                edge_mbmi->ref_frame[1] == GOLDEN_FRAME);
      }
    } else {  // inter/inter
      const int above_has_second = has_second_ref(above_mbmi);
      const int left_has_second = has_second_ref(left_mbmi);
      const MV_REFERENCE_FRAME above0 = above_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME above1 = above_mbmi->ref_frame[1];
      const MV_REFERENCE_FRAME left0 = left_mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME left1 = left_mbmi->ref_frame[1];

      if (above_has_second && left_has_second) {
        if (above0 == left0 && above1 == left1)
          pred_context = 3 * (above0 == GOLDEN_FRAME ||
                              above1 == GOLDEN_FRAME ||
                              left0 == GOLDEN_FRAME ||
                              left1 == GOLDEN_FRAME);
        else
          pred_context = 2;
      } else if (above_has_second || left_has_second) {
        const MV_REFERENCE_FRAME rfs = !above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf1 = above_has_second ? above0 : left0;
        const MV_REFERENCE_FRAME crf2 = above_has_second ? above1 : left1;

        if (rfs == GOLDEN_FRAME)
          pred_context = 3 + (crf1 == GOLDEN_FRAME || crf2 == GOLDEN_FRAME);
        else if (rfs == ALTREF_FRAME)
          pred_context = crf1 == GOLDEN_FRAME || crf2 == GOLDEN_FRAME;
        else
          pred_context = 1 + 2 * (crf1 == GOLDEN_FRAME || crf2 == GOLDEN_FRAME);
      } else {
        if (above0 == LAST_FRAME && left0 == LAST_FRAME) {
          pred_context = 3;
        } else if (above0 == LAST_FRAME || left0 == LAST_FRAME) {
          const MV_REFERENCE_FRAME edge0 = (above0 == LAST_FRAME) ? left0
                                                                  : above0;
          pred_context = 4 * (edge0 == GOLDEN_FRAME);
        } else {
          pred_context = 2 * (above0 == GOLDEN_FRAME) +
                             2 * (left0 == GOLDEN_FRAME);
        }
      }
    }
  } else if (has_above || has_left) {  // one edge available
    const MB_MODE_INFO *edge_mbmi = has_above ? above_mbmi : left_mbmi;

    if (!is_inter_block(edge_mbmi) ||
        (edge_mbmi->ref_frame[0] == LAST_FRAME && !has_second_ref(edge_mbmi)))
      pred_context = 2;
    else if (!has_second_ref(edge_mbmi))
      pred_context = 4 * (edge_mbmi->ref_frame[0] == GOLDEN_FRAME);
    else
      pred_context = 3 * (edge_mbmi->ref_frame[0] == GOLDEN_FRAME ||
                          edge_mbmi->ref_frame[1] == GOLDEN_FRAME);
  } else {  // no edges available (2)
    pred_context = 2;
  }
  assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
  return pred_context;
}

#endif  // CONFIG_EXT_REFS
