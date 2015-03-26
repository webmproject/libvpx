/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <stdlib.h>  // qsort()

#include "./vp9_rtcd.h"
#include "./vpx_scale_rtcd.h"

#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem_ops.h"
#include "vpx_scale/vpx_scale.h"

#include "vp9/common/vp9_alloccommon.h"
#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_idct.h"
#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_quant_common.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_thread.h"
#include "vp9/common/vp9_tile_common.h"

#include "vp9/decoder/vp9_decodeframe.h"
#include "vp9/decoder/vp9_detokenize.h"
#include "vp9/decoder/vp9_decodemv.h"
#include "vp9/decoder/vp9_decoder.h"
#include "vp9/decoder/vp9_dsubexp.h"
#include "vp9/decoder/vp9_dthread.h"
#include "vp9/decoder/vp9_read_bit_buffer.h"
#include "vp9/decoder/vp9_reader.h"

#define MAX_VP9_HEADER_SIZE 80

static int is_compound_reference_allowed(const VP9_COMMON *cm) {
  int i;
  for (i = 1; i < REFS_PER_FRAME; ++i)
    if (cm->ref_frame_sign_bias[i + 1] != cm->ref_frame_sign_bias[1])
      return 1;

  return 0;
}

static void setup_compound_reference_mode(VP9_COMMON *cm) {
  if (cm->ref_frame_sign_bias[LAST_FRAME] ==
          cm->ref_frame_sign_bias[GOLDEN_FRAME]) {
    cm->comp_fixed_ref = ALTREF_FRAME;
    cm->comp_var_ref[0] = LAST_FRAME;
    cm->comp_var_ref[1] = GOLDEN_FRAME;
  } else if (cm->ref_frame_sign_bias[LAST_FRAME] ==
                 cm->ref_frame_sign_bias[ALTREF_FRAME]) {
    cm->comp_fixed_ref = GOLDEN_FRAME;
    cm->comp_var_ref[0] = LAST_FRAME;
    cm->comp_var_ref[1] = ALTREF_FRAME;
  } else {
    cm->comp_fixed_ref = LAST_FRAME;
    cm->comp_var_ref[0] = GOLDEN_FRAME;
    cm->comp_var_ref[1] = ALTREF_FRAME;
  }
}

static int read_is_valid(const uint8_t *start, size_t len, const uint8_t *end) {
  return len != 0 && len <= (size_t)(end - start);
}

static int decode_unsigned_max(struct vp9_read_bit_buffer *rb, int max) {
  const int data = vp9_rb_read_literal(rb, get_unsigned_bits(max));
  return data > max ? max : data;
}

static TX_MODE read_tx_mode(vp9_reader *r) {
  TX_MODE tx_mode = vp9_read_literal(r, 2);
#if CONFIG_TX64X64
  if (tx_mode == 2)
    tx_mode += vp9_read_bit(r);      // ALLOW_16X16 and ALLOW_32X32
  else if (tx_mode == 3)
    tx_mode += 1 + vp9_read_bit(r);  // ALLOW_64X64 and TX_MODE_SELECT
#else
  if (tx_mode == ALLOW_32X32)
    tx_mode += vp9_read_bit(r);
#endif
  return tx_mode;
}

static void read_tx_mode_probs(struct tx_probs *tx_probs, vp9_reader *r) {
  int i, j;

  for (i = 0; i < TX_SIZE_CONTEXTS; ++i)
    for (j = 0; j < 1; ++j)
      vp9_diff_update_prob(r, &tx_probs->p8x8[i][j]);

  for (i = 0; i < TX_SIZE_CONTEXTS; ++i)
    for (j = 0; j < 2; ++j)
      vp9_diff_update_prob(r, &tx_probs->p16x16[i][j]);

  for (i = 0; i < TX_SIZE_CONTEXTS; ++i)
    for (j = 0; j < 3; ++j)
      vp9_diff_update_prob(r, &tx_probs->p32x32[i][j]);

#if CONFIG_TX64X64
  for (i = 0; i < TX_SIZE_CONTEXTS; ++i)
    for (j = 0; j < 4; ++j)
      vp9_diff_update_prob(r, &tx_probs->p64x64[i][j]);
#endif
}

static void read_switchable_interp_probs(FRAME_CONTEXT *fc, vp9_reader *r) {
  int i, j;
  for (j = 0; j < SWITCHABLE_FILTER_CONTEXTS; ++j)
    for (i = 0; i < SWITCHABLE_FILTERS - 1; ++i)
      vp9_diff_update_prob(r, &fc->switchable_interp_prob[j][i]);
}

static void read_inter_mode_probs(FRAME_CONTEXT *fc, vp9_reader *r) {
  int i, j;
  for (i = 0; i < INTER_MODE_CONTEXTS; ++i)
    for (j = 0; j < INTER_MODES - 1; ++j)
      vp9_diff_update_prob(r, &fc->inter_mode_probs[i][j]);
}

static REFERENCE_MODE read_frame_reference_mode(const VP9_COMMON *cm,
                                                vp9_reader *r) {
  if (is_compound_reference_allowed(cm)) {
    return vp9_read_bit(r) ? (vp9_read_bit(r) ? REFERENCE_MODE_SELECT
                                              : COMPOUND_REFERENCE)
                           : SINGLE_REFERENCE;
  } else {
    return SINGLE_REFERENCE;
  }
}

static void read_frame_reference_mode_probs(VP9_COMMON *cm, vp9_reader *r) {
  FRAME_CONTEXT *const fc = &cm->fc;
  int i;

  if (cm->reference_mode == REFERENCE_MODE_SELECT)
    for (i = 0; i < COMP_INTER_CONTEXTS; ++i)
      vp9_diff_update_prob(r, &fc->comp_inter_prob[i]);

  if (cm->reference_mode != COMPOUND_REFERENCE)
    for (i = 0; i < REF_CONTEXTS; ++i) {
      vp9_diff_update_prob(r, &fc->single_ref_prob[i][0]);
      vp9_diff_update_prob(r, &fc->single_ref_prob[i][1]);
    }

  if (cm->reference_mode != SINGLE_REFERENCE)
    for (i = 0; i < REF_CONTEXTS; ++i)
      vp9_diff_update_prob(r, &fc->comp_ref_prob[i]);
}

static void update_mv_probs(vp9_prob *p, int n, vp9_reader *r) {
  int i;
  for (i = 0; i < n; ++i)
    if (vp9_read(r, MV_UPDATE_PROB))
      p[i] = (vp9_read_literal(r, 7) << 1) | 1;
}

static void read_mv_probs(nmv_context *ctx, int allow_hp, vp9_reader *r) {
  int i, j;

  update_mv_probs(ctx->joints, MV_JOINTS - 1, r);

  for (i = 0; i < 2; ++i) {
    nmv_component *const comp_ctx = &ctx->comps[i];
    update_mv_probs(&comp_ctx->sign, 1, r);
    update_mv_probs(comp_ctx->classes, MV_CLASSES - 1, r);
    update_mv_probs(comp_ctx->class0, CLASS0_SIZE - 1, r);
    update_mv_probs(comp_ctx->bits, MV_OFFSET_BITS, r);
  }

  for (i = 0; i < 2; ++i) {
    nmv_component *const comp_ctx = &ctx->comps[i];
    for (j = 0; j < CLASS0_SIZE; ++j)
      update_mv_probs(comp_ctx->class0_fp[j], MV_FP_SIZE - 1, r);
    update_mv_probs(comp_ctx->fp, 3, r);
  }

  if (allow_hp) {
    for (i = 0; i < 2; ++i) {
      nmv_component *const comp_ctx = &ctx->comps[i];
      update_mv_probs(&comp_ctx->class0_hp, 1, r);
      update_mv_probs(&comp_ctx->hp, 1, r);
    }
  }
}

static void setup_plane_dequants(VP9_COMMON *cm, MACROBLOCKD *xd, int q_index) {
  int i;
  xd->plane[0].dequant = cm->y_dequant[q_index];
#if CONFIG_NEW_QUANT
  xd->plane[0].dequant_val_nuq =
      (const dequant_val_type_nuq *)cm->y_dequant_val_nuq[q_index];
#endif  // CONFIG_NEW_QUANT

  for (i = 1; i < MAX_MB_PLANE; i++) {
    xd->plane[i].dequant = cm->uv_dequant[q_index];
#if CONFIG_NEW_QUANT
    xd->plane[i].dequant_val_nuq =
        (const dequant_val_type_nuq *)cm->uv_dequant_val_nuq[q_index];
#endif  // CONFIG_NEW_QUANT
  }
}

#if CONFIG_TX_SKIP
static void vp9_intra_dpcm_add(tran_low_t *dqcoeff, uint8_t *dst, int stride,
                               PREDICTION_MODE mode, int bs, int shift) {
  int r, c, temp;

  switch (mode) {
    case H_PRED:
      for (r = 0; r < bs; r++) {
        temp = dst[r * stride] + (dqcoeff[r * bs] >> shift);
        dst[r * stride] = clip_pixel(temp);
      }
      for (r = 0; r < bs; r++)
        for (c = 1; c < bs; c++) {
          temp = dst[r * stride + c - 1] +
              (dqcoeff[r * bs + c] >> shift);
          dst[r * stride + c] = clip_pixel(temp);
        }
      break;
    case V_PRED:
      for (c = 0; c < bs; c++) {
        temp = dst[c] + (dqcoeff[c] >> shift);
        dst[c] = clip_pixel(temp);
      }
      for (r = 1; r < bs; r++)
        for (c = 0; c < bs; c++) {
          temp = dst[(r - 1) * stride + c] +
              (dqcoeff[r * bs + c] >> shift);
          dst[r * stride + c] = clip_pixel(temp);
        }
      break;
    case TM_PRED:
      for (c = 0; c < bs; c++) {
        temp = dst[c] + (dqcoeff[c] >> shift);
        dst[c] = clip_pixel(temp);
      }
      for (r = 1; r < bs; r++) {
        temp = dst[r * stride] + (dqcoeff[r * bs] >> shift);
        dst[r * stride] = clip_pixel(temp);
      }
      for (r = 1; r < bs; r++)
        for (c = 1; c < bs; c++) {
          temp = dst[stride * r + c - 1] + dst[stride * (r - 1) + c] -
                 dst[stride * (r - 1) + c - 1];
          temp = clip_pixel(temp);
          temp = temp + (dqcoeff[r * bs + c] >> shift);
          dst[stride * r + c] = clip_pixel(temp);
        }
      break;
    default:
      break;
  }
}

static void vp9_intra_dpcm_add_nocoeff(uint8_t *dst, int stride,
                                       PREDICTION_MODE mode, int bs) {
  int r, c, temp;

  switch (mode) {
    case H_PRED:
      for (r = 0; r < bs; r++)
        vpx_memset(dst + r * stride + 1, dst[r * stride], bs - 1);
      break;
    case V_PRED:
      for (r = 1; r < bs; r++)
        vpx_memcpy(dst + r * stride, dst, bs * sizeof(*dst));
      break;
    case TM_PRED:
      for (r = 1; r < bs; r++)
        for (c = 1; c < bs; c++) {
          temp = dst[stride * r + c - 1] + dst[stride * (r - 1) + c] -
              dst[stride * (r - 1) + c - 1];
          dst[stride * r + c] = clip_pixel(temp);
        }
      break;
    default:
      break;
  }
}
#endif  // CONFIG_TX_SKIP

static void inverse_transform_block(MACROBLOCKD* xd, int plane, int block,
                                    TX_SIZE tx_size, uint8_t *dst, int stride,
                                    int eob) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
#if CONFIG_TX_SKIP
  MB_MODE_INFO *mbmi = &xd->mi[0].src_mi->mbmi;
  int shift = mbmi->tx_skip_shift;
  PREDICTION_MODE mode = (plane == 0) ? get_y_mode(xd->mi[0].src_mi, block):
                                        mbmi->uv_mode;
  (void) mode;
#endif
  if (eob > 0) {
    TX_TYPE tx_type = DCT_DCT;
    tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
#if CONFIG_VP9_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
#if CONFIG_TX_SKIP
      if (xd->lossless && !mbmi->tx_skip[plane != 0]) {
#else
      if (xd->lossless) {
#endif  // CONFIG_TX_SKIP
        tx_type = DCT_DCT;
        vp9_highbd_iwht4x4_add(dqcoeff, dst, stride, eob, xd->bd);
      } else {
        const PLANE_TYPE plane_type = pd->plane_type;
        switch (tx_size) {
          case TX_4X4:
            tx_type = get_tx_type_4x4(plane_type, xd, block);
#if CONFIG_TX_SKIP
            if (mbmi->tx_skip[plane != 0]) {
              vp9_tx_identity_add(dqcoeff, dst, stride, 4, shift);
            } else {
              vp9_highbd_iht4x4_add(tx_type, dqcoeff, dst, stride, eob, xd->bd);
            }
#else
            vp9_highbd_iht4x4_add(tx_type, dqcoeff, dst, stride, eob, xd->bd);
#endif
            break;
          case TX_8X8:
            tx_type = get_tx_type(plane_type, xd);
#if CONFIG_TX_SKIP
            if (mbmi->tx_skip[plane != 0]) {
              vp9_tx_identity_add(dqcoeff, dst, stride, 8, shift);
            } else {
              vp9_highbd_iht8x8_add(tx_type, dqcoeff, dst, stride, eob, xd->bd);
            }
#else
            vp9_highbd_iht8x8_add(tx_type, dqcoeff, dst, stride, eob, xd->bd);
#endif
            break;
          case TX_16X16:
            tx_type = get_tx_type(plane_type, xd);
#if CONFIG_TX_SKIP
            if (mbmi->tx_skip[plane != 0]) {
              vp9_tx_identity_add(dqcoeff, dst, stride, 16, shift);
            } else {
              vp9_highbd_iht16x16_add(tx_type, dqcoeff, dst, stride, eob,
                                      xd->bd);
            }
#else
            vp9_highbd_iht16x16_add(tx_type, dqcoeff, dst, stride, eob, xd->bd);
#endif
            break;
          case TX_32X32:
            tx_type = DCT_DCT;
#if CONFIG_TX_SKIP
            if (mbmi->tx_skip[plane != 0]) {
              vp9_tx_identity_add(dqcoeff, dst, stride, 32, shift);
            } else {
              vp9_highbd_idct32x32_add(dqcoeff, dst, stride, eob, xd->bd);
            }
#else
            vp9_highbd_idct32x32_add(dqcoeff, dst, stride, eob, xd->bd);
#endif
            break;
#if CONFIG_TX64X64
          case TX_64X64:
            tx_type = DCT_DCT;
#if CONFIG_TX_SKIP
            if (mbmi->tx_skip[plane != 0]) {
              vp9_tx_identity_add(dqcoeff, dst, stride, 64, shift);
            } else {
              vp9_highbd_idct64x64_add(dqcoeff, dst, stride, eob, xd->bd);
            }
#else
            vp9_highbd_idct64x64_add(dqcoeff, dst, stride, eob, xd->bd);
#endif  // CONFIG_TX_SKIP
            break;
#endif  // CONFIG_TX64X64
          default:
            assert(0 && "Invalid transform size");
        }
      }
    } else {
#if CONFIG_TX_SKIP
      if (xd->lossless && !mbmi->tx_skip[plane != 0]) {
#else
      if (xd->lossless) {
#endif
        tx_type = DCT_DCT;
        vp9_iwht4x4_add(dqcoeff, dst, stride, eob);
      } else {
        const PLANE_TYPE plane_type = pd->plane_type;
        switch (tx_size) {
          case TX_4X4:
            tx_type = get_tx_type_4x4(plane_type, xd, block);
#if CONFIG_TX_SKIP
            if (mbmi->tx_skip[plane != 0]) {
              vp9_tx_identity_add(dqcoeff, dst, stride, 4, shift);
            } else {
              vp9_iht4x4_add(tx_type, dqcoeff, dst, stride, eob);
            }
#else
            vp9_iht4x4_add(tx_type, dqcoeff, dst, stride, eob);
#endif
            break;
          case TX_8X8:
            tx_type = get_tx_type(plane_type, xd);
#if CONFIG_TX_SKIP
            if (mbmi->tx_skip[plane != 0]) {
              vp9_tx_identity_add(dqcoeff, dst, stride, 8, shift);
            } else {
              vp9_iht8x8_add(tx_type, dqcoeff, dst, stride, eob);
            }
#else
            vp9_iht8x8_add(tx_type, dqcoeff, dst, stride, eob);
#endif
            break;
          case TX_16X16:
            tx_type = get_tx_type(plane_type, xd);
#if CONFIG_TX_SKIP
            if (mbmi->tx_skip[plane != 0]) {
              vp9_tx_identity_add(dqcoeff, dst, stride, 16, shift);
            } else {
              vp9_iht16x16_add(tx_type, dqcoeff, dst, stride, eob);
            }
#else
            vp9_iht16x16_add(tx_type, dqcoeff, dst, stride, eob);
#endif
            break;
          case TX_32X32:
            tx_type = DCT_DCT;
#if CONFIG_TX_SKIP
            if (mbmi->tx_skip[plane != 0]) {
              vp9_tx_identity_add(dqcoeff, dst, stride, 32, shift);
            } else {
              vp9_idct32x32_add(dqcoeff, dst, stride, eob);
            }
#else
            vp9_idct32x32_add(dqcoeff, dst, stride, eob);
#endif
            break;
#if CONFIG_TX64X64
          case TX_64X64:
            tx_type = DCT_DCT;
#if CONFIG_TX_SKIP
            if (mbmi->tx_skip[plane != 0]) {
              vp9_tx_identity_add(dqcoeff, dst, stride, 64, shift);
            } else {
              vp9_idct64x64_add(dqcoeff, dst, stride, eob);
            }
#else
              vp9_idct64x64_add(dqcoeff, dst, stride, eob);
#endif  // CONFIG_TX_SKIP
            break;
#endif  // CONFIG_TX64X64
          default:
            assert(0 && "Invalid transform size");
            return;
        }
      }
    }

#else  // CONFIG_VP9_HIGHBITDEPTH

#if CONFIG_TX_SKIP
    if (xd->lossless && !mbmi->tx_skip[plane != 0]) {
#else
    if (xd->lossless) {
#endif
      tx_type = DCT_DCT;
      vp9_iwht4x4_add(dqcoeff, dst, stride, eob);
    } else {
      const PLANE_TYPE plane_type = pd->plane_type;
      switch (tx_size) {
        case TX_4X4:
          tx_type = get_tx_type_4x4(plane_type, xd, block);
#if CONFIG_TX_SKIP
          if (mbmi->tx_skip[plane != 0]) {
            if (mode == V_PRED || mode == H_PRED || mode == TM_PRED)
              vp9_intra_dpcm_add(dqcoeff, dst, stride, mode, 4, shift);
            else
              vp9_tx_identity_add(dqcoeff, dst, stride, 4, shift);
          } else {
            vp9_iht4x4_add(tx_type, dqcoeff, dst, stride, eob);
          }
#else
          vp9_iht4x4_add(tx_type, dqcoeff, dst, stride, eob);
#endif
          break;
        case TX_8X8:
          tx_type = get_tx_type(plane_type, xd);
#if CONFIG_TX_SKIP
          if (mbmi->tx_skip[plane != 0]) {
            if (mode == V_PRED || mode == H_PRED || mode == TM_PRED)
              vp9_intra_dpcm_add(dqcoeff, dst, stride, mode, 8, shift);
            else
              vp9_tx_identity_add(dqcoeff, dst, stride, 8, shift);
          } else {
            vp9_iht8x8_add(tx_type, dqcoeff, dst, stride, eob);
          }
#else
          vp9_iht8x8_add(tx_type, dqcoeff, dst, stride, eob);
#endif
          break;
        case TX_16X16:
          tx_type = get_tx_type(plane_type, xd);
#if CONFIG_TX_SKIP
          if (mbmi->tx_skip[plane != 0]) {
            if (mode == V_PRED || mode == H_PRED || mode == TM_PRED)
              vp9_intra_dpcm_add(dqcoeff, dst, stride, mode, 16, shift);
            else
              vp9_tx_identity_add(dqcoeff, dst, stride, 16, shift);
          } else {
            vp9_iht16x16_add(tx_type, dqcoeff, dst, stride, eob);
          }
#else
          vp9_iht16x16_add(tx_type, dqcoeff, dst, stride, eob);
#endif
          break;
        case TX_32X32:
          tx_type = DCT_DCT;
#if CONFIG_TX_SKIP
          if (mbmi->tx_skip[plane != 0]) {
            if (mode == V_PRED || mode == H_PRED || mode == TM_PRED)
              vp9_intra_dpcm_add(dqcoeff, dst, stride, mode, 32, shift);
            else
              vp9_tx_identity_add(dqcoeff, dst, stride, 32, shift);
          } else {
            vp9_idct32x32_add(dqcoeff, dst, stride, eob);;
          }
#else
          vp9_idct32x32_add(dqcoeff, dst, stride, eob);
#endif
          break;
#if CONFIG_TX64X64
        case TX_64X64:
          tx_type = DCT_DCT;
#if CONFIG_TX_SKIP
          if (mbmi->tx_skip[plane != 0]) {
            vp9_tx_identity_add(dqcoeff, dst, stride, 64, shift);
          } else {
            vp9_idct64x64_add(dqcoeff, dst, stride, eob);;
          }
#else
          vp9_idct64x64_add(dqcoeff, dst, stride, eob);
#endif  // CONFIG_TX_SKIP
          break;
#endif  // CONFIG_TX64X64
        default:
          assert(0 && "Invalid transform size");
          return;
      }
    }
#endif  // CONFIG_VP9_HIGHBITDEPTH

    if (eob == 1) {
      vpx_memset(dqcoeff, 0, 2 * sizeof(dqcoeff[0]));
    } else {
      if (tx_type == DCT_DCT && tx_size <= TX_16X16 && eob <= 10)
        vpx_memset(dqcoeff, 0, 4 * (4 << tx_size) * sizeof(dqcoeff[0]));
      else if (tx_size == TX_32X32 && eob <= 34)
        vpx_memset(dqcoeff, 0, 256 * sizeof(dqcoeff[0]));
      else
        vpx_memset(dqcoeff, 0, (16 << (tx_size << 1)) * sizeof(dqcoeff[0]));
    }
  }
}

struct intra_args {
  VP9_COMMON *cm;
  MACROBLOCKD *xd;
  vp9_reader *r;
};

static void predict_and_reconstruct_intra_block(int plane, int block,
                                                BLOCK_SIZE plane_bsize,
                                                TX_SIZE tx_size, void *arg) {
  struct intra_args *const args = (struct intra_args *)arg;
  VP9_COMMON *const cm = args->cm;
  MACROBLOCKD *const xd = args->xd;
  struct macroblockd_plane *const pd = &xd->plane[plane];
  MODE_INFO *const mi = xd->mi[0].src_mi;
  const PREDICTION_MODE mode = (plane == 0) ? get_y_mode(mi, block)
                                            : mi->mbmi.uv_mode;
  int x, y;
  uint8_t *dst;
#if CONFIG_TX_SKIP
  int no_coeff = 0;
#endif
#if CONFIG_FILTERINTRA
  int fbit;
  if (plane == 0)
    if (mi->mbmi.sb_type < BLOCK_8X8)
      fbit = mi->b_filter_info[block];
    else
      fbit = is_filter_enabled(tx_size) ? mi->mbmi.filterbit : 0;
  else
    fbit = is_filter_enabled(tx_size) ? mi->mbmi.uv_filterbit : 0;
#endif
  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &x, &y);
  dst = &pd->dst.buf[4 * y * pd->dst.stride + 4 * x];

  vp9_predict_intra_block(xd, block >> (tx_size << 1),
                          b_width_log2_lookup[plane_bsize], tx_size, mode,
#if CONFIG_FILTERINTRA
                          fbit,
#endif
                          dst, pd->dst.stride, dst, pd->dst.stride,
                          x, y, plane);
  if (!mi->mbmi.skip) {
    const int eob = vp9_decode_block_tokens(cm, xd, plane, block,
                                            plane_bsize, x, y, tx_size,
                                            args->r);
    inverse_transform_block(xd, plane, block, tx_size, dst, pd->dst.stride,
                            eob);
#if CONFIG_TX_SKIP
    no_coeff = !eob;
#endif
  }

#if CONFIG_TX_SKIP
  if ((mi->mbmi.skip || no_coeff) && mi->mbmi.tx_skip[plane != 0] &&
      mode == TM_PRED && tx_size <= TX_32X32) {
    int bs = 4 * (1 << tx_size);
    vp9_intra_dpcm_add_nocoeff(dst, pd->dst.stride, mode, bs);
  }
#endif

#if CONFIG_TX_SKIP && CONFIG_FILTERINTRA
  if ((mi->mbmi.skip || no_coeff) && mi->mbmi.tx_skip[plane != 0] &&
      (mode == H_PRED || mode == V_PRED) && fbit) {
    int bs = 4 * (1 << tx_size);
    vp9_intra_dpcm_add_nocoeff(dst, pd->dst.stride, mode, bs);
  }
#endif
}

struct inter_args {
  VP9_COMMON *cm;
  MACROBLOCKD *xd;
  vp9_reader *r;
  int *eobtotal;
};

static void reconstruct_inter_block(int plane, int block,
                                    BLOCK_SIZE plane_bsize,
                                    TX_SIZE tx_size, void *arg) {
  struct inter_args *args = (struct inter_args *)arg;
  VP9_COMMON *const cm = args->cm;
  MACROBLOCKD *const xd = args->xd;
  struct macroblockd_plane *const pd = &xd->plane[plane];
  int x, y, eob;
  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &x, &y);
  eob = vp9_decode_block_tokens(cm, xd, plane, block, plane_bsize, x, y,
                                tx_size, args->r);
  inverse_transform_block(xd, plane, block, tx_size,
                          &pd->dst.buf[4 * y * pd->dst.stride + 4 * x],
                          pd->dst.stride, eob);
  *args->eobtotal += eob;
}

static MB_MODE_INFO *set_offsets(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                                 const TileInfo *const tile,
                                 BLOCK_SIZE bsize, int mi_row, int mi_col) {
  const int bw = num_8x8_blocks_wide_lookup[bsize];
  const int bh = num_8x8_blocks_high_lookup[bsize];
  const int x_mis = MIN(bw, cm->mi_cols - mi_col);
  const int y_mis = MIN(bh, cm->mi_rows - mi_row);
  const int offset = mi_row * cm->mi_stride + mi_col;
  int x, y;

  xd->mi = cm->mi + offset;
  xd->mi[0].src_mi = &xd->mi[0];  // Point to self.
  xd->mi[0].mbmi.sb_type = bsize;

  for (y = 0; y < y_mis; ++y)
    for (x = !y; x < x_mis; ++x) {
      xd->mi[y * cm->mi_stride + x].src_mi = &xd->mi[0];
    }

  set_skip_context(xd, mi_row, mi_col);

  // Distance of Mb to the various image edges. These are specified to 8th pel
  // as they are always compared to values that are in 1/8th pel units
  set_mi_row_col(xd, tile, mi_row, bh, mi_col, bw, cm->mi_rows, cm->mi_cols);

  vp9_setup_dst_planes(xd->plane, get_frame_new_buffer(cm), mi_row, mi_col);
  return &xd->mi[0].mbmi;
}

#if CONFIG_SUPERTX
static MB_MODE_INFO *set_offsets_extend(VP9_COMMON *const cm,
                                        MACROBLOCKD *const xd,
                                        const TileInfo *const tile,
                                        BLOCK_SIZE top_bsize,
                                        int mi_row, int mi_col,
                                        int mi_row_ori, int mi_col_ori) {
  const int bw = num_8x8_blocks_wide_lookup[top_bsize];
  const int bh = num_8x8_blocks_high_lookup[top_bsize];
  const int offset = mi_row * cm->mi_stride + mi_col;

  xd->mi = cm->mi + offset;
  xd->mi[0].src_mi = &xd->mi[0];
  set_mi_row_col(xd, tile, mi_row_ori, bh, mi_col_ori, bw,
                 cm->mi_rows, cm->mi_cols);
  return &xd->mi[0].mbmi;
}

static MB_MODE_INFO *set_mb_offsets(VP9_COMMON *const cm,
                                    MACROBLOCKD *const xd,
                                    const TileInfo *const tile,
                                    BLOCK_SIZE bsize,
                                    int mi_row, int mi_col) {
  const int bw = num_8x8_blocks_wide_lookup[bsize];
  const int bh = num_8x8_blocks_high_lookup[bsize];
  const int x_mis = MIN(bw, cm->mi_cols - mi_col);
  const int y_mis = MIN(bh, cm->mi_rows - mi_row);
  const int offset = mi_row * cm->mi_stride + mi_col;
  int x, y;

  xd->mi = cm->mi + offset;
  xd->mi[0].src_mi = &xd->mi[0];
  xd->mi[0].mbmi.sb_type = bsize;
  for (y = 0; y < y_mis; ++y)
    for (x = !y; x < x_mis; ++x)
      xd->mi[y * cm->mi_stride + x] = xd->mi[0];

  set_mi_row_col(xd, tile, mi_row, bh, mi_col, bw, cm->mi_rows, cm->mi_cols);
  return &xd->mi[0].mbmi;
}

static void set_offsets_topblock(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                                 const TileInfo *const tile,
                                 BLOCK_SIZE bsize, int mi_row, int mi_col) {
  const int bw = num_8x8_blocks_wide_lookup[bsize];
  const int bh = num_8x8_blocks_high_lookup[bsize];
  const int offset = mi_row * cm->mi_stride + mi_col;

  xd->mi = cm->mi + offset;
  xd->mi[0].src_mi = &xd->mi[0];

  set_mi_row_col(xd, tile, mi_row, bh, mi_col, bw, cm->mi_rows, cm->mi_cols);

  vp9_setup_dst_planes(xd->plane, get_frame_new_buffer(cm), mi_row, mi_col);
}

static void set_param_topblock(VP9_COMMON *const cm,  MACROBLOCKD *const xd,
                               BLOCK_SIZE bsize, int mi_row, int mi_col,
#if CONFIG_EXT_TX
                               int txfm,
#endif
                               int skip) {
  const int bw = num_8x8_blocks_wide_lookup[bsize];
  const int bh = num_8x8_blocks_high_lookup[bsize];
  const int x_mis = MIN(bw, cm->mi_cols - mi_col);
  const int y_mis = MIN(bh, cm->mi_rows - mi_row);
  const int offset = mi_row * cm->mi_stride + mi_col;
  int x, y;

  xd->mi = cm->mi + offset;
  xd->mi[0].src_mi = &xd->mi[0];

  for (y = 0; y < y_mis; ++y)
    for (x = 0; x < x_mis; ++x) {
      xd->mi[y * cm->mi_stride + x].mbmi.skip = skip;
#if CONFIG_EXT_TX
      xd->mi[y * cm->mi_stride + x].mbmi.ext_txfrm = txfm;
#endif
    }
}

static void set_ref(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                    int idx, int mi_row, int mi_col) {
  MB_MODE_INFO *const mbmi = &xd->mi[0].mbmi;
  RefBuffer *ref_buffer = &cm->frame_refs[mbmi->ref_frame[idx] - LAST_FRAME];
  xd->block_refs[idx] = ref_buffer;
  if (!vp9_is_valid_scale(&ref_buffer->sf))
    vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                       "Invalid scale factors");
  vp9_setup_pre_planes(xd, idx, ref_buffer->buf, mi_row, mi_col,
                       &ref_buffer->sf);
  xd->corrupted |= ref_buffer->buf->corrupted;
}

static void dec_predict_b_extend(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                                 const TileInfo *const tile,
                                 int mi_row, int mi_col,
                                 int mi_row_ori, int mi_col_ori,
                                 BLOCK_SIZE top_bsize) {
  MB_MODE_INFO *mbmi = set_offsets_extend(cm, xd, tile, top_bsize,
                                          mi_row, mi_col,
                                          mi_row_ori, mi_col_ori);
  set_ref(cm, xd, 0, mi_row_ori, mi_col_ori);
  if (has_second_ref(&xd->mi[0].mbmi))
    set_ref(cm, xd, 1, mi_row_ori, mi_col_ori);
  mbmi->tx_size = b_width_log2_lookup[top_bsize];
#if CONFIG_WEDGE_PARTITION
  vp9_dec_build_inter_predictors_sb_extend(xd, mi_row, mi_col,
                                           mi_row_ori, mi_col_ori, top_bsize);
#else
  vp9_dec_build_inter_predictors_sb(xd, mi_row_ori, mi_col_ori, top_bsize);
#endif  // CONFIG_WEDGE_PARTITION
}

static void dec_predict_b_sub8x8_extend(VP9_COMMON *const cm,
                                        MACROBLOCKD *const xd,
                                        const TileInfo *const tile,
                                        int mi_row, int mi_col,
                                        int mi_row_ori, int mi_col_ori,
                                        BLOCK_SIZE top_bsize,
                                        PARTITION_TYPE partition) {
  MB_MODE_INFO *mbmi = set_offsets_extend(cm, xd, tile, top_bsize,
                                          mi_row, mi_col,
                                          mi_row_ori, mi_col_ori);
  set_ref(cm, xd, 0, mi_row_ori, mi_col_ori);
  if (has_second_ref(&xd->mi[0].mbmi))
    set_ref(cm, xd, 1, mi_row_ori, mi_col_ori);
  mbmi->tx_size = b_width_log2_lookup[top_bsize];
  vp9_dec_build_inter_predictors_sby_sub8x8_extend(xd, mi_row, mi_col,
                                                   mi_row_ori, mi_col_ori,
                                                   top_bsize, partition);
  vp9_dec_build_inter_predictors_sbuv_sub8x8_extend(xd,
#if CONFIG_WEDGE_PARTITION
                                                    mi_row, mi_col,
#endif
                                                    mi_row_ori, mi_col_ori,
                                                    top_bsize);
}

static void dec_predict_sb_complex(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                                   const TileInfo *const tile,
                                   int mi_row, int mi_col,
                                   int mi_row_ori, int mi_col_ori,
                                   BLOCK_SIZE bsize, BLOCK_SIZE top_bsize,
                                   uint8_t *dst_buf[3], int dst_stride[3]) {
  const int bsl = b_width_log2_lookup[bsize], hbs = (1 << bsl) / 4;
  PARTITION_TYPE partition;
  BLOCK_SIZE subsize;
  MB_MODE_INFO *mbmi;
  int i, offset = mi_row * cm->mi_stride + mi_col;

  DECLARE_ALIGNED_ARRAY(16, uint8_t, tmp_buf1,
                        MAX_MB_PLANE * MAXTXLEN * MAXTXLEN);
  DECLARE_ALIGNED_ARRAY(16, uint8_t, tmp_buf2,
                        MAX_MB_PLANE * MAXTXLEN * MAXTXLEN);
  DECLARE_ALIGNED_ARRAY(16, uint8_t, tmp_buf3,
                        MAX_MB_PLANE * MAXTXLEN * MAXTXLEN);
  uint8_t *dst_buf1[3] = {
    tmp_buf1,
    tmp_buf1 + MAXTXLEN * MAXTXLEN,
    tmp_buf1 + 2 * MAXTXLEN * MAXTXLEN};
  uint8_t *dst_buf2[3] = {
    tmp_buf2,
    tmp_buf2 + MAXTXLEN * MAXTXLEN,
    tmp_buf2 + 2 * MAXTXLEN * MAXTXLEN};
  uint8_t *dst_buf3[3] = {
    tmp_buf3,
    tmp_buf3 + MAXTXLEN * MAXTXLEN,
    tmp_buf3 + 2 * MAXTXLEN * MAXTXLEN};
  int dst_stride1[3] = {MAXTXLEN, MAXTXLEN, MAXTXLEN};
  int dst_stride2[3] = {MAXTXLEN, MAXTXLEN, MAXTXLEN};
  int dst_stride3[3] = {MAXTXLEN, MAXTXLEN, MAXTXLEN};

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  xd->mi = cm->mi + offset;
  xd->mi[0].src_mi = &xd->mi[0];
  mbmi = &xd->mi[0].mbmi;
  partition = partition_lookup[bsl][mbmi->sb_type];
  subsize = get_subsize(bsize, partition);

  for (i = 0; i < MAX_MB_PLANE; i++) {
    xd->plane[i].dst.buf = dst_buf[i];
    xd->plane[i].dst.stride = dst_stride[i];
  }

  switch (partition) {
    case PARTITION_NONE:
      assert(bsize < top_bsize);
      dec_predict_b_extend(cm, xd, tile, mi_row, mi_col, mi_row_ori, mi_col_ori,
                           top_bsize);
      break;
    case PARTITION_HORZ:
      if (bsize > BLOCK_8X8) {
        dec_predict_b_extend(cm, xd, tile, mi_row, mi_col, mi_row_ori,
                             mi_col_ori, top_bsize);
      } else {
        dec_predict_b_sub8x8_extend(cm, xd, tile, mi_row, mi_col,
                                    mi_row_ori, mi_col_ori,
                                    top_bsize, partition);
      }
      if (mi_row + hbs < cm->mi_rows && bsize > BLOCK_8X8) {
        for (i = 0; i < MAX_MB_PLANE; i++) {
          xd->plane[i].dst.buf = dst_buf1[i];
          xd->plane[i].dst.stride = dst_stride1[i];
        }
        dec_predict_b_extend(cm, xd, tile, mi_row + hbs, mi_col,
                             mi_row_ori, mi_col_ori, top_bsize);
        for (i = 0; i < MAX_MB_PLANE; i++) {
          xd->plane[i].dst.buf = dst_buf[i];
          xd->plane[i].dst.stride = dst_stride[i];
          vp9_build_masked_inter_predictor_complex(xd,
                                                   dst_buf[i], dst_stride[i],
                                                   dst_buf1[i], dst_stride1[i],
                                                   &xd->plane[i],
                                                   mi_row, mi_col,
                                                   mi_row_ori, mi_col_ori,
                                                   bsize, top_bsize,
                                                   PARTITION_HORZ);
        }
      }
      break;
    case PARTITION_VERT:
      if (bsize > BLOCK_8X8) {
        dec_predict_b_extend(cm, xd, tile, mi_row, mi_col, mi_row_ori,
                             mi_col_ori, top_bsize);
      } else {
        dec_predict_b_sub8x8_extend(cm, xd, tile, mi_row, mi_col,
                                    mi_row_ori, mi_col_ori,
                                    top_bsize, partition);
      }
      if (mi_col + hbs < cm->mi_cols && bsize > BLOCK_8X8) {
        for (i = 0; i < MAX_MB_PLANE; i++) {
          xd->plane[i].dst.buf = dst_buf1[i];
          xd->plane[i].dst.stride = dst_stride1[i];
        }
        dec_predict_b_extend(cm, xd, tile, mi_row, mi_col + hbs, mi_row_ori,
                             mi_col_ori, top_bsize);
        for (i = 0; i < MAX_MB_PLANE; i++) {
          xd->plane[i].dst.buf = dst_buf[i];
          xd->plane[i].dst.stride = dst_stride[i];
          vp9_build_masked_inter_predictor_complex(xd,
                                                   dst_buf[i], dst_stride[i],
                                                   dst_buf1[i], dst_stride1[i],
                                                   &xd->plane[i],
                                                   mi_row, mi_col,
                                                   mi_row_ori, mi_col_ori,
                                                   bsize, top_bsize,
                                                   PARTITION_VERT);
        }
      }
      break;
    case PARTITION_SPLIT:
      if (bsize == BLOCK_8X8) {
        dec_predict_b_sub8x8_extend(cm, xd, tile, mi_row, mi_col,
                                    mi_row_ori, mi_col_ori,
                                    top_bsize, partition);
      } else {
        dec_predict_sb_complex(cm, xd, tile, mi_row, mi_col,
                               mi_row_ori, mi_col_ori, subsize, top_bsize,
                               dst_buf, dst_stride);
        if (mi_row < cm->mi_rows && mi_col + hbs < cm->mi_cols)
          dec_predict_sb_complex(cm, xd, tile, mi_row, mi_col + hbs,
                                 mi_row_ori, mi_col_ori, subsize, top_bsize,
                                 dst_buf1, dst_stride1);
        if (mi_row + hbs < cm->mi_rows && mi_col < cm->mi_cols)
          dec_predict_sb_complex(cm, xd, tile, mi_row + hbs, mi_col,
                                 mi_row_ori, mi_col_ori, subsize, top_bsize,
                                 dst_buf2, dst_stride2);
        if (mi_row + hbs < cm->mi_rows && mi_col + hbs < cm->mi_cols)
          dec_predict_sb_complex(cm, xd, tile, mi_row + hbs, mi_col + hbs,
                                 mi_row_ori, mi_col_ori, subsize, top_bsize,
                                 dst_buf3, dst_stride3);
        for (i = 0; i < MAX_MB_PLANE; i++) {
          if (mi_row < cm->mi_rows && mi_col + hbs < cm->mi_cols) {
            vp9_build_masked_inter_predictor_complex(xd,
                                                     dst_buf[i], dst_stride[i],
                                                     dst_buf1[i],
                                                     dst_stride1[i],
                                                     &xd->plane[i],
                                                     mi_row, mi_col,
                                                     mi_row_ori, mi_col_ori,
                                                     bsize, top_bsize,
                                                     PARTITION_VERT);
            if (mi_row + hbs < cm->mi_rows) {
              vp9_build_masked_inter_predictor_complex(xd,
                                                       dst_buf2[i],
                                                       dst_stride2[i],
                                                       dst_buf3[i],
                                                       dst_stride3[i],
                                                       &xd->plane[i],
                                                       mi_row, mi_col,
                                                       mi_row_ori, mi_col_ori,
                                                       bsize, top_bsize,
                                                       PARTITION_VERT);
              vp9_build_masked_inter_predictor_complex(xd,
                                                       dst_buf[i],
                                                       dst_stride[i],
                                                       dst_buf2[i],
                                                       dst_stride2[i],
                                                       &xd->plane[i],
                                                       mi_row, mi_col,
                                                       mi_row_ori, mi_col_ori,
                                                       bsize, top_bsize,
                                                       PARTITION_HORZ);
            }
          } else if (mi_row + hbs < cm->mi_rows && mi_col < cm->mi_cols) {
            vp9_build_masked_inter_predictor_complex(xd,
                                                     dst_buf[i],
                                                     dst_stride[i],
                                                     dst_buf2[i],
                                                     dst_stride2[i],
                                                     &xd->plane[i],
                                                     mi_row, mi_col,
                                                     mi_row_ori, mi_col_ori,
                                                     bsize, top_bsize,
                                                     PARTITION_HORZ);
          }
        }
      }
      break;
    default:
      assert(0);
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
static void dec_predict_sb_complex_highbd(
    VP9_COMMON *const cm, MACROBLOCKD *const xd,
    const TileInfo *const tile,
    int mi_row, int mi_col,
    int mi_row_ori, int mi_col_ori,
    BLOCK_SIZE bsize, BLOCK_SIZE top_bsize,
    uint8_t *dst_buf[3], int dst_stride[3]) {
  const int bsl = b_width_log2_lookup[bsize], hbs = (1 << bsl) / 4;
  PARTITION_TYPE partition;
  BLOCK_SIZE subsize;
  MB_MODE_INFO *mbmi;
  int i, offset = mi_row * cm->mi_stride + mi_col;

  DECLARE_ALIGNED_ARRAY(16, uint8_t, tmp_buf1,
                        MAX_MB_PLANE * MAXTXLEN * MAXTXLEN * sizeof(uint16_t));
  DECLARE_ALIGNED_ARRAY(16, uint8_t, tmp_buf2,
                        MAX_MB_PLANE * MAXTXLEN * MAXTXLEN * sizeof(uint16_t));
  DECLARE_ALIGNED_ARRAY(16, uint8_t, tmp_buf3,
                        MAX_MB_PLANE * MAXTXLEN * MAXTXLEN * sizeof(uint16_t));
  uint8_t *dst_buf1[3] = {
    CONVERT_TO_BYTEPTR(tmp_buf1),
    CONVERT_TO_BYTEPTR(tmp_buf1 + MAXTXLEN * MAXTXLEN * sizeof(uint16_t)),
    CONVERT_TO_BYTEPTR(tmp_buf1 + 2 * MAXTXLEN * MAXTXLEN * sizeof(uint16_t))};
  uint8_t *dst_buf2[3] = {
    CONVERT_TO_BYTEPTR(tmp_buf2),
    CONVERT_TO_BYTEPTR(tmp_buf2 + MAXTXLEN * MAXTXLEN * sizeof(uint16_t)),
    CONVERT_TO_BYTEPTR(tmp_buf2 + 2 * MAXTXLEN * MAXTXLEN * sizeof(uint16_t))};
  uint8_t *dst_buf3[3] = {
    CONVERT_TO_BYTEPTR(tmp_buf3),
    CONVERT_TO_BYTEPTR(tmp_buf3 + MAXTXLEN * MAXTXLEN * sizeof(uint16_t)),
    CONVERT_TO_BYTEPTR(tmp_buf3 + 2 * MAXTXLEN * MAXTXLEN * sizeof(uint16_t))};
  int dst_stride1[3] = {MAXTXLEN, MAXTXLEN, MAXTXLEN};
  int dst_stride2[3] = {MAXTXLEN, MAXTXLEN, MAXTXLEN};
  int dst_stride3[3] = {MAXTXLEN, MAXTXLEN, MAXTXLEN};

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  xd->mi = cm->mi + offset;
  xd->mi[0].src_mi = &xd->mi[0];
  mbmi = &xd->mi[0].mbmi;
  partition = partition_lookup[bsl][mbmi->sb_type];
  subsize = get_subsize(bsize, partition);

  for (i = 0; i < MAX_MB_PLANE; i++) {
    xd->plane[i].dst.buf = dst_buf[i];
    xd->plane[i].dst.stride = dst_stride[i];
  }

  switch (partition) {
    case PARTITION_NONE:
      assert(bsize < top_bsize);
      dec_predict_b_extend(cm, xd, tile, mi_row, mi_col, mi_row_ori, mi_col_ori,
                           top_bsize);
      break;
    case PARTITION_HORZ:
      if (bsize > BLOCK_8X8) {
        dec_predict_b_extend(cm, xd, tile, mi_row, mi_col, mi_row_ori,
                             mi_col_ori, top_bsize);
      } else {
        dec_predict_b_sub8x8_extend(cm, xd, tile, mi_row, mi_col,
                                    mi_row_ori, mi_col_ori,
                                    top_bsize, partition);
      }
      if (mi_row + hbs < cm->mi_rows && bsize > BLOCK_8X8) {
        for (i = 0; i < MAX_MB_PLANE; i++) {
          xd->plane[i].dst.buf = dst_buf1[i];
          xd->plane[i].dst.stride = dst_stride1[i];
        }
        dec_predict_b_extend(cm, xd, tile, mi_row + hbs, mi_col,
                             mi_row_ori, mi_col_ori, top_bsize);
        for (i = 0; i < MAX_MB_PLANE; i++) {
          xd->plane[i].dst.buf = dst_buf[i];
          xd->plane[i].dst.stride = dst_stride[i];
          vp9_build_masked_inter_predictor_complex(
              xd,
              dst_buf[i], dst_stride[i],
              dst_buf1[i], dst_stride1[i],
              &xd->plane[i],
              mi_row, mi_col,
              mi_row_ori, mi_col_ori,
              bsize, top_bsize,
              PARTITION_HORZ);
        }
      }
      break;
    case PARTITION_VERT:
      if (bsize > BLOCK_8X8) {
        dec_predict_b_extend(cm, xd, tile, mi_row, mi_col, mi_row_ori,
                             mi_col_ori, top_bsize);
      } else {
        dec_predict_b_sub8x8_extend(cm, xd, tile, mi_row, mi_col,
                                    mi_row_ori, mi_col_ori,
                                    top_bsize, partition);
      }
      if (mi_col + hbs < cm->mi_cols && bsize > BLOCK_8X8) {
        for (i = 0; i < MAX_MB_PLANE; i++) {
          xd->plane[i].dst.buf = dst_buf1[i];
          xd->plane[i].dst.stride = dst_stride1[i];
        }
        dec_predict_b_extend(cm, xd, tile, mi_row, mi_col + hbs, mi_row_ori,
                             mi_col_ori, top_bsize);
        for (i = 0; i < MAX_MB_PLANE; i++) {
          xd->plane[i].dst.buf = dst_buf[i];
          xd->plane[i].dst.stride = dst_stride[i];
          vp9_build_masked_inter_predictor_complex(
              xd,
              dst_buf[i], dst_stride[i],
              dst_buf1[i], dst_stride1[i],
              &xd->plane[i],
              mi_row, mi_col,
              mi_row_ori, mi_col_ori,
              bsize, top_bsize,
              PARTITION_VERT);
        }
      }
      break;
    case PARTITION_SPLIT:
      if (bsize == BLOCK_8X8) {
        dec_predict_b_sub8x8_extend(cm, xd, tile, mi_row, mi_col,
                                    mi_row_ori, mi_col_ori,
                                    top_bsize, partition);
      } else {
        dec_predict_sb_complex_highbd(cm, xd, tile, mi_row, mi_col,
                                      mi_row_ori, mi_col_ori, subsize,
                                      top_bsize, dst_buf, dst_stride);
        if (mi_row < cm->mi_rows && mi_col + hbs < cm->mi_cols)
          dec_predict_sb_complex_highbd(cm, xd, tile, mi_row, mi_col + hbs,
                                        mi_row_ori, mi_col_ori, subsize,
                                        top_bsize, dst_buf1, dst_stride1);
        if (mi_row + hbs < cm->mi_rows && mi_col < cm->mi_cols)
          dec_predict_sb_complex_highbd(cm, xd, tile, mi_row + hbs, mi_col,
                                        mi_row_ori, mi_col_ori, subsize,
                                        top_bsize, dst_buf2, dst_stride2);
        if (mi_row + hbs < cm->mi_rows && mi_col + hbs < cm->mi_cols)
          dec_predict_sb_complex_highbd(cm, xd, tile,
                                        mi_row + hbs, mi_col + hbs,
                                        mi_row_ori, mi_col_ori, subsize,
                                        top_bsize, dst_buf3, dst_stride3);
        for (i = 0; i < MAX_MB_PLANE; i++) {
          if (mi_row < cm->mi_rows && mi_col + hbs < cm->mi_cols) {
            vp9_build_masked_inter_predictor_complex(
                xd,
                dst_buf[i], dst_stride[i],
                dst_buf1[i],
                dst_stride1[i],
                &xd->plane[i],
                mi_row, mi_col,
                mi_row_ori, mi_col_ori,
                bsize, top_bsize,
                PARTITION_VERT);
            if (mi_row + hbs < cm->mi_rows) {
              vp9_build_masked_inter_predictor_complex(
                  xd,
                  dst_buf2[i],
                  dst_stride2[i],
                  dst_buf3[i],
                  dst_stride3[i],
                  &xd->plane[i],
                  mi_row, mi_col,
                  mi_row_ori, mi_col_ori,
                  bsize, top_bsize,
                  PARTITION_VERT);
              vp9_build_masked_inter_predictor_complex(
                  xd,
                  dst_buf[i],
                  dst_stride[i],
                  dst_buf2[i],
                  dst_stride2[i],
                  &xd->plane[i],
                  mi_row, mi_col,
                  mi_row_ori, mi_col_ori,
                  bsize, top_bsize,
                  PARTITION_HORZ);
            }
          } else if (mi_row + hbs < cm->mi_rows && mi_col < cm->mi_cols) {
            vp9_build_masked_inter_predictor_complex(
                xd,
                dst_buf[i],
                dst_stride[i],
                dst_buf2[i],
                dst_stride2[i],
                &xd->plane[i],
                mi_row, mi_col,
                mi_row_ori, mi_col_ori,
                bsize, top_bsize,
                PARTITION_HORZ);
          }
        }
      }
      break;
    default:
      assert(0);
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_SUPERTX

static void decode_block(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                         const TileInfo *const tile,
#if CONFIG_SUPERTX
                         int supertx_enabled,
#endif
                         int mi_row, int mi_col,
                         vp9_reader *r, BLOCK_SIZE bsize) {
  const int less8x8 = bsize < BLOCK_8X8;
#if CONFIG_TX_SKIP
  int q_idx;
#endif
#if CONFIG_SUPERTX
  MB_MODE_INFO *mbmi;
  if (supertx_enabled) {
    mbmi = set_mb_offsets(cm, xd, tile, bsize, mi_row, mi_col);
  } else {
    mbmi = set_offsets(cm, xd, tile, bsize, mi_row, mi_col);
  }
  vp9_read_mode_info(cm, xd, tile, supertx_enabled, mi_row, mi_col, r);
#else
  MB_MODE_INFO *mbmi = set_offsets(cm, xd, tile, bsize, mi_row, mi_col);
  vp9_read_mode_info(cm, xd, tile, mi_row, mi_col, r);
#endif  // CONFIG_SUPERTX
#if CONFIG_TX_SKIP
  q_idx = vp9_get_qindex(&cm->seg, mbmi->segment_id, cm->base_qindex);
  mbmi->tx_skip_shift = q_idx > TX_SKIP_SHIFT_THRESH ?
                        TX_SKIP_SHIFT_HQ : TX_SKIP_SHIFT_LQ;
#endif

#if CONFIG_SUPERTX
  if (!supertx_enabled) {
#endif
  if (less8x8)
    bsize = BLOCK_8X8;

  if (mbmi->skip) {
    reset_skip_context(xd, bsize);
  } else {
    if (cm->seg.enabled) {
      setup_plane_dequants(cm, xd, vp9_get_qindex(&cm->seg, mbmi->segment_id,
                                                  cm->base_qindex));
    }
  }
  if (!is_inter_block(mbmi)) {
    struct intra_args arg = { cm, xd, r };
    vp9_foreach_transformed_block(xd, bsize,
                                  predict_and_reconstruct_intra_block, &arg);
  } else {
    // Prediction
    vp9_dec_build_inter_predictors_sb(xd, mi_row, mi_col, bsize);

    // Reconstruction
    if (!mbmi->skip) {
      int eobtotal = 0;
      struct inter_args arg = { cm, xd, r, &eobtotal };

      vp9_foreach_transformed_block(xd, bsize, reconstruct_inter_block, &arg);
      if (!less8x8 && eobtotal == 0)
        mbmi->skip = 1;  // skip loopfilter
    }
  }
#if CONFIG_SUPERTX
  }
#endif

  xd->corrupted |= vp9_reader_has_error(r);
}

static PARTITION_TYPE read_partition(VP9_COMMON *cm, MACROBLOCKD *xd, int hbs,
                                     int mi_row, int mi_col, BLOCK_SIZE bsize,
                                     vp9_reader *r) {
  const int ctx = partition_plane_context(xd, mi_row, mi_col, bsize);
  const vp9_prob *const probs = get_partition_probs(cm, ctx);
  const int has_rows = (mi_row + hbs) < cm->mi_rows;
  const int has_cols = (mi_col + hbs) < cm->mi_cols;
  PARTITION_TYPE p;

  if (has_rows && has_cols)
    p = (PARTITION_TYPE)vp9_read_tree(r, vp9_partition_tree, probs);
  else if (!has_rows && has_cols)
    p = vp9_read(r, probs[1]) ? PARTITION_SPLIT : PARTITION_HORZ;
  else if (has_rows && !has_cols)
    p = vp9_read(r, probs[2]) ? PARTITION_SPLIT : PARTITION_VERT;
  else
    p = PARTITION_SPLIT;

  if (!cm->frame_parallel_decoding_mode)
    ++cm->counts.partition[ctx][p];

  return p;
}

#if CONFIG_SUPERTX
static int read_skip_without_seg(VP9_COMMON *cm, const MACROBLOCKD *xd,
                                 vp9_reader *r) {
  const int ctx = vp9_get_skip_context(xd);
  const int skip = vp9_read(r, cm->fc.skip_probs[ctx]);
  if (!cm->frame_parallel_decoding_mode)
    ++cm->counts.skip[ctx][skip];
  return skip;
}
#endif  // CONFIG_SUPERTX

static void decode_partition(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                             const TileInfo *const tile,
#if CONFIG_SUPERTX
                             int supertx_enabled,
#endif
                             int mi_row, int mi_col,
                             vp9_reader* r, BLOCK_SIZE bsize) {
  const int hbs = num_8x8_blocks_wide_lookup[bsize] / 2;
  PARTITION_TYPE partition;
  BLOCK_SIZE subsize, uv_subsize;
#if CONFIG_SUPERTX
  const int read_token = !supertx_enabled;
  int skip = 0;
  TX_SIZE supertx_size = b_width_log2_lookup[bsize];
#if CONFIG_EXT_TX
  int txfm = NORM;
#endif
#endif  // CONFIG_SUPERTX

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  partition = read_partition(cm, xd, hbs, mi_row, mi_col, bsize, r);
  subsize = get_subsize(bsize, partition);
  uv_subsize = ss_size_lookup[subsize][cm->subsampling_x][cm->subsampling_y];
  if (subsize >= BLOCK_8X8 && uv_subsize == BLOCK_INVALID)
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Invalid block size.");
#if CONFIG_SUPERTX
  if (cm->frame_type != KEY_FRAME &&
      partition != PARTITION_NONE &&
      bsize <= MAX_SUPERTX_BLOCK_SIZE &&
      !supertx_enabled && !xd->lossless) {
    const int supertx_context =
        partition_supertx_context_lookup[partition];
    supertx_enabled = vp9_read(
        r, cm->fc.supertx_prob[supertx_context][supertx_size]);
    cm->counts.supertx[supertx_context][supertx_size][supertx_enabled]++;
  }
  if (supertx_enabled && read_token) {
    int offset = mi_row * cm->mi_stride + mi_col;
    xd->mi = cm->mi + offset;
    xd->mi[0].src_mi = &xd->mi[0];
    set_mi_row_col(xd, tile, mi_row, num_8x8_blocks_high_lookup[bsize],
                   mi_col, num_8x8_blocks_wide_lookup[bsize],
                   cm->mi_rows, cm->mi_cols);
    set_skip_context(xd, mi_row, mi_col);
    // Here skip is read without using any segment level feature
    skip = read_skip_without_seg(cm, xd, r);
    if (skip)
      reset_skip_context(xd, bsize);
#if CONFIG_EXT_TX
    if (bsize <= BLOCK_16X16 && !skip) {
      txfm = vp9_read_tree(r, vp9_ext_tx_tree,
                           cm->fc.ext_tx_prob[supertx_size]);
      if (!cm->frame_parallel_decoding_mode)
        ++cm->counts.ext_tx[supertx_size][txfm];
    }
#endif
  }
#endif  // CONFIG_SUPERTX
  if (subsize < BLOCK_8X8) {
    decode_block(cm, xd, tile,
#if CONFIG_SUPERTX
                 supertx_enabled,
#endif
                 mi_row, mi_col, r, subsize);
  } else {
    switch (partition) {
      case PARTITION_NONE:
        decode_block(cm, xd, tile,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row, mi_col, r, subsize);
        break;
      case PARTITION_HORZ:
        decode_block(cm, xd, tile,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row, mi_col, r, subsize);
        if (mi_row + hbs < cm->mi_rows)
          decode_block(cm, xd, tile,
#if CONFIG_SUPERTX
                       supertx_enabled,
#endif
                       mi_row + hbs, mi_col, r, subsize);
        break;
      case PARTITION_VERT:
        decode_block(cm, xd, tile,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row, mi_col, r, subsize);
        if (mi_col + hbs < cm->mi_cols)
          decode_block(cm, xd, tile,
#if CONFIG_SUPERTX
                       supertx_enabled,
#endif
                       mi_row, mi_col + hbs, r, subsize);
        break;
      case PARTITION_SPLIT:
#if CONFIG_SUPERTX
        decode_partition(cm, xd, tile, supertx_enabled,
                         mi_row, mi_col, r, subsize);
        decode_partition(cm, xd, tile, supertx_enabled,
                         mi_row, mi_col + hbs, r, subsize);
        decode_partition(cm, xd, tile, supertx_enabled,
                         mi_row + hbs, mi_col, r, subsize);
        decode_partition(cm, xd, tile, supertx_enabled,
                         mi_row + hbs, mi_col + hbs, r, subsize);
#else
        decode_partition(cm, xd, tile, mi_row,       mi_col,       r, subsize);
        decode_partition(cm, xd, tile, mi_row,       mi_col + hbs, r, subsize);
        decode_partition(cm, xd, tile, mi_row + hbs, mi_col,       r, subsize);
        decode_partition(cm, xd, tile, mi_row + hbs, mi_col + hbs, r, subsize);
#endif
        break;
      default:
        assert(0 && "Invalid partition type");
    }
  }

#if CONFIG_SUPERTX
  if (supertx_enabled && read_token) {
    uint8_t *dst_buf[3];
    int dst_stride[3], i;

    vp9_setup_dst_planes(xd->plane, get_frame_new_buffer(cm), mi_row, mi_col);
    for (i = 0; i < MAX_MB_PLANE; i++) {
      dst_buf[i] = xd->plane[i].dst.buf;
      dst_stride[i] = xd->plane[i].dst.stride;
    }
#if CONFIG_VP9_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
      dec_predict_sb_complex_highbd(cm, xd, tile, mi_row, mi_col, mi_row,
                                    mi_col, bsize, bsize, dst_buf, dst_stride);
    else
#endif  // CONFIG_VP9_HIGHBITDEPTH
      dec_predict_sb_complex(cm, xd, tile, mi_row, mi_col, mi_row, mi_col,
                             bsize, bsize, dst_buf, dst_stride);

    if (!skip) {
      int eobtotal = 0;
      struct inter_args arg = { cm, xd, r, &eobtotal };
      set_offsets_topblock(cm, xd, tile, bsize, mi_row, mi_col);
#if CONFIG_EXT_TX
      xd->mi[0].mbmi.ext_txfrm = txfm;
#endif
      vp9_foreach_transformed_block(xd, bsize, reconstruct_inter_block, &arg);
      if (!(subsize < BLOCK_8X8) && eobtotal == 0)
        skip = 1;
    }
    set_param_topblock(cm, xd, bsize, mi_row, mi_col,
#if CONFIG_EXT_TX
                       txfm,
#endif
                       skip);
  }
#endif  // CONFIG_SUPERTX

  // update partition context
  if (bsize >= BLOCK_8X8 &&
      (bsize == BLOCK_8X8 || partition != PARTITION_SPLIT))
    update_partition_context(xd, mi_row, mi_col, subsize, bsize);
}

static void setup_token_decoder(const uint8_t *data,
                                const uint8_t *data_end,
                                size_t read_size,
                                struct vpx_internal_error_info *error_info,
                                vp9_reader *r,
                                vpx_decrypt_cb decrypt_cb,
                                void *decrypt_state) {
  // Validate the calculated partition length. If the buffer
  // described by the partition can't be fully read, then restrict
  // it to the portion that can be (for EC mode) or throw an error.
  if (!read_is_valid(data, read_size, data_end))
    vpx_internal_error(error_info, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt tile length");

  if (vp9_reader_init(r, data, read_size, decrypt_cb, decrypt_state))
    vpx_internal_error(error_info, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder %d", 1);
}

static void read_coef_probs_common(vp9_coeff_probs_model *coef_probs,
                                   vp9_reader *r) {
  int i, j, k, l, m;

  if (vp9_read_bit(r))
    for (i = 0; i < PLANE_TYPES; ++i)
      for (j = 0; j < REF_TYPES; ++j)
        for (k = 0; k < COEF_BANDS; ++k)
          for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l)
            for (m = 0; m < UNCONSTRAINED_NODES; ++m)
              vp9_diff_update_prob(r, &coef_probs[i][j][k][l][m]);
}

static void read_coef_probs(FRAME_CONTEXT *fc, TX_MODE tx_mode,
                            vp9_reader *r) {
    const TX_SIZE max_tx_size = tx_mode_to_biggest_tx_size[tx_mode];
    TX_SIZE tx_size;
    for (tx_size = TX_4X4; tx_size <= max_tx_size; ++tx_size)
      read_coef_probs_common(fc->coef_probs[tx_size], r);
}

static void setup_segmentation(struct segmentation *seg,
                               struct vp9_read_bit_buffer *rb) {
  int i, j;

  seg->update_map = 0;
  seg->update_data = 0;

  seg->enabled = vp9_rb_read_bit(rb);
  if (!seg->enabled)
    return;

  // Segmentation map update
  seg->update_map = vp9_rb_read_bit(rb);
  if (seg->update_map) {
    for (i = 0; i < SEG_TREE_PROBS; i++)
      seg->tree_probs[i] = vp9_rb_read_bit(rb) ? vp9_rb_read_literal(rb, 8)
                                               : MAX_PROB;

    seg->temporal_update = vp9_rb_read_bit(rb);
    if (seg->temporal_update) {
      for (i = 0; i < PREDICTION_PROBS; i++)
        seg->pred_probs[i] = vp9_rb_read_bit(rb) ? vp9_rb_read_literal(rb, 8)
                                                 : MAX_PROB;
    } else {
      for (i = 0; i < PREDICTION_PROBS; i++)
        seg->pred_probs[i] = MAX_PROB;
    }
  }

  // Segmentation data update
  seg->update_data = vp9_rb_read_bit(rb);
  if (seg->update_data) {
    seg->abs_delta = vp9_rb_read_bit(rb);

    vp9_clearall_segfeatures(seg);

    for (i = 0; i < MAX_SEGMENTS; i++) {
      for (j = 0; j < SEG_LVL_MAX; j++) {
        int data = 0;
        const int feature_enabled = vp9_rb_read_bit(rb);
        if (feature_enabled) {
          vp9_enable_segfeature(seg, i, j);
          data = decode_unsigned_max(rb, vp9_seg_feature_data_max(j));
          if (vp9_is_segfeature_signed(j))
            data = vp9_rb_read_bit(rb) ? -data : data;
        }
        vp9_set_segdata(seg, i, j, data);
      }
    }
  }
}

static void setup_loopfilter(struct loopfilter *lf,
                             struct vp9_read_bit_buffer *rb) {
  lf->filter_level = vp9_rb_read_literal(rb, 6);
  lf->sharpness_level = vp9_rb_read_literal(rb, 3);

  // Read in loop filter deltas applied at the MB level based on mode or ref
  // frame.
  lf->mode_ref_delta_update = 0;

  lf->mode_ref_delta_enabled = vp9_rb_read_bit(rb);
  if (lf->mode_ref_delta_enabled) {
    lf->mode_ref_delta_update = vp9_rb_read_bit(rb);
    if (lf->mode_ref_delta_update) {
      int i;

      for (i = 0; i < MAX_REF_LF_DELTAS; i++)
        if (vp9_rb_read_bit(rb))
          lf->ref_deltas[i] = vp9_rb_read_signed_literal(rb, 6);

      for (i = 0; i < MAX_MODE_LF_DELTAS; i++)
        if (vp9_rb_read_bit(rb))
          lf->mode_deltas[i] = vp9_rb_read_signed_literal(rb, 6);
    }
  }
}

static int read_delta_q(struct vp9_read_bit_buffer *rb, int *delta_q) {
  const int old = *delta_q;
  *delta_q = vp9_rb_read_bit(rb) ? vp9_rb_read_signed_literal(rb, 4) : 0;
  return old != *delta_q;
}

static void setup_quantization(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                               struct vp9_read_bit_buffer *rb) {
  int update = 0;

  cm->base_qindex = vp9_rb_read_literal(rb, QINDEX_BITS);
  update |= read_delta_q(rb, &cm->y_dc_delta_q);
  update |= read_delta_q(rb, &cm->uv_dc_delta_q);
  update |= read_delta_q(rb, &cm->uv_ac_delta_q);
  if (update || cm->bit_depth != cm->dequant_bit_depth) {
    vp9_init_dequantizer(cm);
    cm->dequant_bit_depth = cm->bit_depth;
  }

  xd->lossless = cm->base_qindex == 0 &&
                 cm->y_dc_delta_q == 0 &&
                 cm->uv_dc_delta_q == 0 &&
                 cm->uv_ac_delta_q == 0;
#if CONFIG_VP9_HIGHBITDEPTH
  xd->bd = (int)cm->bit_depth;
#endif
}

static INTERP_FILTER read_interp_filter(struct vp9_read_bit_buffer *rb) {
  const INTERP_FILTER literal_to_filter[] = { EIGHTTAP_SMOOTH,
                                              EIGHTTAP,
                                              EIGHTTAP_SHARP,
                                              BILINEAR };
  return vp9_rb_read_bit(rb) ? SWITCHABLE
                             : literal_to_filter[vp9_rb_read_literal(rb, 2)];
}

void vp9_read_frame_size(struct vp9_read_bit_buffer *rb,
                         int *width, int *height) {
  *width = vp9_rb_read_literal(rb, 16) + 1;
  *height = vp9_rb_read_literal(rb, 16) + 1;
}

static void setup_display_size(VP9_COMMON *cm, struct vp9_read_bit_buffer *rb) {
  cm->display_width = cm->width;
  cm->display_height = cm->height;
  if (vp9_rb_read_bit(rb))
    vp9_read_frame_size(rb, &cm->display_width, &cm->display_height);
}

static void resize_context_buffers(VP9_COMMON *cm, int width, int height) {
#if CONFIG_SIZE_LIMIT
  if (width > DECODE_WIDTH_LIMIT || height > DECODE_HEIGHT_LIMIT)
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Width and height beyond allowed size.");
#endif
  if (cm->width != width || cm->height != height) {
    const int new_mi_rows =
        ALIGN_POWER_OF_TWO(height, MI_SIZE_LOG2) >> MI_SIZE_LOG2;
    const int new_mi_cols =
        ALIGN_POWER_OF_TWO(width,  MI_SIZE_LOG2) >> MI_SIZE_LOG2;

    // Allocations in vp9_alloc_context_buffers() depend on individual
    // dimensions as well as the overall size.
    if (new_mi_cols > cm->mi_cols || new_mi_rows > cm->mi_rows) {
      if (vp9_alloc_context_buffers(cm, width, height))
        vpx_internal_error(&cm->error, VPX_CODEC_MEM_ERROR,
                           "Failed to allocate context buffers");
    } else {
      vp9_set_mb_mi(cm, width, height);
    }
    vp9_init_context_buffers(cm);
    cm->width = width;
    cm->height = height;
  }
}

static void setup_frame_size(VP9_COMMON *cm, struct vp9_read_bit_buffer *rb) {
  int width, height;
  vp9_read_frame_size(rb, &width, &height);
  resize_context_buffers(cm, width, height);
  setup_display_size(cm, rb);

  if (vp9_realloc_frame_buffer(
      get_frame_new_buffer(cm), cm->width, cm->height,
      cm->subsampling_x, cm->subsampling_y,
#if CONFIG_VP9_HIGHBITDEPTH
      cm->use_highbitdepth,
#endif
      VP9_DEC_BORDER_IN_PIXELS,
      &cm->frame_bufs[cm->new_fb_idx].raw_frame_buffer, cm->get_fb_cb,
      cm->cb_priv)) {
    vpx_internal_error(&cm->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate frame buffer");
  }
  cm->frame_bufs[cm->new_fb_idx].buf.subsampling_x = cm->subsampling_x;
  cm->frame_bufs[cm->new_fb_idx].buf.subsampling_y = cm->subsampling_y;
  cm->frame_bufs[cm->new_fb_idx].buf.color_space =
      (vpx_color_space_t)cm->color_space;
  cm->frame_bufs[cm->new_fb_idx].buf.bit_depth = (unsigned int)cm->bit_depth;
}

static INLINE int valid_ref_frame_img_fmt(vpx_bit_depth_t ref_bit_depth,
                                          int ref_xss, int ref_yss,
                                          vpx_bit_depth_t this_bit_depth,
                                          int this_xss, int this_yss) {
  return ref_bit_depth == this_bit_depth && ref_xss == this_xss &&
         ref_yss == this_yss;
}

static void setup_frame_size_with_refs(VP9_COMMON *cm,
                                       struct vp9_read_bit_buffer *rb) {
  int width, height;
  int found = 0, i;
  int has_valid_ref_frame = 0;
  for (i = 0; i < REFS_PER_FRAME; ++i) {
    if (vp9_rb_read_bit(rb)) {
      YV12_BUFFER_CONFIG *const buf = cm->frame_refs[i].buf;
      width = buf->y_crop_width;
      height = buf->y_crop_height;
      if (buf->corrupted) {
        vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                           "Frame reference is corrupt");
      }
      found = 1;
      break;
    }
  }

  if (!found)
    vp9_read_frame_size(rb, &width, &height);

  if (width <= 0 || height <= 0)
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Invalid frame size");

  // Check to make sure at least one of frames that this frame references
  // has valid dimensions.
  for (i = 0; i < REFS_PER_FRAME; ++i) {
    RefBuffer *const ref_frame = &cm->frame_refs[i];
    has_valid_ref_frame |= valid_ref_frame_size(ref_frame->buf->y_crop_width,
                                                ref_frame->buf->y_crop_height,
                                                width, height);
  }
  if (!has_valid_ref_frame)
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Referenced frame has invalid size");
  for (i = 0; i < REFS_PER_FRAME; ++i) {
    RefBuffer *const ref_frame = &cm->frame_refs[i];
    if (!valid_ref_frame_img_fmt(
            ref_frame->buf->bit_depth,
            ref_frame->buf->subsampling_x,
            ref_frame->buf->subsampling_y,
            cm->bit_depth,
            cm->subsampling_x,
            cm->subsampling_y))
      vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                         "Referenced frame has incompatible color space");
  }

  resize_context_buffers(cm, width, height);
  setup_display_size(cm, rb);

  if (vp9_realloc_frame_buffer(
      get_frame_new_buffer(cm), cm->width, cm->height,
      cm->subsampling_x, cm->subsampling_y,
#if CONFIG_VP9_HIGHBITDEPTH
      cm->use_highbitdepth,
#endif
      VP9_DEC_BORDER_IN_PIXELS,
      &cm->frame_bufs[cm->new_fb_idx].raw_frame_buffer, cm->get_fb_cb,
      cm->cb_priv)) {
    vpx_internal_error(&cm->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate frame buffer");
  }
  cm->frame_bufs[cm->new_fb_idx].buf.subsampling_x = cm->subsampling_x;
  cm->frame_bufs[cm->new_fb_idx].buf.subsampling_y = cm->subsampling_y;
  cm->frame_bufs[cm->new_fb_idx].buf.bit_depth = (unsigned int)cm->bit_depth;
}

static void setup_tile_info(VP9_COMMON *cm, struct vp9_read_bit_buffer *rb) {
  int min_log2_tile_cols, max_log2_tile_cols, max_ones;
  vp9_get_tile_n_bits(cm->mi_cols, &min_log2_tile_cols, &max_log2_tile_cols);

  // columns
  max_ones = max_log2_tile_cols - min_log2_tile_cols;
  cm->log2_tile_cols = min_log2_tile_cols;
  while (max_ones-- && vp9_rb_read_bit(rb))
    cm->log2_tile_cols++;

  if (cm->log2_tile_cols > 6)
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Invalid number of tile columns");

  // rows
  cm->log2_tile_rows = vp9_rb_read_bit(rb);
  if (cm->log2_tile_rows)
    cm->log2_tile_rows += vp9_rb_read_bit(rb);
}

typedef struct TileBuffer {
  const uint8_t *data;
  size_t size;
  int col;  // only used with multi-threaded decoding
} TileBuffer;

// Reads the next tile returning its size and adjusting '*data' accordingly
// based on 'is_last'.
static void get_tile_buffer(const uint8_t *const data_end,
                            int is_last,
                            struct vpx_internal_error_info *error_info,
                            const uint8_t **data,
                            vpx_decrypt_cb decrypt_cb, void *decrypt_state,
                            TileBuffer *buf) {
  size_t size;

  if (!is_last) {
    if (!read_is_valid(*data, 4, data_end))
      vpx_internal_error(error_info, VPX_CODEC_CORRUPT_FRAME,
                         "Truncated packet or corrupt tile length");

    if (decrypt_cb) {
      uint8_t be_data[4];
      decrypt_cb(decrypt_state, *data, be_data, 4);
      size = mem_get_be32(be_data);
    } else {
      size = mem_get_be32(*data);
    }
    *data += 4;

    if (size > (size_t)(data_end - *data))
      vpx_internal_error(error_info, VPX_CODEC_CORRUPT_FRAME,
                         "Truncated packet or corrupt tile size");
  } else {
    size = data_end - *data;
  }

  buf->data = *data;
  buf->size = size;

  *data += size;
}

static void get_tile_buffers(VP9Decoder *pbi,
                             const uint8_t *data, const uint8_t *data_end,
                             int tile_cols, int tile_rows,
                             TileBuffer (*tile_buffers)[1 << 6]) {
  int r, c;

  for (r = 0; r < tile_rows; ++r) {
    for (c = 0; c < tile_cols; ++c) {
      const int is_last = (r == tile_rows - 1) && (c == tile_cols - 1);
      TileBuffer *const buf = &tile_buffers[r][c];
      buf->col = c;
      get_tile_buffer(data_end, is_last, &pbi->common.error, &data,
                      pbi->decrypt_cb, pbi->decrypt_state, buf);
    }
  }
}

static const uint8_t *decode_tiles(VP9Decoder *pbi,
                                   const uint8_t *data,
                                   const uint8_t *data_end) {
  VP9_COMMON *const cm = &pbi->common;
  const VP9WorkerInterface *const winterface = vp9_get_worker_interface();
  const int aligned_cols = mi_cols_aligned_to_sb(cm->mi_cols);
  const int tile_cols = 1 << cm->log2_tile_cols;
  const int tile_rows = 1 << cm->log2_tile_rows;
  TileBuffer tile_buffers[4][1 << 6];
  int tile_row, tile_col;
  int mi_row, mi_col;
  TileData *tile_data = NULL;

  if (cm->lf.filter_level && pbi->lf_worker.data1 == NULL) {
    CHECK_MEM_ERROR(cm, pbi->lf_worker.data1,
                    vpx_memalign(32, sizeof(LFWorkerData)));
    pbi->lf_worker.hook = (VP9WorkerHook)vp9_loop_filter_worker;
    if (pbi->max_threads > 1 && !winterface->reset(&pbi->lf_worker)) {
      vpx_internal_error(&cm->error, VPX_CODEC_ERROR,
                         "Loop filter thread creation failed");
    }
  }

  if (cm->lf.filter_level) {
    LFWorkerData *const lf_data = (LFWorkerData*)pbi->lf_worker.data1;
    // Be sure to sync as we might be resuming after a failed frame decode.
    winterface->sync(&pbi->lf_worker);
    lf_data->frame_buffer = get_frame_new_buffer(cm);
    lf_data->cm = cm;
    vp9_copy(lf_data->planes, pbi->mb.plane);
    lf_data->stop = 0;
    lf_data->y_only = 0;
    vp9_loop_filter_frame_init(cm, cm->lf.filter_level);
  }

  assert(tile_rows <= 4);
  assert(tile_cols <= (1 << 6));

  // Note: this memset assumes above_context[0], [1] and [2]
  // are allocated as part of the same buffer.
  vpx_memset(cm->above_context, 0,
             sizeof(*cm->above_context) * MAX_MB_PLANE * 2 * aligned_cols);

  vpx_memset(cm->above_seg_context, 0,
             sizeof(*cm->above_seg_context) * aligned_cols);

  get_tile_buffers(pbi, data, data_end, tile_cols, tile_rows, tile_buffers);

  if (pbi->tile_data == NULL ||
      (tile_cols * tile_rows) != pbi->total_tiles) {
    vpx_free(pbi->tile_data);
    CHECK_MEM_ERROR(
        cm,
        pbi->tile_data,
        vpx_memalign(32, tile_cols * tile_rows * (sizeof(*pbi->tile_data))));
    pbi->total_tiles = tile_rows * tile_cols;
  }

  // Load all tile information into tile_data.
  for (tile_row = 0; tile_row < tile_rows; ++tile_row) {
    for (tile_col = 0; tile_col < tile_cols; ++tile_col) {
      TileInfo tile;
      const TileBuffer *const buf = &tile_buffers[tile_row][tile_col];
      tile_data = pbi->tile_data + tile_cols * tile_row + tile_col;
      tile_data->cm = cm;
      tile_data->xd = pbi->mb;
      tile_data->xd.corrupted = 0;
      vp9_tile_init(&tile, tile_data->cm, tile_row, tile_col);

      setup_token_decoder(buf->data, data_end, buf->size, &cm->error,
                          &tile_data->bit_reader, pbi->decrypt_cb,
                          pbi->decrypt_state);

      init_macroblockd(cm, &tile_data->xd);

      vp9_zero(tile_data->xd.dqcoeff);
    }
  }

  for (tile_row = 0; tile_row < tile_rows; ++tile_row) {
    TileInfo tile;
    vp9_tile_set_row(&tile, cm, tile_row);
    for (mi_row = tile.mi_row_start; mi_row < tile.mi_row_end;
         mi_row += MI_BLOCK_SIZE) {
      for (tile_col = 0; tile_col < tile_cols; ++tile_col) {
        const int col = pbi->inv_tile_order ?
                        tile_cols - tile_col - 1 : tile_col;
        tile_data = pbi->tile_data + tile_cols * tile_row + col;
        vp9_tile_set_col(&tile, tile_data->cm, col);
        vp9_zero(tile_data->xd.left_context);
        vp9_zero(tile_data->xd.left_seg_context);
        for (mi_col = tile.mi_col_start; mi_col < tile.mi_col_end;
             mi_col += MI_BLOCK_SIZE) {
          decode_partition(tile_data->cm, &tile_data->xd, &tile,
#if CONFIG_SUPERTX
                           0,
#endif
                           mi_row, mi_col,
                           &tile_data->bit_reader, BLOCK_64X64);
        }
        pbi->mb.corrupted |= tile_data->xd.corrupted;
      }
      // Loopfilter one row.
      if (cm->lf.filter_level && !pbi->mb.corrupted) {
        const int lf_start = mi_row - MI_BLOCK_SIZE;
        LFWorkerData *const lf_data = (LFWorkerData*)pbi->lf_worker.data1;

        // delay the loopfilter by 1 macroblock row.
        if (lf_start < 0) continue;

        // decoding has completed: finish up the loop filter in this thread.
        if (mi_row + MI_BLOCK_SIZE >= cm->mi_rows) continue;

        winterface->sync(&pbi->lf_worker);
        lf_data->start = lf_start;
        lf_data->stop = mi_row;
        if (pbi->max_threads > 1) {
          winterface->launch(&pbi->lf_worker);
        } else {
          winterface->execute(&pbi->lf_worker);
        }
      }
    }
  }

  // Loopfilter remaining rows in the frame.
  if (cm->lf.filter_level && !pbi->mb.corrupted) {
    LFWorkerData *const lf_data = (LFWorkerData*)pbi->lf_worker.data1;
    winterface->sync(&pbi->lf_worker);
    lf_data->start = lf_data->stop;
    lf_data->stop = cm->mi_rows;
    winterface->execute(&pbi->lf_worker);
  }

  // Get last tile data.
  tile_data = pbi->tile_data + tile_cols * tile_rows - 1;

  return vp9_reader_find_end(&tile_data->bit_reader);
}

static int tile_worker_hook(TileWorkerData *const tile_data,
                            const TileInfo *const tile) {
  int mi_row, mi_col;

  for (mi_row = tile->mi_row_start; mi_row < tile->mi_row_end;
       mi_row += MI_BLOCK_SIZE) {
    vp9_zero(tile_data->xd.left_context);
    vp9_zero(tile_data->xd.left_seg_context);
    for (mi_col = tile->mi_col_start; mi_col < tile->mi_col_end;
         mi_col += MI_BLOCK_SIZE) {
      decode_partition(tile_data->cm, &tile_data->xd, tile,
#if CONFIG_SUPERTX
                       0,
#endif
                       mi_row, mi_col, &tile_data->bit_reader, BLOCK_64X64);
    }
  }
  return !tile_data->xd.corrupted;
}

// sorts in descending order
static int compare_tile_buffers(const void *a, const void *b) {
  const TileBuffer *const buf1 = (const TileBuffer*)a;
  const TileBuffer *const buf2 = (const TileBuffer*)b;
  if (buf1->size < buf2->size) {
    return 1;
  } else if (buf1->size == buf2->size) {
    return 0;
  } else {
    return -1;
  }
}

static const uint8_t *decode_tiles_mt(VP9Decoder *pbi,
                                      const uint8_t *data,
                                      const uint8_t *data_end) {
  VP9_COMMON *const cm = &pbi->common;
  const VP9WorkerInterface *const winterface = vp9_get_worker_interface();
  const uint8_t *bit_reader_end = NULL;
  const int aligned_mi_cols = mi_cols_aligned_to_sb(cm->mi_cols);
  const int tile_cols = 1 << cm->log2_tile_cols;
  const int tile_rows = 1 << cm->log2_tile_rows;
  const int num_workers = MIN(pbi->max_threads & ~1, tile_cols);
  TileBuffer tile_buffers[1][1 << 6];
  int n;
  int final_worker = -1;

  assert(tile_cols <= (1 << 6));
  assert(tile_rows == 1);
  (void)tile_rows;

  // TODO(jzern): See if we can remove the restriction of passing in max
  // threads to the decoder.
  if (pbi->num_tile_workers == 0) {
    const int num_threads = pbi->max_threads & ~1;
    int i;
    // TODO(jzern): Allocate one less worker, as in the current code we only
    // use num_threads - 1 workers.
    CHECK_MEM_ERROR(cm, pbi->tile_workers,
                    vpx_malloc(num_threads * sizeof(*pbi->tile_workers)));
    for (i = 0; i < num_threads; ++i) {
      VP9Worker *const worker = &pbi->tile_workers[i];
      ++pbi->num_tile_workers;

      winterface->init(worker);
      CHECK_MEM_ERROR(cm, worker->data1,
                      vpx_memalign(32, sizeof(TileWorkerData)));
      CHECK_MEM_ERROR(cm, worker->data2, vpx_malloc(sizeof(TileInfo)));
      if (i < num_threads - 1 && !winterface->reset(worker)) {
        vpx_internal_error(&cm->error, VPX_CODEC_ERROR,
                           "Tile decoder thread creation failed");
      }
    }
  }

  // Reset tile decoding hook
  for (n = 0; n < num_workers; ++n) {
    winterface->sync(&pbi->tile_workers[n]);
    pbi->tile_workers[n].hook = (VP9WorkerHook)tile_worker_hook;
  }

  // Note: this memset assumes above_context[0], [1] and [2]
  // are allocated as part of the same buffer.
  vpx_memset(cm->above_context, 0,
             sizeof(*cm->above_context) * MAX_MB_PLANE * 2 * aligned_mi_cols);
  vpx_memset(cm->above_seg_context, 0,
             sizeof(*cm->above_seg_context) * aligned_mi_cols);

  // Load tile data into tile_buffers
  get_tile_buffers(pbi, data, data_end, tile_cols, tile_rows, tile_buffers);

  // Sort the buffers based on size in descending order.
  qsort(tile_buffers[0], tile_cols, sizeof(tile_buffers[0][0]),
        compare_tile_buffers);

  // Rearrange the tile buffers such that per-tile group the largest, and
  // presumably the most difficult, tile will be decoded in the main thread.
  // This should help minimize the number of instances where the main thread is
  // waiting for a worker to complete.
  {
    int group_start = 0;
    while (group_start < tile_cols) {
      const TileBuffer largest = tile_buffers[0][group_start];
      const int group_end = MIN(group_start + num_workers, tile_cols) - 1;
      memmove(tile_buffers[0] + group_start, tile_buffers[0] + group_start + 1,
              (group_end - group_start) * sizeof(tile_buffers[0][0]));
      tile_buffers[0][group_end] = largest;
      group_start = group_end + 1;
    }
  }

  n = 0;
  while (n < tile_cols) {
    int i;
    for (i = 0; i < num_workers && n < tile_cols; ++i) {
      VP9Worker *const worker = &pbi->tile_workers[i];
      TileWorkerData *const tile_data = (TileWorkerData*)worker->data1;
      TileInfo *const tile = (TileInfo*)worker->data2;
      TileBuffer *const buf = &tile_buffers[0][n];

      tile_data->cm = cm;
      tile_data->xd = pbi->mb;
      tile_data->xd.corrupted = 0;
      vp9_tile_init(tile, tile_data->cm, 0, buf->col);
      setup_token_decoder(buf->data, data_end, buf->size, &cm->error,
                          &tile_data->bit_reader, pbi->decrypt_cb,
                          pbi->decrypt_state);
      init_macroblockd(cm, &tile_data->xd);
      vp9_zero(tile_data->xd.dqcoeff);

      worker->had_error = 0;
      if (i == num_workers - 1 || n == tile_cols - 1) {
        winterface->execute(worker);
      } else {
        winterface->launch(worker);
      }

      if (buf->col == tile_cols - 1) {
        final_worker = i;
      }

      ++n;
    }

    for (; i > 0; --i) {
      VP9Worker *const worker = &pbi->tile_workers[i - 1];
      pbi->mb.corrupted |= !winterface->sync(worker);
    }
    if (final_worker > -1) {
      TileWorkerData *const tile_data =
          (TileWorkerData*)pbi->tile_workers[final_worker].data1;
      bit_reader_end = vp9_reader_find_end(&tile_data->bit_reader);
      final_worker = -1;
    }
  }

  return bit_reader_end;
}

static void error_handler(void *data) {
  VP9_COMMON *const cm = (VP9_COMMON *)data;
  vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME, "Truncated packet");
}

int vp9_read_sync_code(struct vp9_read_bit_buffer *const rb) {
  return vp9_rb_read_literal(rb, 8) == VP9_SYNC_CODE_0 &&
         vp9_rb_read_literal(rb, 8) == VP9_SYNC_CODE_1 &&
         vp9_rb_read_literal(rb, 8) == VP9_SYNC_CODE_2;
}

BITSTREAM_PROFILE vp9_read_profile(struct vp9_read_bit_buffer *rb) {
  int profile = vp9_rb_read_bit(rb);
  profile |= vp9_rb_read_bit(rb) << 1;
  if (profile > 2)
    profile += vp9_rb_read_bit(rb);
  return (BITSTREAM_PROFILE) profile;
}

static void read_bitdepth_colorspace_sampling(
    VP9_COMMON *cm, struct vp9_read_bit_buffer *rb) {
  if (cm->profile >= PROFILE_2) {
    cm->bit_depth = vp9_rb_read_bit(rb) ? VPX_BITS_12 : VPX_BITS_10;
#if CONFIG_VP9_HIGHBITDEPTH
    cm->use_highbitdepth = 1;
#endif
  } else {
    cm->bit_depth = VPX_BITS_8;
#if CONFIG_VP9_HIGHBITDEPTH
    cm->use_highbitdepth = 0;
#endif
  }
  cm->color_space = vp9_rb_read_literal(rb, 3);
  if (cm->color_space != VPX_CS_SRGB) {
    vp9_rb_read_bit(rb);  // [16,235] (including xvycc) vs [0,255] range
    if (cm->profile == PROFILE_1 || cm->profile == PROFILE_3) {
      cm->subsampling_x = vp9_rb_read_bit(rb);
      cm->subsampling_y = vp9_rb_read_bit(rb);
      if (cm->subsampling_x == 1 && cm->subsampling_y == 1)
        vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                           "4:2:0 color not supported in profile 1 or 3");
      if (vp9_rb_read_bit(rb))
        vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                           "Reserved bit set");
    } else {
      cm->subsampling_y = cm->subsampling_x = 1;
    }
  } else {
    if (cm->profile == PROFILE_1 || cm->profile == PROFILE_3) {
      // Note if colorspace is SRGB then 4:4:4 chroma sampling is assumed.
      // 4:2:2 or 4:4:0 chroma sampling is not allowed.
      cm->subsampling_y = cm->subsampling_x = 0;
      if (vp9_rb_read_bit(rb))
        vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                           "Reserved bit set");
    } else {
      vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                         "4:4:4 color not supported in profile 0 or 2");
    }
  }
}

static size_t read_uncompressed_header(VP9Decoder *pbi,
                                       struct vp9_read_bit_buffer *rb) {
  VP9_COMMON *const cm = &pbi->common;
  size_t sz;
  int i;

  cm->last_frame_type = cm->frame_type;

  if (vp9_rb_read_literal(rb, 2) != VP9_FRAME_MARKER)
      vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                         "Invalid frame marker");

  cm->profile = vp9_read_profile(rb);

  if (cm->profile >= MAX_PROFILES)
    vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                       "Unsupported bitstream profile");

  cm->show_existing_frame = vp9_rb_read_bit(rb);
  if (cm->show_existing_frame) {
    // Show an existing frame directly.
    const int frame_to_show = cm->ref_frame_map[vp9_rb_read_literal(rb, 3)];

    if (frame_to_show < 0 || cm->frame_bufs[frame_to_show].ref_count < 1)
      vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                         "Buffer %d does not contain a decoded frame",
                         frame_to_show);

    ref_cnt_fb(cm->frame_bufs, &cm->new_fb_idx, frame_to_show);
    pbi->refresh_frame_flags = 0;
    cm->lf.filter_level = 0;
    cm->show_frame = 1;
    return 0;
  }

  cm->frame_type = (FRAME_TYPE) vp9_rb_read_bit(rb);
  cm->show_frame = vp9_rb_read_bit(rb);
  cm->error_resilient_mode = vp9_rb_read_bit(rb);

  if (cm->frame_type == KEY_FRAME) {
    if (!vp9_read_sync_code(rb))
      vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                         "Invalid frame sync code");

    read_bitdepth_colorspace_sampling(cm, rb);
    pbi->refresh_frame_flags = (1 << REF_FRAMES) - 1;

    for (i = 0; i < REFS_PER_FRAME; ++i) {
      cm->frame_refs[i].idx = -1;
      cm->frame_refs[i].buf = NULL;
    }

    setup_frame_size(cm, rb);
    pbi->need_resync = 0;
  } else {
    cm->intra_only = cm->show_frame ? 0 : vp9_rb_read_bit(rb);

    cm->reset_frame_context = cm->error_resilient_mode ?
        0 : vp9_rb_read_literal(rb, 2);

    if (cm->intra_only) {
      if (!vp9_read_sync_code(rb))
        vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                           "Invalid frame sync code");
      if (cm->profile > PROFILE_0) {
        read_bitdepth_colorspace_sampling(cm, rb);
      } else {
        // NOTE: The intra-only frame header does not include the specification
        // of either the color format or color sub-sampling in profile 0. VP9
        // specifies that the default color space should be YUV 4:2:0 in this
        // case (normative).
        cm->color_space = VPX_CS_BT_601;
        cm->subsampling_y = cm->subsampling_x = 1;
        cm->bit_depth = VPX_BITS_8;
#if CONFIG_VP9_HIGHBITDEPTH
        cm->use_highbitdepth = 0;
#endif
      }

      pbi->refresh_frame_flags = vp9_rb_read_literal(rb, REF_FRAMES);
      setup_frame_size(cm, rb);
      pbi->need_resync = 0;
    } else {
      pbi->refresh_frame_flags = vp9_rb_read_literal(rb, REF_FRAMES);
      for (i = 0; i < REFS_PER_FRAME; ++i) {
        const int ref = vp9_rb_read_literal(rb, REF_FRAMES_LOG2);
        const int idx = cm->ref_frame_map[ref];
        RefBuffer *const ref_frame = &cm->frame_refs[i];
        ref_frame->idx = idx;
        ref_frame->buf = &cm->frame_bufs[idx].buf;
        cm->ref_frame_sign_bias[LAST_FRAME + i] = vp9_rb_read_bit(rb);
      }

      setup_frame_size_with_refs(cm, rb);

      cm->allow_high_precision_mv = vp9_rb_read_bit(rb);
      cm->interp_filter = read_interp_filter(rb);

      for (i = 0; i < REFS_PER_FRAME; ++i) {
        RefBuffer *const ref_buf = &cm->frame_refs[i];
#if CONFIG_VP9_HIGHBITDEPTH
        vp9_setup_scale_factors_for_frame(&ref_buf->sf,
                                          ref_buf->buf->y_crop_width,
                                          ref_buf->buf->y_crop_height,
                                          cm->width, cm->height,
                                          cm->use_highbitdepth);
#else
        vp9_setup_scale_factors_for_frame(&ref_buf->sf,
                                          ref_buf->buf->y_crop_width,
                                          ref_buf->buf->y_crop_height,
                                          cm->width, cm->height);
#endif
        if (vp9_is_scaled(&ref_buf->sf))
          vp9_extend_frame_borders(ref_buf->buf);
      }
    }
  }
#if CONFIG_VP9_HIGHBITDEPTH
  get_frame_new_buffer(cm)->bit_depth = cm->bit_depth;
#endif

  if (pbi->need_resync) {
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Keyframe / intra-only frame required to reset decoder"
                       " state");
  }

  if (!cm->error_resilient_mode) {
    cm->refresh_frame_context = vp9_rb_read_bit(rb);
    cm->frame_parallel_decoding_mode = vp9_rb_read_bit(rb);
  } else {
    cm->refresh_frame_context = 0;
    cm->frame_parallel_decoding_mode = 1;
  }

  // This flag will be overridden by the call to vp9_setup_past_independence
  // below, forcing the use of context 0 for those frame types.
  cm->frame_context_idx = vp9_rb_read_literal(rb, FRAME_CONTEXTS_LOG2);

  if (frame_is_intra_only(cm) || cm->error_resilient_mode)
    vp9_setup_past_independence(cm);

  setup_loopfilter(&cm->lf, rb);
  setup_quantization(cm, &pbi->mb, rb);
  setup_segmentation(&cm->seg, rb);

  setup_tile_info(cm, rb);
  sz = vp9_rb_read_literal(rb, 16);

  if (sz == 0)
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Invalid header size");

  return sz;
}

#if CONFIG_EXT_TX
static void read_ext_tx_probs(FRAME_CONTEXT *fc, vp9_reader *r) {
  int i, j;
  if (vp9_read(r, GROUP_DIFF_UPDATE_PROB)) {
    for (j = TX_4X4; j <= TX_16X16; ++j)
      for (i = 0; i < EXT_TX_TYPES - 1; ++i)
        vp9_diff_update_prob(r, &fc->ext_tx_prob[j][i]);
  }
}
#endif  // CONFIG_EXT_TX

#if CONFIG_SUPERTX
static void read_supertx_probs(FRAME_CONTEXT *fc, vp9_reader *r) {
  int i, j;
  if (vp9_read(r, GROUP_DIFF_UPDATE_PROB)) {
    for (i = 0; i < PARTITION_SUPERTX_CONTEXTS; ++i) {
      for (j = 1; j < TX_SIZES; ++j) {
        vp9_diff_update_prob(r, &fc->supertx_prob[i][j]);
      }
    }
  }
}
#endif  // CONFIG_SUPERTX

#if CONFIG_COMPOUND_MODES
static void read_inter_compound_mode_probs(FRAME_CONTEXT *fc, vp9_reader *r) {
  int i, j;
  if (vp9_read(r, GROUP_DIFF_UPDATE_PROB)) {
    for (j = 0; j < INTER_MODE_CONTEXTS; ++j)
      for (i = 0; i < INTER_COMPOUND_MODES - 1; ++i)
        vp9_diff_update_prob(r, &fc->inter_compound_mode_probs[j][i]);
  }
}
#endif  // CONFIG_COMPOUND_MODES

static int read_compressed_header(VP9Decoder *pbi, const uint8_t *data,
                                  size_t partition_size) {
  VP9_COMMON *const cm = &pbi->common;
#if !CONFIG_TX_SKIP || CONFIG_SUPERTX
  MACROBLOCKD *const xd = &pbi->mb;
#endif
  FRAME_CONTEXT *const fc = &cm->fc;
  vp9_reader r;
  int k;

  if (vp9_reader_init(&r, data, partition_size, pbi->decrypt_cb,
                      pbi->decrypt_state))
    vpx_internal_error(&cm->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder 0");

#if CONFIG_TX_SKIP
  cm->tx_mode = read_tx_mode(&r);
#else
  cm->tx_mode = xd->lossless ? ONLY_4X4 : read_tx_mode(&r);
#endif
  if (cm->tx_mode == TX_MODE_SELECT)
    read_tx_mode_probs(&fc->tx_probs, &r);
  read_coef_probs(fc, cm->tx_mode, &r);

  for (k = 0; k < SKIP_CONTEXTS; ++k)
    vp9_diff_update_prob(&r, &fc->skip_probs[k]);

  if (!frame_is_intra_only(cm)) {
    nmv_context *const nmvc = &fc->nmvc;
    int i, j;

    read_inter_mode_probs(fc, &r);
#if CONFIG_COMPOUND_MODES
    read_inter_compound_mode_probs(fc, &r);
#endif

    if (cm->interp_filter == SWITCHABLE)
      read_switchable_interp_probs(fc, &r);

    for (i = 0; i < INTRA_INTER_CONTEXTS; i++)
      vp9_diff_update_prob(&r, &fc->intra_inter_prob[i]);

    cm->reference_mode = read_frame_reference_mode(cm, &r);
    if (cm->reference_mode != SINGLE_REFERENCE)
      setup_compound_reference_mode(cm);
    read_frame_reference_mode_probs(cm, &r);

    for (j = 0; j < BLOCK_SIZE_GROUPS; j++)
      for (i = 0; i < INTRA_MODES - 1; ++i)
        vp9_diff_update_prob(&r, &fc->y_mode_prob[j][i]);

    for (j = 0; j < PARTITION_CONTEXTS; ++j)
      for (i = 0; i < PARTITION_TYPES - 1; ++i)
        vp9_diff_update_prob(&r, &fc->partition_prob[j][i]);

    read_mv_probs(nmvc, cm->allow_high_precision_mv, &r);
#if CONFIG_EXT_TX
    read_ext_tx_probs(fc, &r);
#endif
#if CONFIG_SUPERTX
    if (!xd->lossless)
      read_supertx_probs(fc, &r);
#endif
#if CONFIG_TX_SKIP
  for (i = 0; i < 2; i++)
    vp9_diff_update_prob(&r, &fc->y_tx_skip_prob[i]);
  for (i = 0; i < 2; i++)
    vp9_diff_update_prob(&r, &fc->uv_tx_skip_prob[i]);
#endif
#if CONFIG_COPY_MODE
    for (j = 0; j < COPY_MODE_CONTEXTS; j++) {
      for (i = 0; i < 1; i++)
        vp9_diff_update_prob(&r, &fc->copy_mode_probs_l2[j][i]);
      for (i = 0; i < COPY_MODE_COUNT - 2; i++)
        vp9_diff_update_prob(&r, &fc->copy_mode_probs[j][i]);
    }
#endif
#if CONFIG_INTERINTRA
    if (cm->reference_mode != COMPOUND_REFERENCE) {
      for (i = 0; i < BLOCK_SIZES; i++) {
        if (is_interintra_allowed(i)) {
          vp9_diff_update_prob(&r, &fc->interintra_prob[i]);
        }
      }
#if CONFIG_WEDGE_PARTITION
      for (i = 0; i < BLOCK_SIZES; i++) {
        if (is_interintra_allowed(i) && get_wedge_bits(i))
          vp9_diff_update_prob(&r, &fc->wedge_interintra_prob[i]);
      }
#endif  // CONFIG_WEDGE_PARTITION
    }
#endif  // CONFIG_INTERINTRA
#if CONFIG_WEDGE_PARTITION
    if (cm->reference_mode != SINGLE_REFERENCE) {
      for (i = 0; i < BLOCK_SIZES; i++) {
        if (get_wedge_bits(i))
          vp9_diff_update_prob(&r, &fc->wedge_interinter_prob[i]);
      }
    }
#endif  // CONFIG_WEDGE_PARTITION
  }
#if CONFIG_PALETTE
  if (frame_is_intra_only(cm))
    cm->allow_palette_mode = vp9_read_bit(&r);
#endif

  return vp9_reader_has_error(&r);
}

void vp9_init_dequantizer(VP9_COMMON *cm) {
  int q;

  for (q = 0; q < QINDEX_RANGE; q++) {
    int b;
    cm->y_dequant[q][0] = vp9_dc_quant(q, cm->y_dc_delta_q, cm->bit_depth);
    cm->y_dequant[q][1] = vp9_ac_quant(q, 0, cm->bit_depth);

    cm->uv_dequant[q][0] = vp9_dc_quant(q, cm->uv_dc_delta_q, cm->bit_depth);
    cm->uv_dequant[q][1] = vp9_ac_quant(q, cm->uv_ac_delta_q, cm->bit_depth);

#if CONFIG_NEW_QUANT
    for (b = 0; b < COEF_BANDS; ++b) {
      vp9_get_dequant_val_nuq(
          cm->y_dequant[q][b != 0], b, cm->bit_depth,
          cm->y_dequant_val_nuq[q][b], NULL);
      vp9_get_dequant_val_nuq(
          cm->uv_dequant[q][b != 0], b, cm->bit_depth,
          cm->uv_dequant_val_nuq[q][b], NULL);
    }
#endif  // CONFIG_NEW_QUANT
    (void) b;
  }
}

#ifdef NDEBUG
#define debug_check_frame_counts(cm) (void)0
#else  // !NDEBUG
// Counts should only be incremented when frame_parallel_decoding_mode and
// error_resilient_mode are disabled.
static void debug_check_frame_counts(const VP9_COMMON *const cm) {
  FRAME_COUNTS zero_counts;
  vp9_zero(zero_counts);
  assert(cm->frame_parallel_decoding_mode || cm->error_resilient_mode);
  assert(!memcmp(cm->counts.y_mode, zero_counts.y_mode,
                 sizeof(cm->counts.y_mode)));
  assert(!memcmp(cm->counts.uv_mode, zero_counts.uv_mode,
                 sizeof(cm->counts.uv_mode)));
  assert(!memcmp(cm->counts.partition, zero_counts.partition,
                 sizeof(cm->counts.partition)));
  assert(!memcmp(cm->counts.coef, zero_counts.coef,
                 sizeof(cm->counts.coef)));
  assert(!memcmp(cm->counts.eob_branch, zero_counts.eob_branch,
                 sizeof(cm->counts.eob_branch)));
  assert(!memcmp(cm->counts.switchable_interp, zero_counts.switchable_interp,
                 sizeof(cm->counts.switchable_interp)));
  assert(!memcmp(cm->counts.inter_mode, zero_counts.inter_mode,
                 sizeof(cm->counts.inter_mode)));
  assert(!memcmp(cm->counts.intra_inter, zero_counts.intra_inter,
                 sizeof(cm->counts.intra_inter)));
  assert(!memcmp(cm->counts.comp_inter, zero_counts.comp_inter,
                 sizeof(cm->counts.comp_inter)));
  assert(!memcmp(cm->counts.single_ref, zero_counts.single_ref,
                 sizeof(cm->counts.single_ref)));
  assert(!memcmp(cm->counts.comp_ref, zero_counts.comp_ref,
                 sizeof(cm->counts.comp_ref)));
  assert(!memcmp(&cm->counts.tx, &zero_counts.tx, sizeof(cm->counts.tx)));
  assert(!memcmp(cm->counts.skip, zero_counts.skip, sizeof(cm->counts.skip)));
  assert(!memcmp(&cm->counts.mv, &zero_counts.mv, sizeof(cm->counts.mv)));
#if CONFIG_EXT_TX
  assert(!memcmp(cm->counts.ext_tx, zero_counts.ext_tx,
                 sizeof(cm->counts.ext_tx)));
#endif
#if CONFIG_COMPOUND_MODES
  assert(!memcmp(cm->counts.inter_compound_mode,
                 zero_counts.inter_compound_mode,
                 sizeof(cm->counts.inter_compound_mode)));
#endif
}
#endif  // NDEBUG

static struct vp9_read_bit_buffer* init_read_bit_buffer(
    VP9Decoder *pbi,
    struct vp9_read_bit_buffer *rb,
    const uint8_t *data,
    const uint8_t *data_end,
    uint8_t *clear_data /* buffer size MAX_VP9_HEADER_SIZE */) {
  rb->bit_offset = 0;
  rb->error_handler = error_handler;
  rb->error_handler_data = &pbi->common;
  if (pbi->decrypt_cb) {
    const int n = (int)MIN(MAX_VP9_HEADER_SIZE, data_end - data);
    pbi->decrypt_cb(pbi->decrypt_state, data, clear_data, n);
    rb->bit_buffer = clear_data;
    rb->bit_buffer_end = clear_data + n;
  } else {
    rb->bit_buffer = data;
    rb->bit_buffer_end = data_end;
  }
  return rb;
}

void vp9_decode_frame(VP9Decoder *pbi,
                      const uint8_t *data, const uint8_t *data_end,
                      const uint8_t **p_data_end) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  struct vp9_read_bit_buffer rb = { NULL, NULL, 0, NULL, 0};

  uint8_t clear_data[MAX_VP9_HEADER_SIZE];
  const size_t first_partition_size = read_uncompressed_header(pbi,
      init_read_bit_buffer(pbi, &rb, data, data_end, clear_data));
  const int tile_rows = 1 << cm->log2_tile_rows;
  const int tile_cols = 1 << cm->log2_tile_cols;
  YV12_BUFFER_CONFIG *const new_fb = get_frame_new_buffer(cm);
  xd->cur_buf = new_fb;

  if (!first_partition_size) {
    // showing a frame directly
    *p_data_end = data + (cm->profile <= PROFILE_2 ? 1 : 2);
    return;
  }

  data += vp9_rb_bytes_read(&rb);
  if (!read_is_valid(data, first_partition_size, data_end))
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt header length");

  init_macroblockd(cm, &pbi->mb);

  if (!cm->error_resilient_mode)
    set_prev_mi(cm);
  else
    cm->prev_mi = NULL;

  setup_plane_dequants(cm, xd, cm->base_qindex);
  vp9_setup_block_planes(xd, cm->subsampling_x, cm->subsampling_y);

  cm->fc = cm->frame_contexts[cm->frame_context_idx];
  vp9_zero(cm->counts);
  vp9_zero(xd->dqcoeff);

  xd->corrupted = 0;
  new_fb->corrupted = read_compressed_header(pbi, data, first_partition_size);

  // TODO(jzern): remove frame_parallel_decoding_mode restriction for
  // single-frame tile decoding.
  if (pbi->max_threads > 1 && tile_rows == 1 && tile_cols > 1 &&
      cm->frame_parallel_decoding_mode) {
    *p_data_end = decode_tiles_mt(pbi, data + first_partition_size, data_end);
    if (!xd->corrupted) {
      // If multiple threads are used to decode tiles, then we use those threads
      // to do parallel loopfiltering.
      vp9_loop_filter_frame_mt(new_fb, pbi, cm, cm->lf.filter_level, 0);
    }
  } else {
    *p_data_end = decode_tiles(pbi, data + first_partition_size, data_end);
  }

  new_fb->corrupted |= xd->corrupted;
  if (!new_fb->corrupted) {
    if (!cm->error_resilient_mode && !cm->frame_parallel_decoding_mode) {
      vp9_adapt_coef_probs(cm);

      if (!frame_is_intra_only(cm)) {
        vp9_adapt_mode_probs(cm);
        vp9_adapt_mv_probs(cm, cm->allow_high_precision_mv);
      }
    } else {
      debug_check_frame_counts(cm);
    }
  } else {
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Decode failed. Frame data is corrupted.");
  }

  if (cm->refresh_frame_context)
    cm->frame_contexts[cm->frame_context_idx] = cm->fc;
}
