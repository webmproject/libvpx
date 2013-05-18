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

#include "vp9/decoder/vp9_onyxd_int.h"
#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_header.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/decoder/vp9_decodframe.h"
#include "vp9/decoder/vp9_detokenize.h"
#include "vp9/common/vp9_invtrans.h"
#include "vp9/common/vp9_alloccommon.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_quant_common.h"
#include "vpx_scale/vpx_scale.h"

#include "vp9/decoder/vp9_decodemv.h"
#include "vp9/common/vp9_extend.h"
#include "vp9/common/vp9_modecont.h"
#include "vpx_mem/vpx_mem.h"
#include "vp9/decoder/vp9_dboolhuff.h"

#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_tile_common.h"
#include "vp9_rtcd.h"

// #define DEC_DEBUG
#ifdef DEC_DEBUG
int dec_debug = 0;
#endif

static int read_le16(const uint8_t *p) {
  return (p[1] << 8) | p[0];
}

static int read_le32(const uint8_t *p) {
  return (p[3] << 24) | (p[2] << 16) | (p[1] << 8) | p[0];
}

// len == 0 is not allowed
static int read_is_valid(const uint8_t *start, size_t len,
                         const uint8_t *end) {
  return start + len > start && start + len <= end;
}

static void setup_txfm_mode(VP9_COMMON *pc, int lossless, vp9_reader *r) {
  if (lossless) {
    pc->txfm_mode = ONLY_4X4;
  } else {
    pc->txfm_mode = vp9_read_literal(r, 2);
    if (pc->txfm_mode == ALLOW_32X32)
      pc->txfm_mode += vp9_read_bit(r);

    if (pc->txfm_mode == TX_MODE_SELECT) {
      pc->prob_tx[0] = vp9_read_prob(r);
      pc->prob_tx[1] = vp9_read_prob(r);
      pc->prob_tx[2] = vp9_read_prob(r);
    }
  }
}

static int get_unsigned_bits(unsigned int num_values) {
  int cat = 0;
  if (num_values <= 1)
    return 0;
  num_values--;
  while (num_values > 0) {
    cat++;
    num_values >>= 1;
  }
  return cat;
}

static int inv_recenter_nonneg(int v, int m) {
  if (v > 2 * m)
    return v;

  return v % 2 ? m - (v + 1) / 2 : m + v / 2;
}

static int decode_uniform(vp9_reader *r, int n) {
  int v;
  const int l = get_unsigned_bits(n);
  const int m = (1 << l) - n;
  if (!l)
    return 0;

  v = vp9_read_literal(r, l - 1);
  return v < m ?  v : (v << 1) - m + vp9_read_bit(r);
}

static int decode_term_subexp(vp9_reader *r, int k, int num_syms) {
  int i = 0, mk = 0, word;
  while (1) {
    const int b = i ? k + i - 1 : k;
    const int a = 1 << b;
    if (num_syms <= mk + 3 * a) {
      word = decode_uniform(r, num_syms - mk) + mk;
      break;
    } else {
      if (vp9_read_bit(r)) {
        i++;
        mk += a;
      } else {
        word = vp9_read_literal(r, b) + mk;
        break;
      }
    }
  }
  return word;
}

static int decode_unsigned_max(vp9_reader *r, int max) {
  int data = 0, bit = 0, lmax = max;

  while (lmax) {
    data |= vp9_read_bit(r) << bit++;
    lmax >>= 1;
  }
  return data > max ? max : data;
}

static int merge_index(int v, int n, int modulus) {
  int max1 = (n - 1 - modulus / 2) / modulus + 1;
  if (v < max1) {
    v = v * modulus + modulus / 2;
  } else {
    int w;
    v -= max1;
    w = v;
    v += (v + modulus - modulus / 2) / modulus;
    while (v % modulus == modulus / 2 ||
           w != v - (v + modulus - modulus / 2) / modulus) v++;
  }
  return v;
}

static int inv_remap_prob(int v, int m) {
  const int n = 256;

  v = merge_index(v, n - 1, MODULUS_PARAM);
  if ((m << 1) <= n) {
    return inv_recenter_nonneg(v + 1, m);
  } else {
    return n - 1 - inv_recenter_nonneg(v + 1, n - 1 - m);
  }
}

static vp9_prob read_prob_diff_update(vp9_reader *r, int oldp) {
  int delp = decode_term_subexp(r, SUBEXP_PARAM, 255);
  return (vp9_prob)inv_remap_prob(delp, oldp);
}

void vp9_init_dequantizer(VP9_COMMON *pc) {
  int q;

  for (q = 0; q < QINDEX_RANGE; q++) {
    // DC value
    pc->y_dequant[q][0] = vp9_dc_quant(q, pc->y_dc_delta_q);
    pc->uv_dequant[q][0] = vp9_dc_quant(q, pc->uv_dc_delta_q);

    // AC values
    pc->y_dequant[q][1] = vp9_ac_quant(q, 0);
    pc->uv_dequant[q][1] = vp9_ac_quant(q, pc->uv_ac_delta_q);
  }
}

static void mb_init_dequantizer(VP9_COMMON *pc, MACROBLOCKD *xd) {
  int i;
  const int segment_id = xd->mode_info_context->mbmi.segment_id;
  xd->q_index = vp9_get_qindex(xd, segment_id, pc->base_qindex);

  xd->plane[0].dequant = pc->y_dequant[xd->q_index];
  for (i = 1; i < MAX_MB_PLANE; i++)
    xd->plane[i].dequant = pc->uv_dequant[xd->q_index];
}

static INLINE void dequant_add_y(MACROBLOCKD *xd, TX_TYPE tx_type, int idx,
                                 BLOCK_SIZE_TYPE bsize) {
  struct macroblockd_plane *const y = &xd->plane[0];
  uint8_t* const dst = raster_block_offset_uint8(xd, bsize, 0, idx,
                                                 xd->plane[0].dst.buf,
                                                 xd->plane[0].dst.stride);
  if (tx_type != DCT_DCT) {
    vp9_iht_add_c(tx_type, BLOCK_OFFSET(y->qcoeff, idx, 16),
                  dst, xd->plane[0].dst.stride, y->eobs[idx]);
  } else {
    xd->itxm_add(BLOCK_OFFSET(y->qcoeff, idx, 16),
                 dst, xd->plane[0].dst.stride, y->eobs[idx]);
  }
}

static void decode_block(int plane, int block, BLOCK_SIZE_TYPE bsize,
                         int ss_txfrm_size, void *arg) {
  MACROBLOCKD* const xd = arg;
  int16_t* const qcoeff = BLOCK_OFFSET(xd->plane[plane].qcoeff, block, 16);
  const int stride = xd->plane[plane].dst.stride;
  const int raster_block = txfrm_block_to_raster_block(xd, bsize, plane,
                                                       block, ss_txfrm_size);
  uint8_t* const dst = raster_block_offset_uint8(xd, bsize, plane,
                                                 raster_block,
                                                 xd->plane[plane].dst.buf,
                                                 stride);

  TX_TYPE tx_type;

  switch (ss_txfrm_size / 2) {
    case TX_4X4:
      tx_type = plane == 0 ? get_tx_type_4x4(xd, raster_block) : DCT_DCT;
      if (tx_type == DCT_DCT)
        xd->itxm_add(qcoeff, dst, stride, xd->plane[plane].eobs[block]);
      else
        vp9_iht_add_c(tx_type, qcoeff, dst, stride,
                      xd->plane[plane].eobs[block]);
      break;
    case TX_8X8:
      tx_type = plane == 0 ? get_tx_type_8x8(xd, raster_block) : DCT_DCT;
      vp9_iht_add_8x8_c(tx_type, qcoeff, dst, stride,
                        xd->plane[plane].eobs[block]);
      break;
    case TX_16X16:
      tx_type = plane == 0 ? get_tx_type_16x16(xd, raster_block) : DCT_DCT;
      vp9_iht_add_16x16_c(tx_type, qcoeff, dst, stride,
                          xd->plane[plane].eobs[block]);
      break;
    case TX_32X32:
      vp9_idct_add_32x32(qcoeff, dst, stride, xd->plane[plane].eobs[block]);
      break;
  }
}

static void decode_atom_intra(VP9D_COMP *pbi, MACROBLOCKD *xd,
                              vp9_reader *r,
                              BLOCK_SIZE_TYPE bsize) {
  int i = 0;
  int bwl = b_width_log2(bsize), bhl = b_height_log2(bsize);
  int bc = 1 << (bwl + bhl);
  int tx_type;

  for (i = 0; i < bc; i++) {
    int b_mode = xd->mode_info_context->bmi[i].as_mode.first;

    uint8_t* dst = raster_block_offset_uint8(xd, bsize, 0, i,
                                             xd->plane[0].dst.buf,
                                             xd->plane[0].dst.stride);

    vp9_intra4x4_predict(xd, i, bsize, b_mode, dst, xd->plane[0].dst.stride);
    // TODO(jingning): refactor to use foreach_transformed_block_in_plane_
    tx_type = get_tx_type_4x4(xd, i);
    dequant_add_y(xd, tx_type, i, bsize);
  }

  foreach_transformed_block_uv(xd, bsize, decode_block, xd);
}

static void decode_atom(VP9D_COMP *pbi, MACROBLOCKD *xd,
                        int mi_row, int mi_col,
                        vp9_reader *r, BLOCK_SIZE_TYPE bsize) {
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;

  if (pbi->common.frame_type != KEY_FRAME)
    vp9_setup_interp_filters(xd, mbmi->interp_filter, &pbi->common);

  // prediction
  if (mbmi->ref_frame == INTRA_FRAME)
    vp9_build_intra_predictors_sbuv_s(xd, bsize);
  else
    vp9_build_inter_predictors_sb(xd, mi_row, mi_col, bsize);

  if (mbmi->mb_skip_coeff) {
    vp9_reset_sb_tokens_context(xd, bsize);
  } else {
    // re-initialize macroblock dequantizer before detokenization
    if (xd->segmentation_enabled)
      mb_init_dequantizer(&pbi->common, xd);

    if (!vp9_reader_has_error(r)) {
      vp9_decode_tokens(pbi, xd, r, bsize);
    }
  }

  if (mbmi->ref_frame == INTRA_FRAME)
    decode_atom_intra(pbi, xd, r, bsize);
  else
    foreach_transformed_block(xd, bsize, decode_block, xd);
}

static void decode_sb(VP9D_COMP *pbi, MACROBLOCKD *xd, int mi_row, int mi_col,
                      vp9_reader *r, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mi_width_log2(bsize), bhl = mi_height_log2(bsize);
  const int bw = 1 << bwl, bh = 1 << bhl;
  int n, eobtotal;
  VP9_COMMON *const pc = &pbi->common;
  MODE_INFO *const mi = xd->mode_info_context;
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const int mis = pc->mode_info_stride;

  assert(mbmi->sb_type == bsize);

  if (pbi->common.frame_type != KEY_FRAME)
    vp9_setup_interp_filters(xd, mbmi->interp_filter, pc);

  // generate prediction
  if (mbmi->ref_frame == INTRA_FRAME) {
    vp9_build_intra_predictors_sby_s(xd, bsize);
    vp9_build_intra_predictors_sbuv_s(xd, bsize);
  } else {
    vp9_build_inter_predictors_sb(xd, mi_row, mi_col, bsize);
  }

  if (mbmi->mb_skip_coeff) {
    vp9_reset_sb_tokens_context(xd, bsize);
  } else {
    // re-initialize macroblock dequantizer before detokenization
    if (xd->segmentation_enabled)
      mb_init_dequantizer(pc, xd);

    // dequantization and idct
    eobtotal = vp9_decode_tokens(pbi, xd, r, bsize);
    if (eobtotal == 0) {  // skip loopfilter
      for (n = 0; n < bw * bh; n++) {
        const int x_idx = n & (bw - 1), y_idx = n >> bwl;

        if (mi_col + x_idx < pc->mi_cols && mi_row + y_idx < pc->mi_rows)
          mi[y_idx * mis + x_idx].mbmi.mb_skip_coeff = 1;
      }
    } else {
      foreach_transformed_block(xd, bsize, decode_block, xd);
    }
  }
}

static int get_delta_q(vp9_reader *r, int *dq) {
  const int old_value = *dq;

  if (vp9_read_bit(r)) {  // Update bit
    const int value = vp9_read_literal(r, 4);
    *dq = vp9_read_and_apply_sign(r, value);
  }

  // Trigger a quantizer update if the delta-q value has changed
  return old_value != *dq;
}

static void set_offsets(VP9D_COMP *pbi, BLOCK_SIZE_TYPE bsize,
                        int mi_row, int mi_col) {
  const int bh = 1 << mi_height_log2(bsize);
  const int bw = 1 << mi_width_log2(bsize);
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  const int mi_idx = mi_row * cm->mode_info_stride + mi_col;
  int i;

  xd->mode_info_context = cm->mi + mi_idx;
  xd->mode_info_context->mbmi.sb_type = bsize;
  xd->prev_mode_info_context = cm->prev_mi + mi_idx;

  for (i = 0; i < MAX_MB_PLANE; i++) {
    xd->plane[i].above_context = cm->above_context[i] +
        (mi_col * 2 >> xd->plane[i].subsampling_x);
    xd->plane[i].left_context = cm->left_context[i] +
        (((mi_row * 2) & 15) >> xd->plane[i].subsampling_y);
  }
  xd->above_seg_context = cm->above_seg_context + mi_col;
  xd->left_seg_context  = cm->left_seg_context + (mi_row & MI_MASK);

  // Distance of Mb to the various image edges. These are specified to 8th pel
  // as they are always compared to values that are in 1/8th pel units
  set_mi_row_col(cm, xd, mi_row, bh, mi_col, bw);

  setup_dst_planes(xd, &cm->yv12_fb[cm->new_fb_idx], mi_row, mi_col);
}

static void set_refs(VP9D_COMP *pbi, int mi_row, int mi_col) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;

  if (mbmi->ref_frame > INTRA_FRAME) {
    // Select the appropriate reference frame for this MB
    const int fb_idx = cm->active_ref_idx[mbmi->ref_frame - 1];
    const YV12_BUFFER_CONFIG *cfg = &cm->yv12_fb[fb_idx];
    xd->scale_factor[0]    = cm->active_ref_scale[mbmi->ref_frame - 1];
    xd->scale_factor_uv[0] = cm->active_ref_scale[mbmi->ref_frame - 1];
    setup_pre_planes(xd, cfg, NULL, mi_row, mi_col,
                     xd->scale_factor, xd->scale_factor_uv);
    xd->corrupted |= cfg->corrupted;

    if (mbmi->second_ref_frame > INTRA_FRAME) {
      // Select the appropriate reference frame for this MB
      const int second_fb_idx = cm->active_ref_idx[mbmi->second_ref_frame - 1];
      const YV12_BUFFER_CONFIG *second_cfg = &cm->yv12_fb[second_fb_idx];
      xd->scale_factor[1]    = cm->active_ref_scale[mbmi->second_ref_frame - 1];
      xd->scale_factor_uv[1] = cm->active_ref_scale[mbmi->second_ref_frame - 1];
      setup_pre_planes(xd, NULL, second_cfg, mi_row, mi_col,
                       xd->scale_factor, xd->scale_factor_uv);
      xd->corrupted |= second_cfg->corrupted;
    }
  }
}

static void decode_modes_b(VP9D_COMP *pbi, int mi_row, int mi_col,
                           vp9_reader *r, BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &pbi->mb;

  set_offsets(pbi, bsize, mi_row, mi_col);
  vp9_decode_mb_mode_mv(pbi, xd, mi_row, mi_col, r);
  set_refs(pbi, mi_row, mi_col);

#if CONFIG_AB4X4
  if (bsize < BLOCK_SIZE_SB8X8)
#else
  if (bsize == BLOCK_SIZE_SB8X8 &&
      (xd->mode_info_context->mbmi.mode == SPLITMV ||
       xd->mode_info_context->mbmi.mode == I4X4_PRED))
#endif
    decode_atom(pbi, xd, mi_row, mi_col, r, BLOCK_SIZE_SB8X8);
  else
    decode_sb(pbi, xd, mi_row, mi_col, r, bsize);

  xd->corrupted |= vp9_reader_has_error(r);
}

static void decode_modes_sb(VP9D_COMP *pbi, int mi_row, int mi_col,
                            vp9_reader* r, BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON *const pc = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  int bsl = mi_width_log2(bsize), bs = (1 << bsl) / 2;
  int n;
  PARTITION_TYPE partition = PARTITION_NONE;
  BLOCK_SIZE_TYPE subsize;

  if (mi_row >= pc->mi_rows || mi_col >= pc->mi_cols)
    return;

#if CONFIG_AB4X4
  if (bsize < BLOCK_SIZE_SB8X8)
    if (xd->ab_index != 0)
      return;
#endif

#if CONFIG_AB4X4
  if (bsize >= BLOCK_SIZE_SB8X8) {
#else
  if (bsize > BLOCK_SIZE_SB8X8) {
#endif
    int pl;
    // read the partition information
    xd->left_seg_context = pc->left_seg_context + (mi_row & MI_MASK);
    xd->above_seg_context = pc->above_seg_context + mi_col;
    pl = partition_plane_context(xd, bsize);
    partition = treed_read(r, vp9_partition_tree,
                           pc->fc.partition_prob[pl]);
    pc->fc.partition_counts[pl][partition]++;
  }

  subsize = get_subsize(bsize, partition);

  switch (partition) {
    case PARTITION_NONE:
      decode_modes_b(pbi, mi_row, mi_col, r, subsize);
      break;
    case PARTITION_HORZ:
      decode_modes_b(pbi, mi_row, mi_col, r, subsize);
      if (mi_row + bs < pc->mi_rows)
        decode_modes_b(pbi, mi_row + bs, mi_col, r, subsize);
      break;
    case PARTITION_VERT:
      decode_modes_b(pbi, mi_row, mi_col, r, subsize);
      if (mi_col + bs < pc->mi_cols)
        decode_modes_b(pbi, mi_row, mi_col + bs, r, subsize);
      break;
    case PARTITION_SPLIT:
      for (n = 0; n < 4; n++) {
        int j = n >> 1, i = n & 0x01;
        *(get_sb_index(xd, subsize)) = n;
        decode_modes_sb(pbi, mi_row + j * bs, mi_col + i * bs, r, subsize);
      }
      break;
    default:
      assert(0);
  }
  // update partition context
#if CONFIG_AB4X4
  if (bsize >= BLOCK_SIZE_SB8X8 &&
      (bsize == BLOCK_SIZE_SB8X8 || partition != PARTITION_SPLIT)) {
#else
  if (bsize > BLOCK_SIZE_SB8X8 &&
      (bsize == BLOCK_SIZE_MB16X16 || partition != PARTITION_SPLIT)) {
#endif
    set_partition_seg_context(pc, xd, mi_row, mi_col);
    update_partition_context(xd, subsize, bsize);
  }
}

static void setup_token_decoder(VP9D_COMP *pbi,
                                const uint8_t *data,
                                vp9_reader *r) {
  VP9_COMMON *pc = &pbi->common;
  const uint8_t *data_end = pbi->source + pbi->source_sz;
  const size_t partition_size = data_end - data;

  // Validate the calculated partition length. If the buffer
  // described by the partition can't be fully read, then restrict
  // it to the portion that can be (for EC mode) or throw an error.
  if (!read_is_valid(data, partition_size, data_end))
    vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt partition "
                       "%d length", 1);

  if (vp9_reader_init(r, data, partition_size))
    vpx_internal_error(&pc->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder %d", 1);
}

static void init_frame(VP9D_COMP *pbi) {
  VP9_COMMON *const pc = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;

  if (pc->frame_type == KEY_FRAME) {
    vp9_setup_past_independence(pc, xd);
    // All buffers are implicitly updated on key frames.
    pbi->refresh_frame_flags = (1 << NUM_REF_FRAMES) - 1;
  } else if (pc->error_resilient_mode) {
    vp9_setup_past_independence(pc, xd);
  }

  xd->mode_info_context = pc->mi;
  xd->prev_mode_info_context = pc->prev_mi;
  xd->frame_type = pc->frame_type;
  xd->mode_info_context->mbmi.mode = DC_PRED;
  xd->mode_info_stride = pc->mode_info_stride;
}

static void read_coef_probs_common(vp9_coeff_probs *coef_probs,
                                   TX_SIZE tx_size,
                                   vp9_reader *r) {
#if CONFIG_MODELCOEFPROB && MODEL_BASED_UPDATE
  const int entropy_nodes_update = UNCONSTRAINED_UPDATE_NODES;
#else
  const int entropy_nodes_update = ENTROPY_NODES;
#endif

  int i, j, k, l, m;

  if (vp9_read_bit(r)) {
    for (i = 0; i < BLOCK_TYPES; i++) {
      for (j = 0; j < REF_TYPES; j++) {
        for (k = 0; k < COEF_BANDS; k++) {
          for (l = 0; l < PREV_COEF_CONTEXTS; l++) {
            const int mstart = 0;
            if (l >= 3 && k == 0)
              continue;

            for (m = mstart; m < entropy_nodes_update; m++) {
              vp9_prob *const p = coef_probs[i][j][k][l] + m;

              if (vp9_read(r, vp9_coef_update_prob[m])) {
                *p = read_prob_diff_update(r, *p);
#if CONFIG_MODELCOEFPROB && MODEL_BASED_UPDATE
                if (m == UNCONSTRAINED_NODES - 1)
                  vp9_get_model_distribution(*p, coef_probs[i][j][k][l], i, j);
#endif
              }
            }
          }
        }
      }
    }
  }
}

static void read_coef_probs(VP9D_COMP *pbi, vp9_reader *r) {
  const TXFM_MODE mode = pbi->common.txfm_mode;
  FRAME_CONTEXT *const fc = &pbi->common.fc;

  read_coef_probs_common(fc->coef_probs_4x4, TX_4X4, r);

  if (mode > ONLY_4X4)
    read_coef_probs_common(fc->coef_probs_8x8, TX_8X8, r);

  if (mode > ALLOW_8X8)
    read_coef_probs_common(fc->coef_probs_16x16, TX_16X16, r);

  if (mode > ALLOW_16X16)
    read_coef_probs_common(fc->coef_probs_32x32, TX_32X32, r);
}

static void setup_segmentation(VP9_COMMON *pc, MACROBLOCKD *xd, vp9_reader *r) {
  int i, j;

  xd->update_mb_segmentation_map = 0;
  xd->update_mb_segmentation_data = 0;
#if CONFIG_IMPLICIT_SEGMENTATION
  xd->allow_implicit_segment_update = 0;
#endif

  xd->segmentation_enabled = vp9_read_bit(r);
  if (!xd->segmentation_enabled)
    return;

  // Segmentation map update
  xd->update_mb_segmentation_map = vp9_read_bit(r);
#if CONFIG_IMPLICIT_SEGMENTATION
    xd->allow_implicit_segment_update = vp9_read_bit(r);
#endif
  if (xd->update_mb_segmentation_map) {
    for (i = 0; i < MB_SEG_TREE_PROBS; i++)
      xd->mb_segment_tree_probs[i] = vp9_read_bit(r) ? vp9_read_prob(r)
                                                     : MAX_PROB;

    pc->temporal_update = vp9_read_bit(r);
    if (pc->temporal_update) {
      for (i = 0; i < PREDICTION_PROBS; i++)
        pc->segment_pred_probs[i] = vp9_read_bit(r) ? vp9_read_prob(r)
                                                    : MAX_PROB;
    } else {
      for (i = 0; i < PREDICTION_PROBS; i++)
        pc->segment_pred_probs[i] = MAX_PROB;
    }
  }

  // Segmentation data update
  xd->update_mb_segmentation_data = vp9_read_bit(r);
  if (xd->update_mb_segmentation_data) {
    xd->mb_segment_abs_delta = vp9_read_bit(r);

    vp9_clearall_segfeatures(xd);

    for (i = 0; i < MAX_MB_SEGMENTS; i++) {
      for (j = 0; j < SEG_LVL_MAX; j++) {
        int data = 0;
        const int feature_enabled = vp9_read_bit(r);
        if (feature_enabled) {
          vp9_enable_segfeature(xd, i, j);
          data = decode_unsigned_max(r, vp9_seg_feature_data_max(j));
          if (vp9_is_segfeature_signed(j))
            data = vp9_read_and_apply_sign(r, data);
        }
        vp9_set_segdata(xd, i, j, data);
      }
    }
  }
}

static void setup_pred_probs(VP9_COMMON *pc, vp9_reader *r) {
  // Read common prediction model status flag probability updates for the
  // reference frame
  if (pc->frame_type == KEY_FRAME) {
    // Set the prediction probabilities to defaults
    pc->ref_pred_probs[0] = DEFAULT_PRED_PROB_0;
    pc->ref_pred_probs[1] = DEFAULT_PRED_PROB_1;
    pc->ref_pred_probs[2] = DEFAULT_PRED_PROB_2;
  } else {
    int i;
    for (i = 0; i < PREDICTION_PROBS; ++i)
      if (vp9_read_bit(r))
        pc->ref_pred_probs[i] = vp9_read_prob(r);
  }
}

static void setup_loopfilter(VP9_COMMON *pc, MACROBLOCKD *xd, vp9_reader *r) {
  pc->filter_level = vp9_read_literal(r, 6);
  pc->sharpness_level = vp9_read_literal(r, 3);

#if CONFIG_LOOP_DERING
  if (vp9_read_bit(r))
    pc->dering_enabled = 1 + vp9_read_literal(r, 4);
  else
    pc->dering_enabled = 0;
#endif

  // Read in loop filter deltas applied at the MB level based on mode or ref
  // frame.
  xd->mode_ref_lf_delta_update = 0;

  xd->mode_ref_lf_delta_enabled = vp9_read_bit(r);
  if (xd->mode_ref_lf_delta_enabled) {
    xd->mode_ref_lf_delta_update = vp9_read_bit(r);
    if (xd->mode_ref_lf_delta_update) {
      int i;

      for (i = 0; i < MAX_REF_LF_DELTAS; i++) {
        if (vp9_read_bit(r)) {
          const int value = vp9_read_literal(r, 6);
          xd->ref_lf_deltas[i] = vp9_read_and_apply_sign(r, value);
        }
      }

      for (i = 0; i < MAX_MODE_LF_DELTAS; i++) {
        if (vp9_read_bit(r)) {
          const int value = vp9_read_literal(r, 6);
          xd->mode_lf_deltas[i] = vp9_read_and_apply_sign(r, value);
        }
      }
    }
  }
}

static void setup_quantization(VP9D_COMP *pbi, vp9_reader *r) {
  // Read the default quantizers
  VP9_COMMON *const pc = &pbi->common;

  pc->base_qindex = vp9_read_literal(r, QINDEX_BITS);
  if (get_delta_q(r, &pc->y_dc_delta_q) |
      get_delta_q(r, &pc->uv_dc_delta_q) |
      get_delta_q(r, &pc->uv_ac_delta_q))
    vp9_init_dequantizer(pc);

  mb_init_dequantizer(pc, &pbi->mb);  // MB level dequantizer setup
}

static INTERPOLATIONFILTERTYPE read_mcomp_filter_type(vp9_reader *r) {
  return vp9_read_bit(r) ? SWITCHABLE
                         : vp9_read_literal(r, 2);
}

static const uint8_t *read_frame_size(VP9_COMMON *const pc, const uint8_t *data,
                                      const uint8_t *data_end,
                                      int *width, int *height) {
  if (data + 4 < data_end) {
    const int w = read_le16(data);
    const int h = read_le16(data + 2);
    if (w <= 0)
      vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                         "Invalid frame width");

    if (h <= 0)
      vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                         "Invalid frame height");
    *width = w;
    *height = h;
    data += 4;
  } else {
    vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                       "Failed to read frame size");
  }
  return data;
}

static const uint8_t *setup_frame_size(VP9D_COMP *pbi, int scaling_active,
                                       const uint8_t *data,
                                       const uint8_t *data_end) {
  // If error concealment is enabled we should only parse the new size
  // if we have enough data. Otherwise we will end up with the wrong size.
  VP9_COMMON *const pc = &pbi->common;
  int display_width = pc->display_width;
  int display_height = pc->display_height;
  int width = pc->width;
  int height = pc->height;

  if (scaling_active)
    data = read_frame_size(pc, data, data_end, &display_width, &display_height);

  data = read_frame_size(pc, data, data_end, &width, &height);

  if (pc->width != width || pc->height != height) {
    if (!pbi->initial_width || !pbi->initial_height) {
      if (vp9_alloc_frame_buffers(pc, width, height))
        vpx_internal_error(&pc->error, VPX_CODEC_MEM_ERROR,
                           "Failed to allocate frame buffers");
        pbi->initial_width = width;
        pbi->initial_height = height;
    } else {
      if (width > pbi->initial_width)
        vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                           "Frame width too large");

      if (height > pbi->initial_height)
        vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                           "Frame height too large");
    }

    pc->width = width;
    pc->height = height;
    pc->display_width = scaling_active ? display_width : width;
    pc->display_height = scaling_active ? display_height : height;

    vp9_update_frame_size(pc);
  }

  return data;
}

static void update_frame_context(FRAME_CONTEXT *fc) {
  vp9_copy(fc->pre_coef_probs_4x4, fc->coef_probs_4x4);
  vp9_copy(fc->pre_coef_probs_8x8, fc->coef_probs_8x8);
  vp9_copy(fc->pre_coef_probs_16x16, fc->coef_probs_16x16);
  vp9_copy(fc->pre_coef_probs_32x32, fc->coef_probs_32x32);
  vp9_copy(fc->pre_ymode_prob, fc->ymode_prob);
  vp9_copy(fc->pre_sb_ymode_prob, fc->sb_ymode_prob);
  vp9_copy(fc->pre_uv_mode_prob, fc->uv_mode_prob);
  vp9_copy(fc->pre_bmode_prob, fc->bmode_prob);
  vp9_copy(fc->pre_sub_mv_ref_prob, fc->sub_mv_ref_prob);
  vp9_copy(fc->pre_partition_prob, fc->partition_prob);
  fc->pre_nmvc = fc->nmvc;

  vp9_zero(fc->coef_counts_4x4);
  vp9_zero(fc->coef_counts_8x8);
  vp9_zero(fc->coef_counts_16x16);
  vp9_zero(fc->coef_counts_32x32);
  vp9_zero(fc->eob_branch_counts);
  vp9_zero(fc->ymode_counts);
  vp9_zero(fc->sb_ymode_counts);
  vp9_zero(fc->uv_mode_counts);
  vp9_zero(fc->bmode_counts);
  vp9_zero(fc->sub_mv_ref_counts);
  vp9_zero(fc->NMVcount);
  vp9_zero(fc->mv_ref_ct);
  vp9_zero(fc->partition_counts);
}

static void decode_tile(VP9D_COMP *pbi, vp9_reader *r) {
  VP9_COMMON *const pc = &pbi->common;
  int mi_row, mi_col;

  for (mi_row = pc->cur_tile_mi_row_start;
       mi_row < pc->cur_tile_mi_row_end; mi_row += 64 / MI_SIZE) {
    // For a SB there are 2 left contexts, each pertaining to a MB row within
    vpx_memset(&pc->left_context, 0, sizeof(pc->left_context));
    vpx_memset(pc->left_seg_context, 0, sizeof(pc->left_seg_context));
    for (mi_col = pc->cur_tile_mi_col_start;
         mi_col < pc->cur_tile_mi_col_end; mi_col += 64 / MI_SIZE)
      decode_modes_sb(pbi, mi_row, mi_col, r, BLOCK_SIZE_SB64X64);
  }
}

static void decode_tiles(VP9D_COMP *pbi,
                         const uint8_t *data, int first_partition_size,
                         vp9_reader *header_bc, vp9_reader *residual_bc) {
  VP9_COMMON *const pc = &pbi->common;

  const uint8_t *data_ptr = data + first_partition_size;
  int tile_row, tile_col, delta_log2_tiles;

  vp9_get_tile_n_bits(pc, &pc->log2_tile_columns, &delta_log2_tiles);
  while (delta_log2_tiles--) {
    if (vp9_read_bit(header_bc)) {
      pc->log2_tile_columns++;
    } else {
      break;
    }
  }
  pc->log2_tile_rows = vp9_read_bit(header_bc);
  if (pc->log2_tile_rows)
    pc->log2_tile_rows += vp9_read_bit(header_bc);
  pc->tile_columns = 1 << pc->log2_tile_columns;
  pc->tile_rows    = 1 << pc->log2_tile_rows;

  // Note: this memset assumes above_context[0], [1] and [2]
  // are allocated as part of the same buffer.
  vpx_memset(pc->above_context[0], 0, sizeof(ENTROPY_CONTEXT) * 2 *
                                      MAX_MB_PLANE * mi_cols_aligned_to_sb(pc));

  vpx_memset(pc->above_seg_context, 0, sizeof(PARTITION_CONTEXT) *
                                       mi_cols_aligned_to_sb(pc));

  if (pbi->oxcf.inv_tile_order) {
    const int n_cols = pc->tile_columns;
    const uint8_t *data_ptr2[4][1 << 6];
    vp9_reader bc_bak = {0};

    // pre-initialize the offsets, we're going to read in inverse order
    data_ptr2[0][0] = data_ptr;
    for (tile_row = 0; tile_row < pc->tile_rows; tile_row++) {
      if (tile_row) {
        const int size = read_le32(data_ptr2[tile_row - 1][n_cols - 1]);
        data_ptr2[tile_row - 1][n_cols - 1] += 4;
        data_ptr2[tile_row][0] = data_ptr2[tile_row - 1][n_cols - 1] + size;
      }

      for (tile_col = 1; tile_col < n_cols; tile_col++) {
        const int size = read_le32(data_ptr2[tile_row][tile_col - 1]);
        data_ptr2[tile_row][tile_col - 1] += 4;
        data_ptr2[tile_row][tile_col] =
            data_ptr2[tile_row][tile_col - 1] + size;
      }
    }

    for (tile_row = 0; tile_row < pc->tile_rows; tile_row++) {
      vp9_get_tile_row_offsets(pc, tile_row);
      for (tile_col = n_cols - 1; tile_col >= 0; tile_col--) {
        vp9_get_tile_col_offsets(pc, tile_col);
        setup_token_decoder(pbi, data_ptr2[tile_row][tile_col], residual_bc);
        decode_tile(pbi, residual_bc);
        if (tile_row == pc->tile_rows - 1 && tile_col == n_cols - 1)
          bc_bak = *residual_bc;
      }
    }
    *residual_bc = bc_bak;
  } else {
    int has_more;

    for (tile_row = 0; tile_row < pc->tile_rows; tile_row++) {
      vp9_get_tile_row_offsets(pc, tile_row);
      for (tile_col = 0; tile_col < pc->tile_columns; tile_col++) {
        vp9_get_tile_col_offsets(pc, tile_col);

        has_more = tile_col < pc->tile_columns - 1 ||
                   tile_row < pc->tile_rows - 1;

        setup_token_decoder(pbi, data_ptr + (has_more ? 4 : 0), residual_bc);
        decode_tile(pbi, residual_bc);

        if (has_more) {
          const int size = read_le32(data_ptr);
          data_ptr += 4 + size;
        }
      }
    }
  }
}

int vp9_decode_frame(VP9D_COMP *pbi, const uint8_t **p_data_end) {
  vp9_reader header_bc, residual_bc;
  VP9_COMMON *const pc = &pbi->common;
  MACROBLOCKD *const xd  = &pbi->mb;
  const uint8_t *data = pbi->source;
  const uint8_t *data_end = data + pbi->source_sz;
  size_t first_partition_size = 0;
  YV12_BUFFER_CONFIG *new_fb = &pc->yv12_fb[pc->new_fb_idx];
  int i;

  xd->corrupted = 0;  // start with no corruption of current frame
  new_fb->corrupted = 0;

  if (data_end - data < 3) {
    vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME, "Truncated packet");
  } else {
    int scaling_active;
    pc->last_frame_type = pc->frame_type;
    pc->frame_type = (FRAME_TYPE)(data[0] & 1);
    pc->version = (data[0] >> 1) & 7;
    pc->show_frame = (data[0] >> 4) & 1;
    scaling_active = (data[0] >> 5) & 1;
    pc->subsampling_x = (data[0] >> 6) & 1;
    pc->subsampling_y = (data[0] >> 7) & 1;
    first_partition_size = read_le16(data + 1);

    if (!read_is_valid(data, first_partition_size, data_end))
      vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                         "Truncated packet or corrupt partition 0 length");

    data += 3;

    vp9_setup_version(pc);

    if (pc->frame_type == KEY_FRAME) {
      // When error concealment is enabled we should only check the sync
      // code if we have enough bits available
      if (data + 3 < data_end) {
        if (data[0] != 0x49 || data[1] != 0x83 || data[2] != 0x42)
          vpx_internal_error(&pc->error, VPX_CODEC_UNSUP_BITSTREAM,
                             "Invalid frame sync code");
      }
      data += 3;
    }

    data = setup_frame_size(pbi, scaling_active, data, data_end);
  }

  if ((!pbi->decoded_key_frame && pc->frame_type != KEY_FRAME) ||
      pc->width == 0 || pc->height == 0) {
    return -1;
  }

  init_frame(pbi);

  // Reset the frame pointers to the current frame size
  vp9_realloc_frame_buffer(new_fb, pc->width, pc->height,
                           pc->subsampling_x, pc->subsampling_y,
                           VP9BORDERINPIXELS);

  if (vp9_reader_init(&header_bc, data, first_partition_size))
    vpx_internal_error(&pc->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder 0");

  pc->clr_type = (YUV_TYPE)vp9_read_bit(&header_bc);
  pc->clamp_type = (CLAMP_TYPE)vp9_read_bit(&header_bc);
  pc->error_resilient_mode = vp9_read_bit(&header_bc);

  setup_loopfilter(pc, xd, &header_bc);

  setup_quantization(pbi, &header_bc);

  xd->lossless = pc->base_qindex == 0 &&
                 pc->y_dc_delta_q == 0 &&
                 pc->uv_dc_delta_q == 0 &&
                 pc->uv_ac_delta_q == 0;
  if (xd->lossless) {
    xd->inv_txm4x4_1      = vp9_short_iwalsh4x4_1;
    xd->inv_txm4x4        = vp9_short_iwalsh4x4;
    xd->itxm_add          = vp9_idct_add_lossless_c;
    xd->itxm_add_y_block  = vp9_idct_add_y_block_lossless_c;
    xd->itxm_add_uv_block = vp9_idct_add_uv_block_lossless_c;
  } else {
    xd->inv_txm4x4_1      = vp9_short_idct4x4_1;
    xd->inv_txm4x4        = vp9_short_idct4x4;
    xd->itxm_add          = vp9_idct_add;
    xd->itxm_add_y_block  = vp9_idct_add_y_block;
    xd->itxm_add_uv_block = vp9_idct_add_uv_block;
  }

  // Determine if the golden frame or ARF buffer should be updated and how.
  // For all non key frames the GF and ARF refresh flags and sign bias
  // flags must be set explicitly.
  if (pc->frame_type == KEY_FRAME) {
    for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i)
      pc->active_ref_idx[i] = pc->new_fb_idx;
  } else {
    // Should the GF or ARF be updated from the current frame
    pbi->refresh_frame_flags = vp9_read_literal(&header_bc, NUM_REF_FRAMES);

    // Select active reference frames and calculate scaling factors
    for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i) {
      const int ref = vp9_read_literal(&header_bc, NUM_REF_FRAMES_LG2);
      pc->active_ref_idx[i] = pc->ref_frame_map[ref];
      vp9_setup_scale_factors(pc, i);
    }

    // Read the sign bias for each reference frame buffer.
    for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i)
      pc->ref_frame_sign_bias[i + 1] = vp9_read_bit(&header_bc);

    xd->allow_high_precision_mv = vp9_read_bit(&header_bc);
    pc->mcomp_filter_type = read_mcomp_filter_type(&header_bc);

    // To enable choice of different interpolation filters
    vp9_setup_interp_filters(xd, pc->mcomp_filter_type, pc);
  }

  if (!pc->error_resilient_mode) {
    pc->refresh_frame_context = vp9_read_bit(&header_bc);
    pc->frame_parallel_decoding_mode = vp9_read_bit(&header_bc);
  } else {
    pc->refresh_frame_context = 0;
    pc->frame_parallel_decoding_mode = 1;
  }

  pc->frame_context_idx = vp9_read_literal(&header_bc, NUM_FRAME_CONTEXTS_LG2);
  pc->fc = pc->frame_contexts[pc->frame_context_idx];

  setup_segmentation(pc, xd, &header_bc);

  setup_pred_probs(pc, &header_bc);

  setup_txfm_mode(pc, xd->lossless, &header_bc);

  // Read inter mode probability context updates
  if (pc->frame_type != KEY_FRAME) {
    int i, j;
    for (i = 0; i < INTER_MODE_CONTEXTS; ++i)
      for (j = 0; j < 4; ++j)
        if (vp9_read(&header_bc, 252))
          pc->fc.vp9_mode_contexts[i][j] = vp9_read_prob(&header_bc);
  }
#if CONFIG_MODELCOEFPROB
  if (pc->frame_type == KEY_FRAME)
    vp9_default_coef_probs(pc);
#endif

  update_frame_context(&pc->fc);

  read_coef_probs(pbi, &header_bc);

  // Initialize xd pointers. Any reference should do for xd->pre, so use 0.
  setup_pre_planes(xd, &pc->yv12_fb[pc->active_ref_idx[0]], NULL,
                   0, 0, NULL, NULL);
  setup_dst_planes(xd, new_fb, 0, 0);

  // Create the segmentation map structure and set to 0
  if (!pc->last_frame_seg_map)
    CHECK_MEM_ERROR(pc->last_frame_seg_map,
                    vpx_calloc((pc->mi_rows * pc->mi_cols), 1));

  vp9_setup_block_dptrs(xd, pc->subsampling_x, pc->subsampling_y);

  // clear out the coeff buffer
  for (i = 0; i < MAX_MB_PLANE; ++i)
    vp9_zero(xd->plane[i].qcoeff);

  vp9_decode_mode_mvs_init(pbi, &header_bc);

  decode_tiles(pbi, data, first_partition_size, &header_bc, &residual_bc);

  pc->last_width = pc->width;
  pc->last_height = pc->height;

  new_fb->corrupted = vp9_reader_has_error(&header_bc) | xd->corrupted;

  if (!pbi->decoded_key_frame) {
    if (pc->frame_type == KEY_FRAME && !new_fb->corrupted)
      pbi->decoded_key_frame = 1;
    else
      vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                         "A stream must start with a complete key frame");
  }

  // Adaptation
  if (!pc->error_resilient_mode && !pc->frame_parallel_decoding_mode) {
    vp9_adapt_coef_probs(pc);

    if (pc->frame_type != KEY_FRAME) {
      vp9_adapt_mode_probs(pc);
      vp9_adapt_mode_context(pc);
      vp9_adapt_nmv_probs(pc, xd->allow_high_precision_mv);
    }
  }

#if CONFIG_IMPLICIT_SEGMENTATION
  // If signalled at the frame level apply implicit updates to the segment map.
  if (!pc->error_resilient_mode && xd->allow_implicit_segment_update) {
    vp9_implicit_segment_map_update(pc);
  }
#endif

  if (pc->refresh_frame_context)
    pc->frame_contexts[pc->frame_context_idx] = pc->fc;

  *p_data_end = vp9_reader_find_end(&residual_bc);
  return 0;
}
