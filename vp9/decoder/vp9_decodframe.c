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

#include "./vp9_rtcd.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_scale/vpx_scale.h"

#include "vp9/common/vp9_extend.h"
#include "vp9/common/vp9_modecont.h"
#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_alloccommon.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_quant_common.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_tile_common.h"

#include "vp9/decoder/vp9_dboolhuff.h"
#include "vp9/decoder/vp9_decodframe.h"
#include "vp9/decoder/vp9_detokenize.h"
#include "vp9/decoder/vp9_decodemv.h"
#include "vp9/decoder/vp9_onyxd_int.h"
#include "vp9/decoder/vp9_read_bit_buffer.h"


// #define DEC_DEBUG
#ifdef DEC_DEBUG
int dec_debug = 0;
#endif

static int read_be32(const uint8_t *p) {
  return (p[0] << 24) | (p[1] << 16) | (p[2] << 8) | p[3];
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
      int i, j;
      for (i = 0; i < TX_SIZE_CONTEXTS; ++i) {
        for (j = 0; j < TX_SIZE_MAX_SB - 3; ++j) {
          if (vp9_read(r, VP9_MODE_UPDATE_PROB))
            pc->fc.tx_probs_8x8p[i][j] =
                vp9_read_prob_diff_update(r, pc->fc.tx_probs_8x8p[i][j]);
        }
      }
      for (i = 0; i < TX_SIZE_CONTEXTS; ++i) {
        for (j = 0; j < TX_SIZE_MAX_SB - 2; ++j) {
          if (vp9_read(r, VP9_MODE_UPDATE_PROB))
            pc->fc.tx_probs_16x16p[i][j] =
                vp9_read_prob_diff_update(r, pc->fc.tx_probs_16x16p[i][j]);
        }
      }
      for (i = 0; i < TX_SIZE_CONTEXTS; ++i) {
        for (j = 0; j < TX_SIZE_MAX_SB - 1; ++j) {
          if (vp9_read(r, VP9_MODE_UPDATE_PROB))
            pc->fc.tx_probs_32x32p[i][j] =
                vp9_read_prob_diff_update(r, pc->fc.tx_probs_32x32p[i][j]);
        }
      }
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

static int decode_unsigned_max(struct vp9_read_bit_buffer *rb, int max) {
  const int data = vp9_rb_read_literal(rb, get_unsigned_bits(max));
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
  const int n = 255;

  v = merge_index(v, n - 1, MODULUS_PARAM);
  m--;
  if ((m << 1) <= n) {
    return 1 + inv_recenter_nonneg(v + 1, m);
  } else {
    return n - inv_recenter_nonneg(v + 1, n - 1 - m);
  }
}

vp9_prob vp9_read_prob_diff_update(vp9_reader *r, int oldp) {
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

static void decode_block(int plane, int block, BLOCK_SIZE_TYPE bsize,
                         int ss_txfrm_size, void *arg) {
  MACROBLOCKD* const xd = arg;
  struct macroblockd_plane *pd = &xd->plane[plane];
  int16_t* const qcoeff = BLOCK_OFFSET(pd->qcoeff, block, 16);
  const int stride = pd->dst.stride;
  const int raster_block = txfrm_block_to_raster_block(xd, bsize, plane,
                                                       block, ss_txfrm_size);
  uint8_t* const dst = raster_block_offset_uint8(xd, bsize, plane,
                                                 raster_block,
                                                 pd->dst.buf, stride);

  TX_TYPE tx_type;

  switch (ss_txfrm_size / 2) {
    case TX_4X4:
      tx_type = plane == 0 ? get_tx_type_4x4(xd, raster_block) : DCT_DCT;
      if (tx_type == DCT_DCT)
        xd->itxm_add(qcoeff, dst, stride, pd->eobs[block]);
      else
        vp9_iht_add_c(tx_type, qcoeff, dst, stride, pd->eobs[block]);
      break;
    case TX_8X8:
      tx_type = plane == 0 ? get_tx_type_8x8(xd, raster_block) : DCT_DCT;
      vp9_iht_add_8x8_c(tx_type, qcoeff, dst, stride, pd->eobs[block]);
      break;
    case TX_16X16:
      tx_type = plane == 0 ? get_tx_type_16x16(xd, raster_block) : DCT_DCT;
      vp9_iht_add_16x16_c(tx_type, qcoeff, dst, stride, pd->eobs[block]);
      break;
    case TX_32X32:
      vp9_idct_add_32x32(qcoeff, dst, stride, pd->eobs[block]);
      break;
  }
}

static void decode_block_intra(int plane, int block, BLOCK_SIZE_TYPE bsize,
                               int ss_txfrm_size, void *arg) {
  MACROBLOCKD* const xd = arg;
  struct macroblockd_plane *pd = &xd->plane[plane];

  const int raster_block = txfrm_block_to_raster_block(xd, bsize, plane,
                                                       block, ss_txfrm_size);
  uint8_t* const dst = raster_block_offset_uint8(xd, bsize, plane,
                                                 raster_block,
                                                 pd->dst.buf, pd->dst.stride);
  const TX_SIZE tx_size = (TX_SIZE)(ss_txfrm_size / 2);
  int b_mode;
  int plane_b_size;
  const int tx_ib = raster_block >> tx_size;
  const int mode = plane == 0 ? xd->mode_info_context->mbmi.mode
                              : xd->mode_info_context->mbmi.uv_mode;


  if (plane == 0 && xd->mode_info_context->mbmi.sb_type < BLOCK_SIZE_SB8X8) {
    assert(bsize == BLOCK_SIZE_SB8X8);
    b_mode = xd->mode_info_context->bmi[raster_block].as_mode.first;
  } else {
    b_mode = mode;
  }

  if (xd->mb_to_right_edge < 0 || xd->mb_to_bottom_edge < 0)
    extend_for_intra(xd, plane, block, bsize, ss_txfrm_size);

  plane_b_size = b_width_log2(bsize) - pd->subsampling_x;
  vp9_predict_intra_block(xd, tx_ib, plane_b_size, tx_size, b_mode,
                          dst, pd->dst.stride);

  // Early exit if there are no coefficients
  if (xd->mode_info_context->mbmi.mb_skip_coeff)
    return;

  decode_block(plane, block, bsize, ss_txfrm_size, arg);
}

static void decode_atom(VP9D_COMP *pbi, MACROBLOCKD *xd,
                        int mi_row, int mi_col,
                        vp9_reader *r, BLOCK_SIZE_TYPE bsize) {
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;

  assert(mbmi->ref_frame[0] != INTRA_FRAME);
  vp9_setup_interp_filters(xd, mbmi->interp_filter, &pbi->common);

  // prediction
  vp9_build_inter_predictors_sb(xd, mi_row, mi_col, bsize);

  if (mbmi->mb_skip_coeff) {
    vp9_reset_sb_tokens_context(xd, bsize);
  } else {
    if (xd->segmentation_enabled)
      mb_init_dequantizer(&pbi->common, xd);

    if (!vp9_reader_has_error(r))
      vp9_decode_tokens(pbi, r, bsize);

    foreach_transformed_block(xd, bsize, decode_block, xd);
  }
}

static void decode_sb_intra(VP9D_COMP *pbi, MACROBLOCKD *xd,
                          int mi_row, int mi_col,
                          vp9_reader *r, BLOCK_SIZE_TYPE bsize) {
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;
  if (mbmi->mb_skip_coeff) {
    vp9_reset_sb_tokens_context(xd, bsize);
  } else {
    if (xd->segmentation_enabled)
      mb_init_dequantizer(&pbi->common, xd);

    if (!vp9_reader_has_error(r))
      vp9_decode_tokens(pbi, r, bsize);
  }

  foreach_transformed_block(xd, bsize, decode_block_intra, xd);
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
  assert(mbmi->ref_frame[0] != INTRA_FRAME);

  vp9_setup_interp_filters(xd, mbmi->interp_filter, pc);

  // generate prediction
  vp9_build_inter_predictors_sb(xd, mi_row, mi_col, bsize);

  if (mbmi->mb_skip_coeff) {
    vp9_reset_sb_tokens_context(xd, bsize);
  } else {
    // re-initialize macroblock dequantizer before detokenization
    if (xd->segmentation_enabled)
      mb_init_dequantizer(pc, xd);

    // dequantization and idct
    eobtotal = vp9_decode_tokens(pbi, r, bsize);
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

static void set_offsets(VP9D_COMP *pbi, BLOCK_SIZE_TYPE bsize,
                        int mi_row, int mi_col) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  const int bh = 1 << mi_height_log2(bsize);
  const int bw = 1 << mi_width_log2(bsize);
  const int mi_idx = mi_row * cm->mode_info_stride + mi_col;
  int i;

  xd->mode_info_context = cm->mi + mi_idx;
  xd->mode_info_context->mbmi.sb_type = bsize;
  // Special case: if prev_mi is NULL, the previous mode info context
  // cannot be used.
  xd->prev_mode_info_context = cm->prev_mi ? cm->prev_mi + mi_idx : NULL;

  for (i = 0; i < MAX_MB_PLANE; i++) {
    struct macroblockd_plane *pd = &xd->plane[i];
    pd->above_context = cm->above_context[i] +
                            (mi_col * 2 >> pd->subsampling_x);
    pd->left_context = cm->left_context[i] +
                           (((mi_row * 2) & 15) >> pd->subsampling_y);
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

  // Select the appropriate reference frame for this MB
  const int fb_idx = cm->active_ref_idx[mbmi->ref_frame[0] - 1];
  const YV12_BUFFER_CONFIG *cfg = &cm->yv12_fb[fb_idx];
  xd->scale_factor[0] = cm->active_ref_scale[mbmi->ref_frame[0] - 1];
  xd->scale_factor_uv[0] = cm->active_ref_scale[mbmi->ref_frame[0] - 1];
  setup_pre_planes(xd, cfg, NULL, mi_row, mi_col, xd->scale_factor,
                   xd->scale_factor_uv);
  xd->corrupted |= cfg->corrupted;

  if (mbmi->ref_frame[1] > INTRA_FRAME) {
    // Select the appropriate reference frame for this MB
    const int second_fb_idx = cm->active_ref_idx[mbmi->ref_frame[1] - 1];
    const YV12_BUFFER_CONFIG *second_cfg = &cm->yv12_fb[second_fb_idx];
    xd->scale_factor[1] = cm->active_ref_scale[mbmi->ref_frame[1] - 1];
    xd->scale_factor_uv[1] = cm->active_ref_scale[mbmi->ref_frame[1] - 1];
    setup_pre_planes(xd, NULL, second_cfg, mi_row, mi_col, xd->scale_factor,
                     xd->scale_factor_uv);
    xd->corrupted |= second_cfg->corrupted;
  }
}

static void decode_modes_b(VP9D_COMP *pbi, int mi_row, int mi_col,
                           vp9_reader *r, BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &pbi->mb;

  if (bsize < BLOCK_SIZE_SB8X8)
    if (xd->ab_index > 0)
      return;
  set_offsets(pbi, bsize, mi_row, mi_col);
  vp9_decode_mb_mode_mv(pbi, xd, mi_row, mi_col, r);

  if (xd->mode_info_context->mbmi.ref_frame[0] == INTRA_FRAME) {
    decode_sb_intra(pbi, xd, mi_row, mi_col, r, (bsize < BLOCK_SIZE_SB8X8) ?
                                     BLOCK_SIZE_SB8X8 : bsize);
  } else {
    set_refs(pbi, mi_row, mi_col);
    if (bsize < BLOCK_SIZE_SB8X8)
      decode_atom(pbi, xd, mi_row, mi_col, r, BLOCK_SIZE_SB8X8);
    else
      decode_sb(pbi, xd, mi_row, mi_col, r, bsize);
  }
  xd->corrupted |= vp9_reader_has_error(r);
}

static void decode_modes_sb(VP9D_COMP *pbi, int mi_row, int mi_col,
                            vp9_reader* r, BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON *const pc = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  int bs = (1 << mi_width_log2(bsize)) / 2, n;
  PARTITION_TYPE partition = PARTITION_NONE;
  BLOCK_SIZE_TYPE subsize;

  if (mi_row >= pc->mi_rows || mi_col >= pc->mi_cols)
    return;

  if (bsize < BLOCK_SIZE_SB8X8)
    if (xd->ab_index != 0)
      return;

  if (bsize >= BLOCK_SIZE_SB8X8) {
    int pl;
    int idx = check_bsize_coverage(pc, xd, mi_row, mi_col, bsize);
    // read the partition information
    xd->left_seg_context = pc->left_seg_context + (mi_row & MI_MASK);
    xd->above_seg_context = pc->above_seg_context + mi_col;
    pl = partition_plane_context(xd, bsize);

    if (idx == 0)
      partition = treed_read(r, vp9_partition_tree,
                             pc->fc.partition_prob[pc->frame_type][pl]);
    else if (idx > 0 &&
        !vp9_read(r, pc->fc.partition_prob[pc->frame_type][pl][idx]))
      partition = (idx == 1) ? PARTITION_HORZ : PARTITION_VERT;
    else
      partition = PARTITION_SPLIT;

    pc->fc.partition_counts[pl][partition]++;
  }

  subsize = get_subsize(bsize, partition);
  *(get_sb_index(xd, subsize)) = 0;

  switch (partition) {
    case PARTITION_NONE:
      decode_modes_b(pbi, mi_row, mi_col, r, subsize);
      break;
    case PARTITION_HORZ:
      decode_modes_b(pbi, mi_row, mi_col, r, subsize);
      *(get_sb_index(xd, subsize)) = 1;
      if (mi_row + bs < pc->mi_rows)
        decode_modes_b(pbi, mi_row + bs, mi_col, r, subsize);
      break;
    case PARTITION_VERT:
      decode_modes_b(pbi, mi_row, mi_col, r, subsize);
      *(get_sb_index(xd, subsize)) = 1;
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
  if (bsize >= BLOCK_SIZE_SB8X8 &&
      (bsize == BLOCK_SIZE_SB8X8 || partition != PARTITION_SPLIT)) {
    set_partition_seg_context(pc, xd, mi_row, mi_col);
    update_partition_context(xd, subsize, bsize);
  }
}

static void setup_token_decoder(VP9D_COMP *pbi,
                                const uint8_t *data, size_t read_size,
                                vp9_reader *r) {
  VP9_COMMON *pc = &pbi->common;
  const uint8_t *data_end = pbi->source + pbi->source_sz;

  // Validate the calculated partition length. If the buffer
  // described by the partition can't be fully read, then restrict
  // it to the portion that can be (for EC mode) or throw an error.
  if (!read_is_valid(data, read_size, data_end))
    vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt tile length");

  if (vp9_reader_init(r, data, read_size))
    vpx_internal_error(&pc->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder %d", 1);
}

static void read_coef_probs_common(FRAME_CONTEXT *fc, TX_SIZE tx_size,
                                   vp9_reader *r) {
  vp9_coeff_probs_model *coef_probs = fc->coef_probs[tx_size];

  if (vp9_read_bit(r)) {
    int i, j, k, l, m;
    for (i = 0; i < BLOCK_TYPES; i++) {
      for (j = 0; j < REF_TYPES; j++) {
        for (k = 0; k < COEF_BANDS; k++) {
          for (l = 0; l < PREV_COEF_CONTEXTS; l++) {
            if (l >= 3 && k == 0)
              continue;

            for (m = 0; m < UNCONSTRAINED_NODES; m++) {
              vp9_prob *const p = coef_probs[i][j][k][l] + m;

              if (vp9_read(r, VP9_COEF_UPDATE_PROB))
                *p = vp9_read_prob_diff_update(r, *p);
            }
          }
        }
      }
    }
  }
}

static void read_coef_probs(VP9D_COMP *pbi, vp9_reader *r) {
  const TXFM_MODE txfm_mode = pbi->common.txfm_mode;
  FRAME_CONTEXT *const fc = &pbi->common.fc;

  read_coef_probs_common(fc, TX_4X4, r);

  if (txfm_mode > ONLY_4X4)
    read_coef_probs_common(fc, TX_8X8, r);

  if (txfm_mode > ALLOW_8X8)
    read_coef_probs_common(fc, TX_16X16, r);

  if (txfm_mode > ALLOW_16X16)
    read_coef_probs_common(fc, TX_32X32, r);
}

static void setup_segmentation(VP9D_COMP *pbi, struct vp9_read_bit_buffer *rb) {
  int i, j;

  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;

  xd->update_mb_segmentation_map = 0;
  xd->update_mb_segmentation_data = 0;

  xd->segmentation_enabled = vp9_rb_read_bit(rb);
  if (!xd->segmentation_enabled)
    return;

  // Segmentation map update
  xd->update_mb_segmentation_map = vp9_rb_read_bit(rb);
  if (xd->update_mb_segmentation_map) {
    for (i = 0; i < MB_SEG_TREE_PROBS; i++)
      xd->mb_segment_tree_probs[i] = vp9_rb_read_bit(rb) ?
                                         vp9_rb_read_literal(rb, 8) : MAX_PROB;

    cm->temporal_update = vp9_rb_read_bit(rb);
    if (cm->temporal_update) {
      for (i = 0; i < PREDICTION_PROBS; i++)
        cm->segment_pred_probs[i] = vp9_rb_read_bit(rb) ?
                                        vp9_rb_read_literal(rb, 8) : MAX_PROB;
    } else {
      for (i = 0; i < PREDICTION_PROBS; i++)
        cm->segment_pred_probs[i] = MAX_PROB;
    }
  }

  // Segmentation data update
  xd->update_mb_segmentation_data = vp9_rb_read_bit(rb);
  if (xd->update_mb_segmentation_data) {
    xd->mb_segment_abs_delta = vp9_rb_read_bit(rb);

    vp9_clearall_segfeatures(xd);

    for (i = 0; i < MAX_MB_SEGMENTS; i++) {
      for (j = 0; j < SEG_LVL_MAX; j++) {
        int data = 0;
        const int feature_enabled = vp9_rb_read_bit(rb);
        if (feature_enabled) {
          vp9_enable_segfeature(xd, i, j);
          data = decode_unsigned_max(rb, vp9_seg_feature_data_max(j));
          if (vp9_is_segfeature_signed(j))
            data = vp9_rb_read_bit(rb) ? -data : data;
        }
        vp9_set_segdata(xd, i, j, data);
      }
    }
  }
}

static void setup_loopfilter(VP9D_COMP *pbi, struct vp9_read_bit_buffer *rb) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;

  cm->filter_level = vp9_rb_read_literal(rb, 6);
  cm->sharpness_level = vp9_rb_read_literal(rb, 3);

  // Read in loop filter deltas applied at the MB level based on mode or ref
  // frame.
  xd->mode_ref_lf_delta_update = 0;

  xd->mode_ref_lf_delta_enabled = vp9_rb_read_bit(rb);
  if (xd->mode_ref_lf_delta_enabled) {
    xd->mode_ref_lf_delta_update = vp9_rb_read_bit(rb);
    if (xd->mode_ref_lf_delta_update) {
      int i;

      for (i = 0; i < MAX_REF_LF_DELTAS; i++) {
        if (vp9_rb_read_bit(rb)) {
          const int value = vp9_rb_read_literal(rb, 6);
          xd->ref_lf_deltas[i] = vp9_rb_read_bit(rb) ? -value : value;
        }
      }

      for (i = 0; i < MAX_MODE_LF_DELTAS; i++) {
        if (vp9_rb_read_bit(rb)) {
          const int value = vp9_rb_read_literal(rb, 6);
          xd->mode_lf_deltas[i] = vp9_rb_read_bit(rb) ? -value : value;
        }
      }
    }
  }
}

static int read_delta_q(struct vp9_read_bit_buffer *rb, int *delta_q) {
  const int old = *delta_q;
  if (vp9_rb_read_bit(rb)) {
    const int value = vp9_rb_read_literal(rb, 4);
    *delta_q = vp9_rb_read_bit(rb) ? -value : value;
  }
  return old != *delta_q;
}

static void setup_quantization(VP9D_COMP *pbi, struct vp9_read_bit_buffer *rb) {
  MACROBLOCKD *const xd = &pbi->mb;
  VP9_COMMON *const cm = &pbi->common;
  int update = 0;

  cm->base_qindex = vp9_rb_read_literal(rb, QINDEX_BITS);
  update |= read_delta_q(rb, &cm->y_dc_delta_q);
  update |= read_delta_q(rb, &cm->uv_dc_delta_q);
  update |= read_delta_q(rb, &cm->uv_ac_delta_q);
  if (update)
    vp9_init_dequantizer(cm);

  xd->lossless = cm->base_qindex == 0 &&
                 cm->y_dc_delta_q == 0 &&
                 cm->uv_dc_delta_q == 0 &&
                 cm->uv_ac_delta_q == 0;
  if (xd->lossless) {
    xd->itxm_add          = vp9_idct_add_lossless_c;
  } else {
    xd->itxm_add          = vp9_idct_add;
  }
}

static INTERPOLATIONFILTERTYPE read_interp_filter_type(
    struct vp9_read_bit_buffer *rb) {
  return vp9_rb_read_bit(rb) ? SWITCHABLE
                             : vp9_rb_read_literal(rb, 2);
}

static void read_frame_size(VP9_COMMON *cm, struct vp9_read_bit_buffer *rb,
                            int *width, int *height) {
  const int w = vp9_rb_read_literal(rb, 16) + 1;
  const int h = vp9_rb_read_literal(rb, 16) + 1;
  *width = w;
  *height = h;
}

static void setup_display_size(VP9D_COMP *pbi, struct vp9_read_bit_buffer *rb) {
  VP9_COMMON *const cm = &pbi->common;
  cm->display_width = cm->width;
  cm->display_height = cm->height;
  if (vp9_rb_read_bit(rb))
    read_frame_size(cm, rb, &cm->display_width, &cm->display_height);
}

static void apply_frame_size(VP9D_COMP *pbi, int width, int height) {
  VP9_COMMON *cm = &pbi->common;

  if (cm->width != width || cm->height != height) {
    if (!pbi->initial_width || !pbi->initial_height) {
      if (vp9_alloc_frame_buffers(cm, width, height))
        vpx_internal_error(&cm->error, VPX_CODEC_MEM_ERROR,
                           "Failed to allocate frame buffers");
      pbi->initial_width = width;
      pbi->initial_height = height;
    } else {
      if (width > pbi->initial_width)
        vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                           "Frame width too large");

      if (height > pbi->initial_height)
        vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                           "Frame height too large");
    }

    cm->width = width;
    cm->height = height;

    vp9_update_frame_size(cm);
  }

  vp9_realloc_frame_buffer(&cm->yv12_fb[cm->new_fb_idx], cm->width, cm->height,
                           cm->subsampling_x, cm->subsampling_y,
                           VP9BORDERINPIXELS);
}

static void setup_frame_size(VP9D_COMP *pbi,
                             struct vp9_read_bit_buffer *rb) {
  VP9_COMMON *const cm = &pbi->common;
  int width, height;
  read_frame_size(cm, rb, &width, &height);
  setup_display_size(pbi, rb);
  apply_frame_size(pbi, width, height);
}

static void setup_frame_size_with_refs(VP9D_COMP *pbi,
                                       struct vp9_read_bit_buffer *rb) {
  VP9_COMMON *const cm = &pbi->common;

  int width, height;
  int found = 0, i;
  for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i) {
    if (vp9_rb_read_bit(rb)) {
      YV12_BUFFER_CONFIG *cfg = &cm->yv12_fb[cm->active_ref_idx[i]];
      width = cfg->y_crop_width;
      height = cfg->y_crop_height;
      found = 1;
      break;
    }
  }

  if (!found)
    read_frame_size(cm, rb, &width, &height);

  if (!width || !height)
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Referenced frame with invalid size");

  setup_display_size(pbi, rb);
  apply_frame_size(pbi, width, height);
}

static void update_frame_context(FRAME_CONTEXT *fc) {
  vp9_copy(fc->pre_coef_probs, fc->coef_probs);
  vp9_copy(fc->pre_y_mode_prob, fc->y_mode_prob);
  vp9_copy(fc->pre_uv_mode_prob, fc->uv_mode_prob);
  vp9_copy(fc->pre_partition_prob, fc->partition_prob[1]);
  vp9_copy(fc->pre_intra_inter_prob, fc->intra_inter_prob);
  vp9_copy(fc->pre_comp_inter_prob, fc->comp_inter_prob);
  vp9_copy(fc->pre_single_ref_prob, fc->single_ref_prob);
  vp9_copy(fc->pre_comp_ref_prob, fc->comp_ref_prob);
  fc->pre_nmvc = fc->nmvc;
  vp9_copy(fc->pre_switchable_interp_prob, fc->switchable_interp_prob);
  vp9_copy(fc->pre_inter_mode_probs, fc->inter_mode_probs);
  vp9_copy(fc->pre_tx_probs_8x8p, fc->tx_probs_8x8p);
  vp9_copy(fc->pre_tx_probs_16x16p, fc->tx_probs_16x16p);
  vp9_copy(fc->pre_tx_probs_32x32p, fc->tx_probs_32x32p);
  vp9_copy(fc->pre_mbskip_probs, fc->mbskip_probs);

  vp9_zero(fc->coef_counts);
  vp9_zero(fc->eob_branch_counts);
  vp9_zero(fc->y_mode_counts);
  vp9_zero(fc->uv_mode_counts);
  vp9_zero(fc->NMVcount);
  vp9_zero(fc->inter_mode_counts);
  vp9_zero(fc->partition_counts);
  vp9_zero(fc->switchable_interp_count);
  vp9_zero(fc->intra_inter_count);
  vp9_zero(fc->comp_inter_count);
  vp9_zero(fc->single_ref_count);
  vp9_zero(fc->comp_ref_count);
  vp9_zero(fc->tx_count_8x8p);
  vp9_zero(fc->tx_count_16x16p);
  vp9_zero(fc->tx_count_32x32p);
  vp9_zero(fc->mbskip_count);
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

static void setup_tile_info(VP9_COMMON *cm, struct vp9_read_bit_buffer *rb) {
  int delta_log2_tiles;

  vp9_get_tile_n_bits(cm, &cm->log2_tile_columns, &delta_log2_tiles);
  while (delta_log2_tiles--) {
    if (vp9_rb_read_bit(rb)) {
      cm->log2_tile_columns++;
    } else {
      break;
    }
  }

  cm->log2_tile_rows = vp9_rb_read_bit(rb);
  if (cm->log2_tile_rows)
    cm->log2_tile_rows += vp9_rb_read_bit(rb);

  cm->tile_columns = 1 << cm->log2_tile_columns;
  cm->tile_rows    = 1 << cm->log2_tile_rows;
}

static void decode_tiles(VP9D_COMP *pbi,
                         const uint8_t *data, size_t first_partition_size,
                         vp9_reader *residual_bc) {
  VP9_COMMON *const pc = &pbi->common;

  const uint8_t *data_ptr = data + first_partition_size;
  const uint8_t* const data_end = pbi->source + pbi->source_sz;
  int tile_row, tile_col;

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
        const int size = read_be32(data_ptr2[tile_row - 1][n_cols - 1]);
        data_ptr2[tile_row - 1][n_cols - 1] += 4;
        data_ptr2[tile_row][0] = data_ptr2[tile_row - 1][n_cols - 1] + size;
      }

      for (tile_col = 1; tile_col < n_cols; tile_col++) {
        const int size = read_be32(data_ptr2[tile_row][tile_col - 1]);
        data_ptr2[tile_row][tile_col - 1] += 4;
        data_ptr2[tile_row][tile_col] =
            data_ptr2[tile_row][tile_col - 1] + size;
      }
    }

    for (tile_row = 0; tile_row < pc->tile_rows; tile_row++) {
      vp9_get_tile_row_offsets(pc, tile_row);
      for (tile_col = n_cols - 1; tile_col >= 0; tile_col--) {
        vp9_get_tile_col_offsets(pc, tile_col);
        setup_token_decoder(pbi, data_ptr2[tile_row][tile_col],
                            data_end - data_ptr2[tile_row][tile_col],
                            residual_bc);
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
        size_t size;

        vp9_get_tile_col_offsets(pc, tile_col);

        has_more = tile_col < pc->tile_columns - 1 ||
                   tile_row < pc->tile_rows - 1;
        if (has_more) {
          if (!read_is_valid(data_ptr, 4, data_end))
            vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                         "Truncated packet or corrupt tile length");

          size = read_be32(data_ptr);
          data_ptr += 4;
        } else {
          size = data_end - data_ptr;
        }

        setup_token_decoder(pbi, data_ptr, size, residual_bc);
        decode_tile(pbi, residual_bc);
        data_ptr += size;
      }
    }
  }
}

static void check_sync_code(VP9_COMMON *cm, struct vp9_read_bit_buffer *rb) {
  if (vp9_rb_read_literal(rb, 8) != SYNC_CODE_0 ||
      vp9_rb_read_literal(rb, 8) != SYNC_CODE_1 ||
      vp9_rb_read_literal(rb, 8) != SYNC_CODE_2) {
    vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                       "Invalid frame sync code");
  }
}

static void error_handler(void *data, size_t bit_offset) {
  VP9_COMMON *const cm = (VP9_COMMON *)data;
  vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME, "Truncated packet");
}

static void setup_inter_inter(VP9_COMMON *cm) {
  int i;

  cm->allow_comp_inter_inter = 0;
  for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i) {
    cm->allow_comp_inter_inter |= i > 0 &&
        cm->ref_frame_sign_bias[i + 1] != cm->ref_frame_sign_bias[1];
  }

  if (cm->allow_comp_inter_inter) {
    // which one is always-on in comp inter-inter?
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
}

#define RESERVED \
  if (vp9_rb_read_bit(rb)) \
      vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM, \
                         "Reserved bit must be unset")

static size_t read_uncompressed_header(VP9D_COMP *pbi,
                                       struct vp9_read_bit_buffer *rb) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  int i;

  cm->last_frame_type = cm->frame_type;

  if (vp9_rb_read_literal(rb, 2) != 0x2)
      vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                         "Invalid frame marker");

  cm->version = vp9_rb_read_bit(rb);
  RESERVED;

  if (vp9_rb_read_bit(rb)) {
    // show an existing frame directly
    int frame_to_show = cm->ref_frame_map[vp9_rb_read_literal(rb, 3)];
    ref_cnt_fb(cm->fb_idx_ref_cnt, &cm->new_fb_idx, frame_to_show);
    pbi->refresh_frame_flags = 0;
    cm->filter_level = 0;
    return 0;
  }

  cm->frame_type = (FRAME_TYPE) vp9_rb_read_bit(rb);
  cm->show_frame = vp9_rb_read_bit(rb);
  cm->error_resilient_mode = vp9_rb_read_bit(rb);

  if (cm->frame_type == KEY_FRAME) {
    int csp;

    check_sync_code(cm, rb);

    csp = vp9_rb_read_literal(rb, 3);  // colorspace
    if (csp != 7) {  // != sRGB
      vp9_rb_read_bit(rb);  // [16,235] (including xvycc) vs [0,255] range
      if (cm->version == 1) {
        cm->subsampling_x = vp9_rb_read_bit(rb);
        cm->subsampling_y = vp9_rb_read_bit(rb);
        vp9_rb_read_bit(rb);  // has extra plane
      } else {
        cm->subsampling_y = cm->subsampling_x = 1;
      }
    } else {
      if (cm->version == 1) {
        cm->subsampling_y = cm->subsampling_x = 0;
        vp9_rb_read_bit(rb);  // has extra plane
      } else {
        vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                           "RGB not supported in profile 0");
      }
    }

    pbi->refresh_frame_flags = (1 << NUM_REF_FRAMES) - 1;

    for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i)
      cm->active_ref_idx[i] = cm->new_fb_idx;

    setup_frame_size(pbi, rb);
  } else {
    cm->intra_only = cm->show_frame ? 0 : vp9_rb_read_bit(rb);

    cm->reset_frame_context = cm->error_resilient_mode ?
        0 : vp9_rb_read_literal(rb, 2);

    if (cm->intra_only) {
      check_sync_code(cm, rb);

      pbi->refresh_frame_flags = vp9_rb_read_literal(rb, NUM_REF_FRAMES);
      setup_frame_size(pbi, rb);
    } else {
       pbi->refresh_frame_flags = vp9_rb_read_literal(rb, NUM_REF_FRAMES);

      for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i) {
        const int ref = vp9_rb_read_literal(rb, NUM_REF_FRAMES_LG2);
        cm->active_ref_idx[i] = cm->ref_frame_map[ref];
        cm->ref_frame_sign_bias[LAST_FRAME + i] = vp9_rb_read_bit(rb);
      }

      setup_frame_size_with_refs(pbi, rb);

      xd->allow_high_precision_mv = vp9_rb_read_bit(rb);
      cm->mcomp_filter_type = read_interp_filter_type(rb);

      for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i)
        vp9_setup_scale_factors(cm, i);

      setup_inter_inter(cm);
    }
  }

  if (!cm->error_resilient_mode) {
    cm->refresh_frame_context = vp9_rb_read_bit(rb);
    cm->frame_parallel_decoding_mode = vp9_rb_read_bit(rb);
  } else {
    cm->refresh_frame_context = 0;
    cm->frame_parallel_decoding_mode = 1;
  }

  cm->frame_context_idx = vp9_rb_read_literal(rb, NUM_FRAME_CONTEXTS_LG2);

  if (cm->frame_type == KEY_FRAME || cm->error_resilient_mode || cm->intra_only)
    vp9_setup_past_independence(cm, xd);

  setup_loopfilter(pbi, rb);
  setup_quantization(pbi, rb);
  setup_segmentation(pbi, rb);

  setup_tile_info(cm, rb);

  return vp9_rb_read_literal(rb, 16);
}

int vp9_decode_frame(VP9D_COMP *pbi, const uint8_t **p_data_end) {
  int i;
  vp9_reader header_bc, residual_bc;
  VP9_COMMON *const pc = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;

  const uint8_t *data = pbi->source;
  const uint8_t *data_end = pbi->source + pbi->source_sz;

  struct vp9_read_bit_buffer rb = { data, data_end, 0,
                                    pc, error_handler };
  const size_t first_partition_size = read_uncompressed_header(pbi, &rb);
  const int keyframe = pc->frame_type == KEY_FRAME;
  YV12_BUFFER_CONFIG *new_fb = &pc->yv12_fb[pc->new_fb_idx];

  if (!first_partition_size) {
    // showing a frame directly
    *p_data_end = data + 1;
    return 0;
  }
  data += vp9_rb_bytes_read(&rb);
  xd->corrupted = 0;
  new_fb->corrupted = 0;

  if (!pbi->decoded_key_frame && !keyframe)
    return -1;

  if (!read_is_valid(data, first_partition_size, data_end))
    vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt header length");

  xd->mode_info_context = pc->mi;
  xd->prev_mode_info_context = pc->prev_mi;
  xd->frame_type = pc->frame_type;
  xd->mode_info_stride = pc->mode_info_stride;

  if (vp9_reader_init(&header_bc, data, first_partition_size))
    vpx_internal_error(&pc->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder 0");

  mb_init_dequantizer(pc, &pbi->mb);  // MB level dequantizer setup

  if (!keyframe)
    vp9_setup_interp_filters(xd, pc->mcomp_filter_type, pc);

  pc->fc = pc->frame_contexts[pc->frame_context_idx];

  update_frame_context(&pc->fc);

  setup_txfm_mode(pc, xd->lossless, &header_bc);

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

  set_prev_mi(pc);

  vp9_decode_mode_mvs_init(pbi, &header_bc);

  decode_tiles(pbi, data, first_partition_size, &residual_bc);

  pc->last_width = pc->width;
  pc->last_height = pc->height;

  new_fb->corrupted = vp9_reader_has_error(&header_bc) | xd->corrupted;

  if (!pbi->decoded_key_frame) {
    if (keyframe && !new_fb->corrupted)
      pbi->decoded_key_frame = 1;
    else
      vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                         "A stream must start with a complete key frame");
  }

  // Adaptation
  if (!pc->error_resilient_mode && !pc->frame_parallel_decoding_mode) {
    vp9_adapt_coef_probs(pc);

    if ((!keyframe) && (!pc->intra_only)) {
      vp9_adapt_mode_probs(pc);
      vp9_adapt_mode_context(pc);
      vp9_adapt_mv_probs(pc, xd->allow_high_precision_mv);
    }
  }

  if (pc->refresh_frame_context)
    pc->frame_contexts[pc->frame_context_idx] = pc->fc;

  *p_data_end = vp9_reader_find_end(&residual_bc);
  return 0;
}
