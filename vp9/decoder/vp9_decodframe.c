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

#include "vp9/common/vp9_alloccommon.h"
#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_extend.h"
#include "vp9/common/vp9_idct.h"
#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_quant_common.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_tile_common.h"

#include "vp9/decoder/vp9_dboolhuff.h"
#include "vp9/decoder/vp9_decodframe.h"
#include "vp9/decoder/vp9_detokenize.h"
#include "vp9/decoder/vp9_decodemv.h"
#include "vp9/decoder/vp9_dsubexp.h"
#include "vp9/decoder/vp9_onyxd_int.h"
#include "vp9/decoder/vp9_read_bit_buffer.h"
#include "vp9/decoder/vp9_thread.h"
#include "vp9/decoder/vp9_treereader.h"

static int read_be32(const uint8_t *p) {
  return (p[0] << 24) | (p[1] << 16) | (p[2] << 8) | p[3];
}

static int is_compound_prediction_allowed(const VP9_COMMON *cm) {
  int i;
  for (i = 1; i < ALLOWED_REFS_PER_FRAME; ++i)
    if  (cm->ref_frame_sign_bias[i + 1] != cm->ref_frame_sign_bias[1])
      return 1;

  return 0;
}

static void setup_compound_prediction(VP9_COMMON *cm) {
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

// len == 0 is not allowed
static int read_is_valid(const uint8_t *start, size_t len, const uint8_t *end) {
  return start + len > start && start + len <= end;
}

static int decode_unsigned_max(struct vp9_read_bit_buffer *rb, int max) {
  const int data = vp9_rb_read_literal(rb, get_unsigned_bits(max));
  return data > max ? max : data;
}

static TX_MODE read_tx_mode(vp9_reader *r) {
  TX_MODE tx_mode = vp9_read_literal(r, 2);
  if (tx_mode == ALLOW_32X32)
    tx_mode += vp9_read_bit(r);
  return tx_mode;
}

static void read_tx_probs(struct tx_probs *tx_probs, vp9_reader *r) {
  int i, j;

  for (i = 0; i < TX_SIZE_CONTEXTS; ++i)
    for (j = 0; j < TX_SIZES - 3; ++j)
      vp9_diff_update_prob(r, &tx_probs->p8x8[i][j]);

  for (i = 0; i < TX_SIZE_CONTEXTS; ++i)
    for (j = 0; j < TX_SIZES - 2; ++j)
      vp9_diff_update_prob(r, &tx_probs->p16x16[i][j]);

  for (i = 0; i < TX_SIZE_CONTEXTS; ++i)
    for (j = 0; j < TX_SIZES - 1; ++j)
      vp9_diff_update_prob(r, &tx_probs->p32x32[i][j]);
}

static void read_switchable_interp_probs(FRAME_CONTEXT *fc, vp9_reader *r) {
  int i, j;
  for (j = 0; j < SWITCHABLE_FILTERS + 1; ++j)
    for (i = 0; i < SWITCHABLE_FILTERS - 1; ++i)
      vp9_diff_update_prob(r, &fc->switchable_interp_prob[j][i]);
}

static void read_inter_mode_probs(FRAME_CONTEXT *fc, vp9_reader *r) {
  int i, j;
  for (i = 0; i < INTER_MODE_CONTEXTS; ++i)
    for (j = 0; j < INTER_MODES - 1; ++j)
      vp9_diff_update_prob(r, &fc->inter_mode_probs[i][j]);
}

static INLINE COMPPREDMODE_TYPE read_comp_pred_mode(vp9_reader *r) {
  COMPPREDMODE_TYPE mode = vp9_read_bit(r);
  if (mode)
    mode += vp9_read_bit(r);
  return mode;
}

static void read_comp_pred(VP9_COMMON *cm, vp9_reader *r) {
  int i;

  const int compound_allowed = is_compound_prediction_allowed(cm);
  cm->comp_pred_mode = compound_allowed ? read_comp_pred_mode(r)
                                        : SINGLE_PREDICTION_ONLY;
  if (compound_allowed)
    setup_compound_prediction(cm);

  if (cm->comp_pred_mode == HYBRID_PREDICTION)
    for (i = 0; i < COMP_INTER_CONTEXTS; i++)
      vp9_diff_update_prob(r, &cm->fc.comp_inter_prob[i]);

  if (cm->comp_pred_mode != COMP_PREDICTION_ONLY)
    for (i = 0; i < REF_CONTEXTS; i++) {
      vp9_diff_update_prob(r, &cm->fc.single_ref_prob[i][0]);
      vp9_diff_update_prob(r, &cm->fc.single_ref_prob[i][1]);
    }

  if (cm->comp_pred_mode != SINGLE_PREDICTION_ONLY)
    for (i = 0; i < REF_CONTEXTS; i++)
      vp9_diff_update_prob(r, &cm->fc.comp_ref_prob[i]);
}

static void update_mv(vp9_reader *r, vp9_prob *p) {
  if (vp9_read(r, NMV_UPDATE_PROB))
    *p = (vp9_read_literal(r, 7) << 1) | 1;
}

static void read_mv_probs(vp9_reader *r, nmv_context *mvc, int allow_hp) {
  int i, j, k;

  for (j = 0; j < MV_JOINTS - 1; ++j)
    update_mv(r, &mvc->joints[j]);

  for (i = 0; i < 2; ++i) {
    nmv_component *const comp = &mvc->comps[i];

    update_mv(r, &comp->sign);

    for (j = 0; j < MV_CLASSES - 1; ++j)
      update_mv(r, &comp->classes[j]);

    for (j = 0; j < CLASS0_SIZE - 1; ++j)
      update_mv(r, &comp->class0[j]);

    for (j = 0; j < MV_OFFSET_BITS; ++j)
      update_mv(r, &comp->bits[j]);
  }

  for (i = 0; i < 2; ++i) {
    nmv_component *const comp = &mvc->comps[i];

    for (j = 0; j < CLASS0_SIZE; ++j)
      for (k = 0; k < 3; ++k)
        update_mv(r, &comp->class0_fp[j][k]);

    for (j = 0; j < 3; ++j)
      update_mv(r, &comp->fp[j]);
  }

  if (allow_hp) {
    for (i = 0; i < 2; ++i) {
      update_mv(r, &mvc->comps[i].class0_hp);
      update_mv(r, &mvc->comps[i].hp);
    }
  }
}

static void setup_plane_dequants(VP9_COMMON *cm, MACROBLOCKD *xd, int q_index) {
  int i;
  xd->plane[0].dequant = cm->y_dequant[q_index];

  for (i = 1; i < MAX_MB_PLANE; i++)
    xd->plane[i].dequant = cm->uv_dequant[q_index];
}

// Allocate storage for each tile column.
// TODO(jzern): when max_threads <= 1 the same storage could be used for each
// tile.
static void alloc_tile_storage(VP9D_COMP *pbi, int tile_cols) {
  VP9_COMMON *const cm = &pbi->common;
  int tile_col;

  CHECK_MEM_ERROR(cm, pbi->mi_streams,
                  vpx_realloc(pbi->mi_streams, tile_cols *
                              sizeof(*pbi->mi_streams)));
  for (tile_col = 0; tile_col < tile_cols; ++tile_col) {
    vp9_get_tile_col_offsets(cm, tile_col);
    pbi->mi_streams[tile_col] =
        &cm->mi[cm->mi_rows * cm->cur_tile_mi_col_start];
  }
}

static void decode_block(int plane, int block, BLOCK_SIZE plane_bsize,
                         TX_SIZE tx_size, void *arg) {
  MACROBLOCKD* const xd = arg;
  struct macroblockd_plane *const pd = &xd->plane[plane];
  int16_t* const qcoeff = BLOCK_OFFSET(pd->qcoeff, block);
  const int stride = pd->dst.stride;
  const int eob = pd->eobs[block];
  if (eob > 0) {
    TX_TYPE tx_type;
    const int raster_block = txfrm_block_to_raster_block(plane_bsize, tx_size,
                                                         block);
    uint8_t* const dst = raster_block_offset_uint8(plane_bsize, raster_block,
                                                   pd->dst.buf, stride);
    switch (tx_size) {
      case TX_4X4:
        tx_type = get_tx_type_4x4(pd->plane_type, xd, raster_block);
        if (tx_type == DCT_DCT)
          xd->itxm_add(qcoeff, dst, stride, eob);
        else
          vp9_iht4x4_add(tx_type, qcoeff, dst, stride, eob);
        break;
      case TX_8X8:
        tx_type = get_tx_type_8x8(pd->plane_type, xd);
        vp9_iht8x8_add(tx_type, qcoeff, dst, stride, eob);
        break;
      case TX_16X16:
        tx_type = get_tx_type_16x16(pd->plane_type, xd);
        vp9_iht16x16_add(tx_type, qcoeff, dst, stride, eob);
        break;
      case TX_32X32:
        tx_type = DCT_DCT;
        vp9_idct32x32_add(qcoeff, dst, stride, eob);
        break;
      default:
        assert(!"Invalid transform size");
    }

    if (eob == 1) {
      *((int32_t *)qcoeff) = 0;
    } else {
      if (tx_type == DCT_DCT && tx_size <= TX_16X16 && eob <= 10)
        vpx_memset(qcoeff, 0, 4 * (4 << tx_size) * sizeof(qcoeff[0]));
      else
        vpx_memset(qcoeff, 0, (16 << (tx_size << 1)) * sizeof(qcoeff[0]));
    }
  }
}

static void decode_block_intra(int plane, int block, BLOCK_SIZE plane_bsize,
                               TX_SIZE tx_size, void *arg) {
  MACROBLOCKD* const xd = arg;
  struct macroblockd_plane *const pd = &xd->plane[plane];
  MODE_INFO *const mi = xd->mi_8x8[0];
  const int raster_block = txfrm_block_to_raster_block(plane_bsize, tx_size,
                                                       block);
  uint8_t* const dst = raster_block_offset_uint8(plane_bsize, raster_block,
                                                 pd->dst.buf, pd->dst.stride);
  const MB_PREDICTION_MODE mode = (plane == 0)
        ? ((mi->mbmi.sb_type < BLOCK_8X8) ? mi->bmi[raster_block].as_mode
                                          : mi->mbmi.mode)
        : mi->mbmi.uv_mode;

  if (xd->mb_to_right_edge < 0 || xd->mb_to_bottom_edge < 0)
    extend_for_intra(xd, plane_bsize, plane, block, tx_size);

  vp9_predict_intra_block(xd, raster_block >> tx_size,
                          b_width_log2(plane_bsize), tx_size, mode,
                          dst, pd->dst.stride, dst, pd->dst.stride);

  if (!mi->mbmi.skip_coeff)
    decode_block(plane, block, plane_bsize, tx_size, arg);
}

static int decode_tokens(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                         BLOCK_SIZE bsize, vp9_reader *r) {
  MB_MODE_INFO *const mbmi = &xd->mi_8x8[0]->mbmi;

  if (mbmi->skip_coeff) {
    reset_skip_context(xd, bsize);
    return -1;
  } else {
    if (cm->seg.enabled)
      setup_plane_dequants(cm, xd, vp9_get_qindex(&cm->seg, mbmi->segment_id,
                                                  cm->base_qindex));

    // TODO(dkovalev) if (!vp9_reader_has_error(r))
    return vp9_decode_tokens(cm, xd, &cm->seg, r, bsize);
  }
}

static void set_offsets(VP9D_COMP *pbi, BLOCK_SIZE bsize,
                        int mi_row, int mi_col) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  const int bh = num_8x8_blocks_high_lookup[bsize];
  const int bw = num_8x8_blocks_wide_lookup[bsize];
  const int offset = mi_row * cm->mode_info_stride + mi_col;

  xd->mode_info_stride = cm->mode_info_stride;

  xd->mi_8x8 = cm->mi_grid_visible + offset;
  xd->prev_mi_8x8 = cm->prev_mi_grid_visible + offset;

  // we are using the mode info context stream here
  xd->mi_8x8[0] = xd->mi_stream;
  xd->mi_8x8[0]->mbmi.sb_type = bsize;
  ++xd->mi_stream;

  // Special case: if prev_mi is NULL, the previous mode info context
  // cannot be used.
  xd->last_mi = cm->prev_mi ? xd->prev_mi_8x8[0] : NULL;

  set_skip_context(cm, xd, mi_row, mi_col);

  // Distance of Mb to the various image edges. These are specified to 8th pel
  // as they are always compared to values that are in 1/8th pel units
  set_mi_row_col(cm, xd, mi_row, bh, mi_col, bw);

  setup_dst_planes(xd, &cm->yv12_fb[cm->new_fb_idx], mi_row, mi_col);
}

static void set_ref(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                    int idx, int mi_row, int mi_col) {
  MB_MODE_INFO *const mbmi = &xd->mi_8x8[0]->mbmi;
  const int ref = mbmi->ref_frame[idx] - LAST_FRAME;
  const YV12_BUFFER_CONFIG *cfg = get_frame_ref_buffer(cm, ref);
  const struct scale_factors_common *sfc = &cm->active_ref_scale_comm[ref];
  if (!vp9_is_valid_scale(sfc))
    vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                       "Invalid scale factors");

  xd->scale_factor[idx].sfc = sfc;
  setup_pre_planes(xd, idx, cfg, mi_row, mi_col, &xd->scale_factor[idx]);
  xd->corrupted |= cfg->corrupted;
}

static void decode_modes_b(VP9D_COMP *pbi, int mi_row, int mi_col,
                           vp9_reader *r, BLOCK_SIZE bsize, int index) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  const int less8x8 = bsize < BLOCK_8X8;
  MB_MODE_INFO *mbmi;
  int eobtotal;

  if (less8x8)
    if (index > 0)
      return;

  set_offsets(pbi, bsize, mi_row, mi_col);
  vp9_read_mode_info(cm, xd, mi_row, mi_col, r);

  if (less8x8)
    bsize = BLOCK_8X8;

  // Has to be called after set_offsets
  mbmi = &xd->mi_8x8[0]->mbmi;
  eobtotal = decode_tokens(cm, xd, bsize, r);

  if (!is_inter_block(mbmi)) {
    // Intra reconstruction
    foreach_transformed_block(xd, bsize, decode_block_intra, xd);
  } else {
    // Inter reconstruction
    const int decode_blocks = (eobtotal > 0);

    if (!less8x8) {
      assert(mbmi->sb_type == bsize);
      if (eobtotal == 0)
        mbmi->skip_coeff = 1;  // skip loopfilter
    }

    set_ref(cm, xd, 0, mi_row, mi_col);
    if (has_second_ref(mbmi))
      set_ref(cm, xd, 1, mi_row, mi_col);

    xd->subpix.filter_x = xd->subpix.filter_y =
        vp9_get_filter_kernel(mbmi->interp_filter);
    vp9_build_inter_predictors_sb(xd, mi_row, mi_col, bsize);

    if (decode_blocks)
      foreach_transformed_block(xd, bsize, decode_block, xd);
  }
  xd->corrupted |= vp9_reader_has_error(r);
}

static void decode_modes_sb(VP9D_COMP *pbi, int mi_row, int mi_col,
                            vp9_reader* r, BLOCK_SIZE bsize, int index) {
  VP9_COMMON *const cm = &pbi->common;
  const int hbs = num_8x8_blocks_wide_lookup[bsize] / 2;
  PARTITION_TYPE partition = PARTITION_NONE;
  BLOCK_SIZE subsize;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  if (bsize < BLOCK_8X8) {
    if (index > 0)
      return;
  } else {
    int pl;
    const int idx = check_bsize_coverage(hbs, cm->mi_rows, cm->mi_cols,
                                         mi_row, mi_col);
    pl = partition_plane_context(cm, mi_row, mi_col, bsize);

    if (idx == 0)
      partition = treed_read(r, vp9_partition_tree,
                             cm->fc.partition_prob[cm->frame_type][pl]);
    else if (idx > 0 &&
        !vp9_read(r, cm->fc.partition_prob[cm->frame_type][pl][idx]))
      partition = (idx == 1) ? PARTITION_HORZ : PARTITION_VERT;
    else
      partition = PARTITION_SPLIT;

    if (!cm->frame_parallel_decoding_mode)
      ++cm->counts.partition[pl][partition];
  }

  subsize = get_subsize(bsize, partition);

  switch (partition) {
    case PARTITION_NONE:
      decode_modes_b(pbi, mi_row, mi_col, r, subsize, 0);
      break;
    case PARTITION_HORZ:
      decode_modes_b(pbi, mi_row, mi_col, r, subsize, 0);
      if (mi_row + hbs < cm->mi_rows)
        decode_modes_b(pbi, mi_row + hbs, mi_col, r, subsize, 1);
      break;
    case PARTITION_VERT:
      decode_modes_b(pbi, mi_row, mi_col, r, subsize, 0);
      if (mi_col + hbs < cm->mi_cols)
        decode_modes_b(pbi, mi_row, mi_col + hbs, r, subsize, 1);
      break;
    case PARTITION_SPLIT: {
      int n;
      for (n = 0; n < 4; n++) {
        const int j = n >> 1, i = n & 1;
        decode_modes_sb(pbi, mi_row + j * hbs, mi_col + i * hbs,
                        r, subsize, n);
      }
    } break;
    default:
      assert(!"Invalid partition type");
  }

  // update partition context
  if (bsize >= BLOCK_8X8 &&
      (bsize == BLOCK_8X8 || partition != PARTITION_SPLIT))
    update_partition_context(cm, mi_row, mi_col, subsize, bsize);
}

static void setup_token_decoder(const uint8_t *data,
                                const uint8_t *data_end,
                                size_t read_size,
                                struct vpx_internal_error_info *error_info,
                                vp9_reader *r) {
  // Validate the calculated partition length. If the buffer
  // described by the partition can't be fully read, then restrict
  // it to the portion that can be (for EC mode) or throw an error.
  if (!read_is_valid(data, read_size, data_end))
    vpx_internal_error(error_info, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt tile length");

  if (vp9_reader_init(r, data, read_size))
    vpx_internal_error(error_info, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder %d", 1);
}

static void read_coef_probs_common(vp9_coeff_probs_model *coef_probs,
                                   vp9_reader *r) {
  int i, j, k, l, m;

  if (vp9_read_bit(r))
    for (i = 0; i < BLOCK_TYPES; i++)
      for (j = 0; j < REF_TYPES; j++)
        for (k = 0; k < COEF_BANDS; k++)
          for (l = 0; l < PREV_COEF_CONTEXTS; l++)
            if (k > 0 || l < 3)
              for (m = 0; m < UNCONSTRAINED_NODES; m++)
                vp9_diff_update_prob(r, &coef_probs[i][j][k][l][m]);
}

static void read_coef_probs(FRAME_CONTEXT *fc, TX_MODE tx_mode,
                            vp9_reader *r) {
  read_coef_probs_common(fc->coef_probs[TX_4X4], r);

  if (tx_mode > ONLY_4X4)
    read_coef_probs_common(fc->coef_probs[TX_8X8], r);

  if (tx_mode > ALLOW_8X8)
    read_coef_probs_common(fc->coef_probs[TX_16X16], r);

  if (tx_mode > ALLOW_16X16)
    read_coef_probs_common(fc->coef_probs[TX_32X32], r);
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

  xd->itxm_add = xd->lossless ? vp9_iwht4x4_add : vp9_idct4x4_add;
}

static INTERPOLATION_TYPE read_interp_filter_type(
                              struct vp9_read_bit_buffer *rb) {
  const INTERPOLATION_TYPE literal_to_type[] = { EIGHTTAP_SMOOTH,
                                                 EIGHTTAP,
                                                 EIGHTTAP_SHARP,
                                                 BILINEAR };
  return vp9_rb_read_bit(rb) ? SWITCHABLE
                             : literal_to_type[vp9_rb_read_literal(rb, 2)];
}

static void read_frame_size(struct vp9_read_bit_buffer *rb,
                            int *width, int *height) {
  const int w = vp9_rb_read_literal(rb, 16) + 1;
  const int h = vp9_rb_read_literal(rb, 16) + 1;
  *width = w;
  *height = h;
}

static void setup_display_size(VP9_COMMON *cm, struct vp9_read_bit_buffer *rb) {
  cm->display_width = cm->width;
  cm->display_height = cm->height;
  if (vp9_rb_read_bit(rb))
    read_frame_size(rb, &cm->display_width, &cm->display_height);
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
  int width, height;
  read_frame_size(rb, &width, &height);
  apply_frame_size(pbi, width, height);
  setup_display_size(&pbi->common, rb);
}

static void setup_frame_size_with_refs(VP9D_COMP *pbi,
                                       struct vp9_read_bit_buffer *rb) {
  VP9_COMMON *const cm = &pbi->common;

  int width, height;
  int found = 0, i;
  for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i) {
    if (vp9_rb_read_bit(rb)) {
      YV12_BUFFER_CONFIG *const cfg = get_frame_ref_buffer(cm, i);
      width = cfg->y_crop_width;
      height = cfg->y_crop_height;
      found = 1;
      break;
    }
  }

  if (!found)
    read_frame_size(rb, &width, &height);

  if (!width || !height)
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Referenced frame with invalid size");

  apply_frame_size(pbi, width, height);
  setup_display_size(cm, rb);
}

static void decode_tile(VP9D_COMP *pbi, vp9_reader *r, int tile_col) {
  const int num_threads = pbi->oxcf.max_threads;
  VP9_COMMON *const cm = &pbi->common;
  int mi_row, mi_col;
  YV12_BUFFER_CONFIG *const fb = &cm->yv12_fb[cm->new_fb_idx];
  MACROBLOCKD *xd = &pbi->mb;

  xd->mi_stream = pbi->mi_streams[tile_col];

  if (pbi->do_loopfilter_inline) {
    LFWorkerData *const lf_data = (LFWorkerData*)pbi->lf_worker.data1;
    lf_data->frame_buffer = fb;
    lf_data->cm = cm;
    lf_data->xd = pbi->mb;
    lf_data->stop = 0;
    lf_data->y_only = 0;
    vp9_loop_filter_frame_init(cm, cm->lf.filter_level);
  }

  for (mi_row = cm->cur_tile_mi_row_start; mi_row < cm->cur_tile_mi_row_end;
       mi_row += MI_BLOCK_SIZE) {
    // For a SB there are 2 left contexts, each pertaining to a MB row within
    vp9_zero(cm->left_context);
    vp9_zero(cm->left_seg_context);
    for (mi_col = cm->cur_tile_mi_col_start; mi_col < cm->cur_tile_mi_col_end;
         mi_col += MI_BLOCK_SIZE)
      decode_modes_sb(pbi, mi_row, mi_col, r, BLOCK_64X64, 0);

    if (pbi->do_loopfilter_inline) {
      const int lf_start = mi_row - MI_BLOCK_SIZE;
      LFWorkerData *const lf_data = (LFWorkerData*)pbi->lf_worker.data1;

      // delay the loopfilter by 1 macroblock row.
      if (lf_start < 0) continue;

      // decoding has completed: finish up the loop filter in this thread.
      if (mi_row + MI_BLOCK_SIZE >= cm->cur_tile_mi_row_end) continue;

      vp9_worker_sync(&pbi->lf_worker);
      lf_data->start = lf_start;
      lf_data->stop = mi_row;
      if (num_threads > 1) {
        vp9_worker_launch(&pbi->lf_worker);
      } else {
        vp9_worker_execute(&pbi->lf_worker);
      }
    }
  }

  if (pbi->do_loopfilter_inline) {
    LFWorkerData *const lf_data = (LFWorkerData*)pbi->lf_worker.data1;

    vp9_worker_sync(&pbi->lf_worker);
    lf_data->start = lf_data->stop;
    lf_data->stop = cm->mi_rows;
    vp9_worker_execute(&pbi->lf_worker);
  }
}

static void setup_tile_info(VP9_COMMON *cm, struct vp9_read_bit_buffer *rb) {
  int min_log2_tile_cols, max_log2_tile_cols, max_ones;
  vp9_get_tile_n_bits(cm->mi_cols, &min_log2_tile_cols, &max_log2_tile_cols);

  // columns
  max_ones = max_log2_tile_cols - min_log2_tile_cols;
  cm->log2_tile_cols = min_log2_tile_cols;
  while (max_ones-- && vp9_rb_read_bit(rb))
    cm->log2_tile_cols++;

  // rows
  cm->log2_tile_rows = vp9_rb_read_bit(rb);
  if (cm->log2_tile_rows)
    cm->log2_tile_rows += vp9_rb_read_bit(rb);
}

static const uint8_t *decode_tiles(VP9D_COMP *pbi, const uint8_t *data) {
  vp9_reader residual_bc;

  VP9_COMMON *const cm = &pbi->common;

  const uint8_t *const data_end = pbi->source + pbi->source_sz;
  const int aligned_mi_cols = mi_cols_aligned_to_sb(cm->mi_cols);
  const int tile_cols = 1 << cm->log2_tile_cols;
  const int tile_rows = 1 << cm->log2_tile_rows;
  int tile_row, tile_col;

  // Note: this memset assumes above_context[0], [1] and [2]
  // are allocated as part of the same buffer.
  vpx_memset(cm->above_context[0], 0,
             sizeof(ENTROPY_CONTEXT) * MAX_MB_PLANE * (2 * aligned_mi_cols));

  vpx_memset(cm->above_seg_context, 0,
             sizeof(PARTITION_CONTEXT) * aligned_mi_cols);

  if (pbi->oxcf.inv_tile_order) {
    const uint8_t *data_ptr2[4][1 << 6];
    vp9_reader bc_bak = {0};

    // pre-initialize the offsets, we're going to read in inverse order
    data_ptr2[0][0] = data;
    for (tile_row = 0; tile_row < tile_rows; tile_row++) {
      if (tile_row) {
        const int size = read_be32(data_ptr2[tile_row - 1][tile_cols - 1]);
        data_ptr2[tile_row - 1][tile_cols - 1] += 4;
        data_ptr2[tile_row][0] = data_ptr2[tile_row - 1][tile_cols - 1] + size;
      }

      for (tile_col = 1; tile_col < tile_cols; tile_col++) {
        const int size = read_be32(data_ptr2[tile_row][tile_col - 1]);
        data_ptr2[tile_row][tile_col - 1] += 4;
        data_ptr2[tile_row][tile_col] =
            data_ptr2[tile_row][tile_col - 1] + size;
      }
    }

    for (tile_row = 0; tile_row < tile_rows; tile_row++) {
      vp9_get_tile_row_offsets(cm, tile_row);
      for (tile_col = tile_cols - 1; tile_col >= 0; tile_col--) {
        vp9_get_tile_col_offsets(cm, tile_col);
        setup_token_decoder(data_ptr2[tile_row][tile_col], data_end,
                            data_end - data_ptr2[tile_row][tile_col],
                            &cm->error, &residual_bc);
        decode_tile(pbi, &residual_bc, tile_col);
        if (tile_row == tile_rows - 1 && tile_col == tile_cols - 1)
          bc_bak = residual_bc;
      }
    }
    residual_bc = bc_bak;
  } else {
    int has_more;

    for (tile_row = 0; tile_row < tile_rows; tile_row++) {
      vp9_get_tile_row_offsets(cm, tile_row);
      for (tile_col = 0; tile_col < tile_cols; tile_col++) {
        size_t size;

        vp9_get_tile_col_offsets(cm, tile_col);

        has_more = tile_col < tile_cols - 1 || tile_row < tile_rows - 1;
        if (has_more) {
          if (!read_is_valid(data, 4, data_end))
            vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                         "Truncated packet or corrupt tile length");

          size = read_be32(data);
          data += 4;
        } else {
          size = data_end - data;
        }

        setup_token_decoder(data, data_end, size, &cm->error, &residual_bc);
        decode_tile(pbi, &residual_bc, tile_col);
        data += size;
      }
    }
  }

  return vp9_reader_find_end(&residual_bc);
}

static void check_sync_code(VP9_COMMON *cm, struct vp9_read_bit_buffer *rb) {
  if (vp9_rb_read_literal(rb, 8) != VP9_SYNC_CODE_0 ||
      vp9_rb_read_literal(rb, 8) != VP9_SYNC_CODE_1 ||
      vp9_rb_read_literal(rb, 8) != VP9_SYNC_CODE_2) {
    vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                       "Invalid frame sync code");
  }
}

static void error_handler(void *data, size_t bit_offset) {
  VP9_COMMON *const cm = (VP9_COMMON *)data;
  vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME, "Truncated packet");
}

#define RESERVED \
  if (vp9_rb_read_bit(rb)) \
      vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM, \
                         "Reserved bit must be unset")

static size_t read_uncompressed_header(VP9D_COMP *pbi,
                                       struct vp9_read_bit_buffer *rb) {
  VP9_COMMON *const cm = &pbi->common;
  size_t sz;
  int i;

  cm->last_frame_type = cm->frame_type;

  if (vp9_rb_read_literal(rb, 2) != VP9_FRAME_MARKER)
      vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                         "Invalid frame marker");

  cm->version = vp9_rb_read_bit(rb);
  RESERVED;

  if (vp9_rb_read_bit(rb)) {
    // show an existing frame directly
    int frame_to_show = cm->ref_frame_map[vp9_rb_read_literal(rb, 3)];
    ref_cnt_fb(cm->fb_idx_ref_cnt, &cm->new_fb_idx, frame_to_show);
    pbi->refresh_frame_flags = 0;
    cm->lf.filter_level = 0;
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
        const int ref = vp9_rb_read_literal(rb, NUM_REF_FRAMES_LOG2);
        cm->active_ref_idx[i] = cm->ref_frame_map[ref];
        cm->ref_frame_sign_bias[LAST_FRAME + i] = vp9_rb_read_bit(rb);
      }

      setup_frame_size_with_refs(pbi, rb);

      cm->allow_high_precision_mv = vp9_rb_read_bit(rb);
      cm->mcomp_filter_type = read_interp_filter_type(rb);

      for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i)
        vp9_setup_scale_factors(cm, i);
    }
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
  cm->frame_context_idx = vp9_rb_read_literal(rb, NUM_FRAME_CONTEXTS_LOG2);

  if (frame_is_intra_only(cm) || cm->error_resilient_mode)
    vp9_setup_past_independence(cm);

  setup_loopfilter(&cm->lf, rb);
  setup_quantization(pbi, rb);
  setup_segmentation(&cm->seg, rb);

  setup_tile_info(cm, rb);
  sz = vp9_rb_read_literal(rb, 16);

  return sz > 0 ? sz : -1;
}

static int read_compressed_header(VP9D_COMP *pbi, const uint8_t *data,
                                  size_t partition_size) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  FRAME_CONTEXT *const fc = &cm->fc;
  vp9_reader r;
  int k;

  if (vp9_reader_init(&r, data, partition_size))
    vpx_internal_error(&cm->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder 0");

  cm->tx_mode = xd->lossless ? ONLY_4X4 : read_tx_mode(&r);
  if (cm->tx_mode == TX_MODE_SELECT)
    read_tx_probs(&fc->tx_probs, &r);
  read_coef_probs(fc, cm->tx_mode, &r);

  for (k = 0; k < MBSKIP_CONTEXTS; ++k)
    vp9_diff_update_prob(&r, &fc->mbskip_probs[k]);

  if (!frame_is_intra_only(cm)) {
    nmv_context *const nmvc = &fc->nmvc;
    int i, j;

    read_inter_mode_probs(fc, &r);

    if (cm->mcomp_filter_type == SWITCHABLE)
      read_switchable_interp_probs(fc, &r);

    for (i = 0; i < INTRA_INTER_CONTEXTS; i++)
      vp9_diff_update_prob(&r, &fc->intra_inter_prob[i]);

    read_comp_pred(cm, &r);

    for (j = 0; j < BLOCK_SIZE_GROUPS; j++)
      for (i = 0; i < INTRA_MODES - 1; ++i)
        vp9_diff_update_prob(&r, &fc->y_mode_prob[j][i]);

    for (j = 0; j < PARTITION_CONTEXTS; ++j)
      for (i = 0; i < PARTITION_TYPES - 1; ++i)
        vp9_diff_update_prob(&r, &fc->partition_prob[INTER_FRAME][j][i]);

    read_mv_probs(&r, nmvc, cm->allow_high_precision_mv);
  }

  return vp9_reader_has_error(&r);
}

void vp9_init_dequantizer(VP9_COMMON *cm) {
  int q;

  for (q = 0; q < QINDEX_RANGE; q++) {
    cm->y_dequant[q][0] = vp9_dc_quant(q, cm->y_dc_delta_q);
    cm->y_dequant[q][1] = vp9_ac_quant(q, 0);

    cm->uv_dequant[q][0] = vp9_dc_quant(q, cm->uv_dc_delta_q);
    cm->uv_dequant[q][1] = vp9_ac_quant(q, cm->uv_ac_delta_q);
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
  assert(!memcmp(cm->counts.mbskip, zero_counts.mbskip,
                 sizeof(cm->counts.mbskip)));
  assert(!memcmp(&cm->counts.mv, &zero_counts.mv, sizeof(cm->counts.mv)));
}
#endif  // NDEBUG

int vp9_decode_frame(VP9D_COMP *pbi, const uint8_t **p_data_end) {
  int i;
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;

  const uint8_t *data = pbi->source;
  const uint8_t *data_end = pbi->source + pbi->source_sz;

  struct vp9_read_bit_buffer rb = { data, data_end, 0,
                                    cm, error_handler };
  const size_t first_partition_size = read_uncompressed_header(pbi, &rb);
  const int keyframe = cm->frame_type == KEY_FRAME;
  YV12_BUFFER_CONFIG *new_fb = &cm->yv12_fb[cm->new_fb_idx];
  const int tile_cols = 1 << cm->log2_tile_cols;

  if (!first_partition_size) {
    if (!keyframe) {
      // showing a frame directly
      *p_data_end = data + 1;
      return 0;
    } else {
      vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                         "Invalid key frame");
      return -1;
    }
  }
  data += vp9_rb_bytes_read(&rb);
  xd->corrupted = 0;
  new_fb->corrupted = 0;
  pbi->do_loopfilter_inline =
      (cm->log2_tile_rows | cm->log2_tile_cols) == 0 && cm->lf.filter_level;

  if (!pbi->decoded_key_frame && !keyframe)
    return -1;

  if (!read_is_valid(data, first_partition_size, data_end))
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt header length");

  setup_plane_dequants(cm, &pbi->mb, cm->base_qindex);

  xd->mi_8x8 = cm->mi_grid_visible;
  xd->mode_info_stride = cm->mode_info_stride;

  alloc_tile_storage(pbi, tile_cols);

  cm->fc = cm->frame_contexts[cm->frame_context_idx];

  vp9_zero(cm->counts);

  new_fb->corrupted |= read_compressed_header(pbi, data, first_partition_size);

  setup_block_dptrs(xd, cm->subsampling_x, cm->subsampling_y);

  // clear out the coeff buffer
  for (i = 0; i < MAX_MB_PLANE; ++i)
    vp9_zero(xd->plane[i].qcoeff);

  set_prev_mi(cm);

  *p_data_end = decode_tiles(pbi, data + first_partition_size);

  cm->last_width = cm->width;
  cm->last_height = cm->height;

  new_fb->corrupted |= xd->corrupted;

  if (!pbi->decoded_key_frame) {
    if (keyframe && !new_fb->corrupted)
      pbi->decoded_key_frame = 1;
    else
      vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                         "A stream must start with a complete key frame");
  }

  if (!cm->error_resilient_mode && !cm->frame_parallel_decoding_mode) {
    vp9_adapt_coef_probs(cm);

    if (!frame_is_intra_only(cm)) {
      vp9_adapt_mode_probs(cm);
      vp9_adapt_mv_probs(cm, cm->allow_high_precision_mv);
    }
  } else {
    debug_check_frame_counts(cm);
  }

  if (cm->refresh_frame_context)
    cm->frame_contexts[cm->frame_context_idx] = cm->fc;

  return 0;
}
