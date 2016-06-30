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

#include "./vp10_rtcd.h"
#include "./vpx_dsp_rtcd.h"
#include "./vpx_scale_rtcd.h"

#include "vpx_dsp/bitreader_buffer.h"
#include "vp10/decoder/bitreader.h"
#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem.h"
#include "vpx_ports/mem_ops.h"
#include "vpx_scale/vpx_scale.h"
#include "vpx_util/vpx_thread.h"

#include "vp10/common/alloccommon.h"
#include "vp10/common/common.h"
#include "vp10/common/entropy.h"
#include "vp10/common/entropymode.h"
#include "vp10/common/idct.h"
#include "vp10/common/thread_common.h"
#include "vp10/common/pred_common.h"
#include "vp10/common/quant_common.h"
#include "vp10/common/reconintra.h"
#include "vp10/common/reconinter.h"
#include "vp10/common/seg_common.h"
#include "vp10/common/tile_common.h"

#include "vp10/decoder/decodeframe.h"
#include "vp10/decoder/detokenize.h"
#include "vp10/decoder/decodemv.h"
#include "vp10/decoder/decoder.h"
#include "vp10/decoder/dsubexp.h"

#define MAX_VP9_HEADER_SIZE 80

static int is_compound_reference_allowed(const VP10_COMMON *cm) {
  int i;
  if (frame_is_intra_only(cm))
    return 0;
  for (i = 1; i < REFS_PER_FRAME; ++i)
    if (cm->ref_frame_sign_bias[i + 1] != cm->ref_frame_sign_bias[1])
      return 1;

  return 0;
}

static void setup_compound_reference_mode(VP10_COMMON *cm) {
#if CONFIG_EXT_REFS
  cm->comp_fwd_ref[0] = LAST_FRAME;
  cm->comp_fwd_ref[1] = LAST2_FRAME;
  cm->comp_fwd_ref[2] = LAST3_FRAME;
  cm->comp_fwd_ref[3] = GOLDEN_FRAME;

  cm->comp_bwd_ref[0] = BWDREF_FRAME;
  cm->comp_bwd_ref[1] = ALTREF_FRAME;
#else
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
#endif  // CONFIG_EXT_REFS
}

static int read_is_valid(const uint8_t *start, size_t len, const uint8_t *end) {
  return len != 0 && len <= (size_t)(end - start);
}

static int decode_unsigned_max(struct vpx_read_bit_buffer *rb, int max) {
  const int data = vpx_rb_read_literal(rb, get_unsigned_bits(max));
  return data > max ? max : data;
}

static TX_MODE read_tx_mode(struct vpx_read_bit_buffer *rb) {
  return vpx_rb_read_bit(rb) ? TX_MODE_SELECT : vpx_rb_read_literal(rb, 2);
}

static void read_switchable_interp_probs(FRAME_CONTEXT *fc, vp10_reader *r) {
  int i, j;
  for (j = 0; j < SWITCHABLE_FILTER_CONTEXTS; ++j)
    for (i = 0; i < SWITCHABLE_FILTERS - 1; ++i)
      vp10_diff_update_prob(r, &fc->switchable_interp_prob[j][i]);
}

static void read_inter_mode_probs(FRAME_CONTEXT *fc, vp10_reader *r) {
  int i;
#if CONFIG_REF_MV
  for (i = 0; i < NEWMV_MODE_CONTEXTS; ++i)
    vp10_diff_update_prob(r, &fc->newmv_prob[i]);
  for (i = 0; i < ZEROMV_MODE_CONTEXTS; ++i)
    vp10_diff_update_prob(r, &fc->zeromv_prob[i]);
  for (i = 0; i < REFMV_MODE_CONTEXTS; ++i)
    vp10_diff_update_prob(r, &fc->refmv_prob[i]);
  for (i = 0; i < DRL_MODE_CONTEXTS; ++i)
    vp10_diff_update_prob(r, &fc->drl_prob[i]);
#if CONFIG_EXT_INTER
  vp10_diff_update_prob(r, &fc->new2mv_prob);
#endif  // CONFIG_EXT_INTER
#else
  int j;
  for (i = 0; i < INTER_MODE_CONTEXTS; ++i)
    for (j = 0; j < INTER_MODES - 1; ++j)
      vp10_diff_update_prob(r, &fc->inter_mode_probs[i][j]);
#endif
}

#if CONFIG_EXT_INTER
static void read_inter_compound_mode_probs(FRAME_CONTEXT *fc, vp10_reader *r) {
  int i, j;
  if (vp10_read(r, GROUP_DIFF_UPDATE_PROB)) {
    for (j = 0; j < INTER_MODE_CONTEXTS; ++j) {
      for (i = 0; i < INTER_COMPOUND_MODES - 1; ++i) {
        vp10_diff_update_prob(r, &fc->inter_compound_mode_probs[j][i]);
      }
    }
  }
}
#endif  // CONFIG_EXT_INTER

static REFERENCE_MODE read_frame_reference_mode(const VP10_COMMON *cm,
    struct vpx_read_bit_buffer *rb) {
  if (is_compound_reference_allowed(cm)) {
    return vpx_rb_read_bit(rb) ? REFERENCE_MODE_SELECT
                               : (vpx_rb_read_bit(rb) ? COMPOUND_REFERENCE
                                                      : SINGLE_REFERENCE);
  } else {
    return SINGLE_REFERENCE;
  }
}

static void read_frame_reference_mode_probs(VP10_COMMON *cm, vp10_reader *r) {
  FRAME_CONTEXT *const fc = cm->fc;
  int i, j;

  if (cm->reference_mode == REFERENCE_MODE_SELECT)
    for (i = 0; i < COMP_INTER_CONTEXTS; ++i)
      vp10_diff_update_prob(r, &fc->comp_inter_prob[i]);

  if (cm->reference_mode != COMPOUND_REFERENCE) {
    for (i = 0; i < REF_CONTEXTS; ++i) {
      for (j = 0; j < (SINGLE_REFS - 1); ++j) {
        vp10_diff_update_prob(r, &fc->single_ref_prob[i][j]);
      }
    }
  }

  if (cm->reference_mode != SINGLE_REFERENCE) {
    for (i = 0; i < REF_CONTEXTS; ++i) {
#if CONFIG_EXT_REFS
      for (j = 0; j < (FWD_REFS - 1); ++j)
        vp10_diff_update_prob(r, &fc->comp_ref_prob[i][j]);
      for (j = 0; j < (BWD_REFS - 1); ++j)
        vp10_diff_update_prob(r, &fc->comp_bwdref_prob[i][j]);
#else
      for (j = 0; j < (COMP_REFS - 1); ++j)
        vp10_diff_update_prob(r, &fc->comp_ref_prob[i][j]);
#endif  // CONFIG_EXT_REFS
    }
  }
}

static void update_mv_probs(vpx_prob *p, int n, vp10_reader *r) {
  int i;
  for (i = 0; i < n; ++i)
    vp10_diff_update_prob(r, &p[i]);
}

static void read_mv_probs(nmv_context *ctx, int allow_hp, vp10_reader *r) {
  int i, j;

  update_mv_probs(ctx->joints, MV_JOINTS - 1, r);

#if CONFIG_REF_MV
  vp10_diff_update_prob(r, &ctx->zero_rmv);
#endif

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

static void inverse_transform_block(MACROBLOCKD* xd, int plane,
                                    const TX_TYPE tx_type,
                                    const TX_SIZE tx_size,
                                    uint8_t *dst, int stride,
                                    int eob) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  if (eob > 0) {
    tran_low_t *const dqcoeff = pd->dqcoeff;
    INV_TXFM_PARAM inv_txfm_param;
    inv_txfm_param.tx_type = tx_type;
    inv_txfm_param.tx_size = tx_size;
    inv_txfm_param.eob = eob;
    inv_txfm_param.lossless = xd->lossless[xd->mi[0]->mbmi.segment_id];

#if CONFIG_VP9_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
      inv_txfm_param.bd = xd->bd;
      highbd_inv_txfm_add(dqcoeff, dst, stride, &inv_txfm_param);
    } else {
#endif  // CONFIG_VP9_HIGHBITDEPTH
      inv_txfm_add(dqcoeff, dst, stride, &inv_txfm_param);
#if CONFIG_VP9_HIGHBITDEPTH
    }
#endif  // CONFIG_VP9_HIGHBITDEPTH

    if (eob == 1) {
      dqcoeff[0] = 0;
    } else {
      if (tx_type == DCT_DCT && tx_size <= TX_16X16 && eob <= 10)
        memset(dqcoeff, 0, 4 * (4 << tx_size) * sizeof(dqcoeff[0]));
#if CONFIG_EXT_TX
      else
        memset(dqcoeff, 0, (16 << (tx_size << 1)) * sizeof(dqcoeff[0]));
#else
      else if (tx_size == TX_32X32 && eob <= 34)
        memset(dqcoeff, 0, 256 * sizeof(dqcoeff[0]));
      else
        memset(dqcoeff, 0, (16 << (tx_size << 1)) * sizeof(dqcoeff[0]));
#endif
    }
  }
}

static void predict_and_reconstruct_intra_block(MACROBLOCKD *const xd,
#if CONFIG_ANS
                                                struct AnsDecoder *const r,
#else
                                                vp10_reader *r,
#endif  // CONFIG_ANS
                                                MB_MODE_INFO *const mbmi,
                                                int plane,
                                                int row, int col,
                                                TX_SIZE tx_size) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  PREDICTION_MODE mode = (plane == 0) ? mbmi->mode : mbmi->uv_mode;
  PLANE_TYPE plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
  uint8_t *dst;
  int block_idx = (row << 1) + col;
  dst = &pd->dst.buf[4 * row * pd->dst.stride + 4 * col];

  if (mbmi->sb_type < BLOCK_8X8)
    if (plane == 0)
      mode = xd->mi[0]->bmi[(row << 1) + col].as_mode;

  vp10_predict_intra_block(xd, pd->n4_wl, pd->n4_hl, tx_size, mode,
                          dst, pd->dst.stride, dst, pd->dst.stride,
                          col, row, plane);

  if (!mbmi->skip) {
    TX_TYPE tx_type = get_tx_type(plane_type, xd, block_idx, tx_size);
    const scan_order *sc = get_scan(tx_size, tx_type, 0);
    const int eob = vp10_decode_block_tokens(xd,
                                             plane, sc, col, row, tx_size,
                                             tx_type, r, mbmi->segment_id);
    inverse_transform_block(xd, plane, tx_type, tx_size,
                            dst, pd->dst.stride, eob);
  }
}

#if CONFIG_VAR_TX
static void decode_reconstruct_tx(MACROBLOCKD *const xd, vp10_reader *r,
                                  MB_MODE_INFO *const mbmi,
                                  int plane, BLOCK_SIZE plane_bsize,
                                  int block, int blk_row, int blk_col,
                                  TX_SIZE tx_size, int *eob_total) {
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const BLOCK_SIZE bsize = txsize_to_bsize[tx_size];
  const int tx_row = blk_row >> (1 - pd->subsampling_y);
  const int tx_col = blk_col >> (1 - pd->subsampling_x);
  const TX_SIZE plane_tx_size = plane ?
      get_uv_tx_size_impl(mbmi->inter_tx_size[tx_row][tx_col], bsize, 0, 0) :
      mbmi->inter_tx_size[tx_row][tx_col];
  int max_blocks_high = num_4x4_blocks_high_lookup[plane_bsize];
  int max_blocks_wide = num_4x4_blocks_wide_lookup[plane_bsize];

  if (xd->mb_to_bottom_edge < 0)
    max_blocks_high += xd->mb_to_bottom_edge >> (5 + pd->subsampling_y);
  if (xd->mb_to_right_edge < 0)
    max_blocks_wide += xd->mb_to_right_edge >> (5 + pd->subsampling_x);

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide)
    return;

  if (tx_size == plane_tx_size) {
    PLANE_TYPE plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
    TX_TYPE tx_type = get_tx_type(plane_type, xd, block, tx_size);
    const scan_order *sc = get_scan(tx_size, tx_type, 1);
    const int eob = vp10_decode_block_tokens(xd, plane, sc,
                                             blk_col, blk_row, tx_size,
                                             tx_type, r, mbmi->segment_id);
    inverse_transform_block(xd, plane, tx_type, tx_size,
        &pd->dst.buf[4 * blk_row * pd->dst.stride + 4 * blk_col],
        pd->dst.stride, eob);
    *eob_total += eob;
  } else {
    int bsl = b_width_log2_lookup[bsize];
    int i;

    assert(bsl > 0);
    --bsl;

    for (i = 0; i < 4; ++i) {
      const int offsetr = blk_row + ((i >> 1) << bsl);
      const int offsetc = blk_col + ((i & 0x01) << bsl);
      int step = 1 << (2 * (tx_size - 1));

      if (offsetr >= max_blocks_high || offsetc >= max_blocks_wide)
        continue;

      decode_reconstruct_tx(xd, r, mbmi, plane, plane_bsize, block + i * step,
                            offsetr, offsetc, tx_size - 1, eob_total);
    }
  }
}
#endif  // CONFIG_VAR_TX

#if !CONFIG_VAR_TX || CONFIG_SUPERTX
static int reconstruct_inter_block(MACROBLOCKD *const xd,
#if CONFIG_ANS
                                   struct AnsDecoder *const r,
#else
                                   vp10_reader *r,
#endif
                                   int segment_id, int plane,
                                   int row, int col, TX_SIZE tx_size) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  PLANE_TYPE plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
  int block_idx = (row << 1) + col;
  TX_TYPE tx_type = get_tx_type(plane_type, xd, block_idx, tx_size);
  const scan_order *sc = get_scan(tx_size, tx_type, 1);
  const int eob = vp10_decode_block_tokens(xd,
                                           plane, sc, col, row,
                                           tx_size, tx_type, r,
                                           segment_id);

  inverse_transform_block(xd, plane, tx_type, tx_size,
                          &pd->dst.buf[4 * row * pd->dst.stride + 4 * col],
                          pd->dst.stride, eob);
  return eob;
}
#endif  // !CONFIG_VAR_TX || CONFIG_SUPER_TX

static INLINE TX_SIZE dec_get_uv_tx_size(const MB_MODE_INFO *mbmi,
                                         int n4_wl, int n4_hl) {
  // get minimum log2 num4x4s dimension
  const int x = VPXMIN(n4_wl, n4_hl);
  return VPXMIN(mbmi->tx_size,  x);
}

static INLINE void dec_reset_skip_context(MACROBLOCKD *xd) {
  int i;
  for (i = 0; i < MAX_MB_PLANE; i++) {
    struct macroblockd_plane *const pd = &xd->plane[i];
    memset(pd->above_context, 0, sizeof(ENTROPY_CONTEXT) * pd->n4_w);
    memset(pd->left_context, 0, sizeof(ENTROPY_CONTEXT) * pd->n4_h);
  }
}

static void set_plane_n4(MACROBLOCKD *const xd, int bw, int bh, int bwl,
                         int bhl) {
  int i;
  for (i = 0; i < MAX_MB_PLANE; i++) {
    xd->plane[i].n4_w = (bw << 1) >> xd->plane[i].subsampling_x;
    xd->plane[i].n4_h = (bh << 1) >> xd->plane[i].subsampling_y;
    xd->plane[i].n4_wl = bwl - xd->plane[i].subsampling_x;
    xd->plane[i].n4_hl = bhl - xd->plane[i].subsampling_y;
  }
}

static MB_MODE_INFO *set_offsets(VP10_COMMON *const cm, MACROBLOCKD *const xd,
                                 BLOCK_SIZE bsize, int mi_row, int mi_col,
                                 int bw, int bh, int x_mis, int y_mis,
                                 int bwl, int bhl) {
  const int offset = mi_row * cm->mi_stride + mi_col;
  int x, y;
  const TileInfo *const tile = &xd->tile;

  xd->mi = cm->mi_grid_visible + offset;
  xd->mi[0] = &cm->mi[offset];
  // TODO(slavarnway): Generate sb_type based on bwl and bhl, instead of
  // passing bsize from decode_partition().
  xd->mi[0]->mbmi.sb_type = bsize;
  for (y = 0; y < y_mis; ++y)
    for (x = !y; x < x_mis; ++x) {
      xd->mi[y * cm->mi_stride + x] = xd->mi[0];
    }

  set_plane_n4(xd, bw, bh, bwl, bhl);

  set_skip_context(xd, mi_row, mi_col);


#if CONFIG_VAR_TX
  xd->max_tx_size = max_txsize_lookup[bsize];
#endif

  // Distance of Mb to the various image edges. These are specified to 8th pel
  // as they are always compared to values that are in 1/8th pel units
  set_mi_row_col(xd, tile, mi_row, bh, mi_col, bw, cm->mi_rows, cm->mi_cols);

  vp10_setup_dst_planes(xd->plane, get_frame_new_buffer(cm), mi_row, mi_col);
  return &xd->mi[0]->mbmi;
}

#if CONFIG_SUPERTX
static MB_MODE_INFO *set_offsets_extend(VP10_COMMON *const cm,
                                        MACROBLOCKD *const xd,
                                        const TileInfo *const tile,
                                        BLOCK_SIZE bsize_pred,
                                        int mi_row_pred, int mi_col_pred,
                                        int mi_row_ori, int mi_col_ori) {
  // Used in supertx
  // (mi_row_ori, mi_col_ori): location for mv
  // (mi_row_pred, mi_col_pred, bsize_pred): region to predict
  const int bw = num_8x8_blocks_wide_lookup[bsize_pred];
  const int bh = num_8x8_blocks_high_lookup[bsize_pred];
  const int offset = mi_row_ori * cm->mi_stride + mi_col_ori;
  const int bwl = b_width_log2_lookup[bsize_pred];
  const int bhl = b_height_log2_lookup[bsize_pred];
  xd->mi = cm->mi_grid_visible + offset;
  xd->mi[0] = cm->mi + offset;
  set_mi_row_col(xd, tile, mi_row_pred, bh, mi_col_pred, bw,
                 cm->mi_rows, cm->mi_cols);

  xd->up_available    = (mi_row_ori > tile->mi_row_start);
  xd->left_available  = (mi_col_ori > tile->mi_col_start);

  set_plane_n4(xd, bw, bh, bwl, bhl);

  return &xd->mi[0]->mbmi;
}

static MB_MODE_INFO *set_mb_offsets(VP10_COMMON *const cm,
                                    MACROBLOCKD *const xd,
                                    BLOCK_SIZE bsize,
                                    int mi_row, int mi_col,
                                    int bw, int bh,
                                    int x_mis, int y_mis) {
  const int offset = mi_row * cm->mi_stride + mi_col;
  const TileInfo *const tile = &xd->tile;
  int x, y;

  xd->mi = cm->mi_grid_visible + offset;
  xd->mi[0] = cm->mi + offset;
  xd->mi[0]->mbmi.sb_type = bsize;
  for (y = 0; y < y_mis; ++y)
    for (x = !y; x < x_mis; ++x)
      xd->mi[y * cm->mi_stride + x] = xd->mi[0];

  set_mi_row_col(xd, tile, mi_row, bh, mi_col, bw, cm->mi_rows, cm->mi_cols);
  return &xd->mi[0]->mbmi;
}

static void set_offsets_topblock(VP10_COMMON *const cm, MACROBLOCKD *const xd,
                                 const TileInfo *const tile,
                                 BLOCK_SIZE bsize, int mi_row, int mi_col) {
  const int bw = num_8x8_blocks_wide_lookup[bsize];
  const int bh = num_8x8_blocks_high_lookup[bsize];
  const int offset = mi_row * cm->mi_stride + mi_col;
  const int bwl = b_width_log2_lookup[bsize];
  const int bhl = b_height_log2_lookup[bsize];

  xd->mi = cm->mi_grid_visible + offset;
  xd->mi[0] = cm->mi + offset;

  set_plane_n4(xd, bw, bh, bwl, bhl);

  set_mi_row_col(xd, tile, mi_row, bh, mi_col, bw, cm->mi_rows, cm->mi_cols);

  vp10_setup_dst_planes(xd->plane, get_frame_new_buffer(cm), mi_row, mi_col);
}

static void set_param_topblock(VP10_COMMON *const cm,  MACROBLOCKD *const xd,
                               BLOCK_SIZE bsize, int mi_row, int mi_col,
                               int txfm, int skip) {
  const int bw = num_8x8_blocks_wide_lookup[bsize];
  const int bh = num_8x8_blocks_high_lookup[bsize];
  const int x_mis = VPXMIN(bw, cm->mi_cols - mi_col);
  const int y_mis = VPXMIN(bh, cm->mi_rows - mi_row);
  const int offset = mi_row * cm->mi_stride + mi_col;
  int x, y;

  xd->mi = cm->mi_grid_visible + offset;
  xd->mi[0] = cm->mi + offset;

  for (y = 0; y < y_mis; ++y)
    for (x = 0; x < x_mis; ++x) {
      xd->mi[y * cm->mi_stride + x]->mbmi.skip = skip;
      xd->mi[y * cm->mi_stride + x]->mbmi.tx_type = txfm;
    }
#if CONFIG_VAR_TX
  xd->above_txfm_context = cm->above_txfm_context + mi_col;
  xd->left_txfm_context =
    xd->left_txfm_context_buffer + (mi_row & MAX_MIB_MASK);
  set_txfm_ctx(xd->left_txfm_context, xd->mi[0]->mbmi.tx_size, bh);
  set_txfm_ctx(xd->above_txfm_context, xd->mi[0]->mbmi.tx_size, bw);
#endif
}

static void set_ref(VP10_COMMON *const cm, MACROBLOCKD *const xd,
                    int idx, int mi_row, int mi_col) {
  MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;
  RefBuffer *ref_buffer = &cm->frame_refs[mbmi->ref_frame[idx] - LAST_FRAME];
  xd->block_refs[idx] = ref_buffer;
  if (!vp10_is_valid_scale(&ref_buffer->sf))
    vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                       "Invalid scale factors");
  vp10_setup_pre_planes(xd, idx, ref_buffer->buf, mi_row, mi_col,
                        &ref_buffer->sf);
  xd->corrupted |= ref_buffer->buf->corrupted;
}

static void dec_predict_b_extend(
    VP10Decoder *const pbi, MACROBLOCKD *const xd,
    const TileInfo *const tile, int block,
    int mi_row_ori, int mi_col_ori,
    int mi_row_pred, int mi_col_pred,
    int mi_row_top, int mi_col_top,
    uint8_t * dst_buf[3], int dst_stride[3],
    BLOCK_SIZE bsize_top,
    BLOCK_SIZE bsize_pred,
    int b_sub8x8, int bextend) {
  // Used in supertx
  // (mi_row_ori, mi_col_ori): location for mv
  // (mi_row_pred, mi_col_pred, bsize_pred): region to predict
  // (mi_row_top, mi_col_top, bsize_top): region of the top partition size
  // block: sub location of sub8x8 blocks
  // b_sub8x8: 1: ori is sub8x8; 0: ori is not sub8x8
  // bextend: 1: region to predict is an extension of ori; 0: not
  int r = (mi_row_pred - mi_row_top) * MI_SIZE;
  int c = (mi_col_pred - mi_col_top) * MI_SIZE;
  const int mi_width_top = num_8x8_blocks_wide_lookup[bsize_top];
  const int mi_height_top = num_8x8_blocks_high_lookup[bsize_top];
  MB_MODE_INFO *mbmi;
  VP10_COMMON *const cm = &pbi->common;

  if (mi_row_pred < mi_row_top || mi_col_pred < mi_col_top ||
      mi_row_pred >= mi_row_top + mi_height_top ||
      mi_col_pred >= mi_col_top + mi_width_top ||
      mi_row_pred >= cm->mi_rows || mi_col_pred >= cm->mi_cols)
    return;

  mbmi = set_offsets_extend(cm, xd, tile, bsize_pred,
                            mi_row_pred, mi_col_pred,
                            mi_row_ori, mi_col_ori);
  set_ref(cm, xd, 0, mi_row_pred, mi_col_pred);
  if (has_second_ref(&xd->mi[0]->mbmi))
    set_ref(cm, xd, 1, mi_row_pred, mi_col_pred);

  if (!bextend) {
    mbmi->tx_size = b_width_log2_lookup[bsize_top];
  }

  xd->plane[0].dst.stride = dst_stride[0];
  xd->plane[1].dst.stride = dst_stride[1];
  xd->plane[2].dst.stride = dst_stride[2];
  xd->plane[0].dst.buf = dst_buf[0] +
                         (r >> xd->plane[0].subsampling_y) * dst_stride[0] +
                         (c >> xd->plane[0].subsampling_x);
  xd->plane[1].dst.buf = dst_buf[1] +
                         (r >> xd->plane[1].subsampling_y) * dst_stride[1] +
                         (c >> xd->plane[1].subsampling_x);
  xd->plane[2].dst.buf = dst_buf[2] +
                         (r >> xd->plane[2].subsampling_y) * dst_stride[2] +
                         (c >> xd->plane[2].subsampling_x);

  if (!b_sub8x8)
    vp10_build_inter_predictors_sb_extend(
        xd,
#if CONFIG_EXT_INTER
        mi_row_ori, mi_col_ori,
#endif  // CONFIG_EXT_INTER
        mi_row_pred, mi_col_pred, bsize_pred);
  else
    vp10_build_inter_predictors_sb_sub8x8_extend(
        xd,
#if CONFIG_EXT_INTER
        mi_row_ori, mi_col_ori,
#endif  // CONFIG_EXT_INTER
        mi_row_pred, mi_col_pred, bsize_pred, block);
}

static void dec_extend_dir(VP10Decoder *const pbi, MACROBLOCKD *const xd,
                           const TileInfo *const tile, int block,
                           BLOCK_SIZE bsize, BLOCK_SIZE top_bsize,
                           int mi_row, int mi_col,
                           int mi_row_top, int mi_col_top,
                           uint8_t * dst_buf[3], int dst_stride[3], int dir) {
  // dir: 0-lower, 1-upper, 2-left, 3-right
  //      4-lowerleft, 5-upperleft, 6-lowerright, 7-upperright
  const int mi_width = num_8x8_blocks_wide_lookup[bsize];
  const int mi_height = num_8x8_blocks_high_lookup[bsize];
  int xss = xd->plane[1].subsampling_x;
  int yss = xd->plane[1].subsampling_y;
  int b_sub8x8 = (bsize < BLOCK_8X8) ? 1 : 0;
  BLOCK_SIZE extend_bsize;
  int unit, mi_row_pred, mi_col_pred;

  if (dir == 0 || dir == 1) {
    extend_bsize = (mi_width == 1 || bsize < BLOCK_8X8 || xss < yss) ?
                    BLOCK_8X8 : BLOCK_16X8;
    unit = num_8x8_blocks_wide_lookup[extend_bsize];
    mi_row_pred = mi_row + ((dir == 0) ? mi_height : -1);
    mi_col_pred = mi_col;

    dec_predict_b_extend(pbi, xd, tile, block, mi_row, mi_col,
                         mi_row_pred, mi_col_pred,
                         mi_row_top, mi_col_top,
                         dst_buf, dst_stride,
                         top_bsize, extend_bsize, b_sub8x8, 1);

    if (mi_width > unit) {
      int i;
      assert(!b_sub8x8);
      for (i = 0; i < mi_width/unit - 1; i++) {
        mi_col_pred += unit;
        dec_predict_b_extend(pbi, xd, tile, block, mi_row, mi_col,
                             mi_row_pred, mi_col_pred,
                             mi_row_top, mi_col_top,
                             dst_buf, dst_stride,
                             top_bsize, extend_bsize, b_sub8x8, 1);
      }
    }
  } else if (dir == 2 || dir == 3) {
    extend_bsize = (mi_height == 1 || bsize < BLOCK_8X8 || yss < xss) ?
                    BLOCK_8X8 : BLOCK_8X16;
    unit = num_8x8_blocks_high_lookup[extend_bsize];
    mi_row_pred = mi_row;
    mi_col_pred = mi_col + ((dir == 3) ? mi_width : -1);

    dec_predict_b_extend(pbi, xd, tile, block, mi_row, mi_col,
                         mi_row_pred, mi_col_pred,
                         mi_row_top, mi_col_top,
                         dst_buf, dst_stride,
                         top_bsize, extend_bsize, b_sub8x8, 1);

    if (mi_height > unit) {
      int i;
      for (i = 0; i < mi_height/unit - 1; i++) {
        mi_row_pred += unit;
        dec_predict_b_extend(pbi, xd, tile, block, mi_row, mi_col,
                             mi_row_pred, mi_col_pred,
                             mi_row_top, mi_col_top,
                             dst_buf, dst_stride,
                             top_bsize, extend_bsize, b_sub8x8, 1);
      }
    }
  } else {
    extend_bsize = BLOCK_8X8;
    mi_row_pred = mi_row + ((dir == 4 || dir == 6) ? mi_height : -1);
    mi_col_pred = mi_col + ((dir == 6 || dir == 7) ? mi_width : -1);
    dec_predict_b_extend(pbi, xd, tile, block, mi_row, mi_col,
                         mi_row_pred, mi_col_pred,
                         mi_row_top, mi_col_top,
                         dst_buf, dst_stride,
                         top_bsize, extend_bsize, b_sub8x8, 1);
  }
}

static void dec_extend_all(VP10Decoder *const pbi, MACROBLOCKD *const xd,
                           const TileInfo *const tile, int block,
                           BLOCK_SIZE bsize, BLOCK_SIZE top_bsize,
                           int mi_row, int mi_col,
                           int mi_row_top, int mi_col_top,
                           uint8_t * dst_buf[3], int dst_stride[3]) {
  dec_extend_dir(pbi, xd, tile, block, bsize, top_bsize, mi_row, mi_col,
                 mi_row_top, mi_col_top, dst_buf, dst_stride, 0);
  dec_extend_dir(pbi, xd, tile, block, bsize, top_bsize, mi_row, mi_col,
                 mi_row_top, mi_col_top, dst_buf, dst_stride, 1);
  dec_extend_dir(pbi, xd, tile, block, bsize, top_bsize, mi_row, mi_col,
                 mi_row_top, mi_col_top, dst_buf, dst_stride, 2);
  dec_extend_dir(pbi, xd, tile, block, bsize, top_bsize, mi_row, mi_col,
                 mi_row_top, mi_col_top, dst_buf, dst_stride, 3);
  dec_extend_dir(pbi, xd, tile, block, bsize, top_bsize, mi_row, mi_col,
                 mi_row_top, mi_col_top, dst_buf, dst_stride, 4);
  dec_extend_dir(pbi, xd, tile, block, bsize, top_bsize, mi_row, mi_col,
                 mi_row_top, mi_col_top, dst_buf, dst_stride, 5);
  dec_extend_dir(pbi, xd, tile, block, bsize, top_bsize, mi_row, mi_col,
                 mi_row_top, mi_col_top, dst_buf, dst_stride, 6);
  dec_extend_dir(pbi, xd, tile, block, bsize, top_bsize, mi_row, mi_col,
                 mi_row_top, mi_col_top, dst_buf, dst_stride, 7);
}

static void dec_predict_sb_complex(VP10Decoder *const pbi,
                                   MACROBLOCKD *const xd,
                                   const TileInfo *const tile,
                                   int mi_row, int mi_col,
                                   int mi_row_top, int mi_col_top,
                                   BLOCK_SIZE bsize, BLOCK_SIZE top_bsize,
                                   uint8_t *dst_buf[3], int dst_stride[3]) {
  const VP10_COMMON *const cm = &pbi->common;
  const int hbs = num_8x8_blocks_wide_lookup[bsize] / 2;
  const PARTITION_TYPE partition = get_partition(cm, mi_row, mi_col, bsize);
  const BLOCK_SIZE subsize = get_subsize(bsize, partition);
#if CONFIG_EXT_PARTITION_TYPES
  const BLOCK_SIZE bsize2 = get_subsize(bsize, PARTITION_SPLIT);
#endif
  int i;
  const int mi_offset = mi_row * cm->mi_stride + mi_col;
  uint8_t *dst_buf1[3], *dst_buf2[3], *dst_buf3[3];

  DECLARE_ALIGNED(16, uint8_t,
                  tmp_buf1[MAX_MB_PLANE * MAX_TX_SQUARE * 2]);
  DECLARE_ALIGNED(16, uint8_t,
                  tmp_buf2[MAX_MB_PLANE * MAX_TX_SQUARE * 2]);
  DECLARE_ALIGNED(16, uint8_t,
                  tmp_buf3[MAX_MB_PLANE * MAX_TX_SQUARE * 2]);
  int dst_stride1[3] = {MAX_TX_SIZE, MAX_TX_SIZE, MAX_TX_SIZE};
  int dst_stride2[3] = {MAX_TX_SIZE, MAX_TX_SIZE, MAX_TX_SIZE};
  int dst_stride3[3] = {MAX_TX_SIZE, MAX_TX_SIZE, MAX_TX_SIZE};

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    int len = sizeof(uint16_t);
    dst_buf1[0] = CONVERT_TO_BYTEPTR(tmp_buf1);
    dst_buf1[1] = CONVERT_TO_BYTEPTR(tmp_buf1 + MAX_TX_SQUARE * len);
    dst_buf1[2] = CONVERT_TO_BYTEPTR(tmp_buf1 + 2 * MAX_TX_SQUARE * len);
    dst_buf2[0] = CONVERT_TO_BYTEPTR(tmp_buf2);
    dst_buf2[1] = CONVERT_TO_BYTEPTR(tmp_buf2 + MAX_TX_SQUARE * len);
    dst_buf2[2] = CONVERT_TO_BYTEPTR(tmp_buf2 + 2 * MAX_TX_SQUARE * len);
    dst_buf3[0] = CONVERT_TO_BYTEPTR(tmp_buf3);
    dst_buf3[1] = CONVERT_TO_BYTEPTR(tmp_buf3 + MAX_TX_SQUARE * len);
    dst_buf3[2] = CONVERT_TO_BYTEPTR(tmp_buf3 + 2 * MAX_TX_SQUARE * len);
  } else {
#endif
    dst_buf1[0] = tmp_buf1;
    dst_buf1[1] = tmp_buf1 + MAX_TX_SQUARE;
    dst_buf1[2] = tmp_buf1 + 2 * MAX_TX_SQUARE;
    dst_buf2[0] = tmp_buf2;
    dst_buf2[1] = tmp_buf2 + MAX_TX_SQUARE;
    dst_buf2[2] = tmp_buf2 + 2 * MAX_TX_SQUARE;
    dst_buf3[0] = tmp_buf3;
    dst_buf3[1] = tmp_buf3 + MAX_TX_SQUARE;
    dst_buf3[2] = tmp_buf3 + 2 * MAX_TX_SQUARE;
#if CONFIG_VP9_HIGHBITDEPTH
  }
#endif

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  xd->mi = cm->mi_grid_visible + mi_offset;
  xd->mi[0] = cm->mi + mi_offset;

  for (i = 0; i < MAX_MB_PLANE; i++) {
    xd->plane[i].dst.buf = dst_buf[i];
    xd->plane[i].dst.stride = dst_stride[i];
  }

  switch (partition) {
    case PARTITION_NONE:
      assert(bsize < top_bsize);
      dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col, mi_row, mi_col,
                           mi_row_top, mi_col_top, dst_buf, dst_stride,
                           top_bsize, bsize, 0, 0);
      dec_extend_all(pbi, xd, tile, 0, bsize, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, dst_buf, dst_stride);
      break;
    case PARTITION_HORZ:
      if (bsize == BLOCK_8X8) {
        // For sub8x8, predict in 8x8 unit
        // First half
        dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col, mi_row, mi_col,
                             mi_row_top, mi_col_top, dst_buf, dst_stride,
                             top_bsize, BLOCK_8X8, 1, 0);
        if (bsize < top_bsize)
          dec_extend_all(pbi, xd, tile, 0, subsize, top_bsize, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf, dst_stride);

        // Second half
        dec_predict_b_extend(pbi, xd, tile, 2, mi_row, mi_col, mi_row, mi_col,
                             mi_row_top, mi_col_top, dst_buf1, dst_stride1,
                             top_bsize, BLOCK_8X8, 1, 1);
        if (bsize < top_bsize)
          dec_extend_all(pbi, xd, tile, 2, subsize, top_bsize, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf1, dst_stride1);

        // weighted average to smooth the boundary
        xd->plane[0].dst.buf = dst_buf[0];
        xd->plane[0].dst.stride = dst_stride[0];
        vp10_build_masked_inter_predictor_complex(xd,
                                                  dst_buf[0], dst_stride[0],
                                                  dst_buf1[0], dst_stride1[0],
                                                  mi_row, mi_col,
                                                  mi_row_top, mi_col_top,
                                                  bsize, top_bsize,
                                                  PARTITION_HORZ, 0);
      } else {
        // First half
        dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col, mi_row, mi_col,
                             mi_row_top, mi_col_top, dst_buf, dst_stride,
                             top_bsize, subsize, 0, 0);
        if (bsize < top_bsize)
          dec_extend_all(pbi, xd, tile, 0, subsize, top_bsize, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf, dst_stride);
        else
          dec_extend_dir(pbi, xd, tile, 0, subsize, top_bsize, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf, dst_stride, 0);

        if (mi_row + hbs < cm->mi_rows) {
          // Second half
          dec_predict_b_extend(pbi, xd, tile, 0, mi_row + hbs, mi_col,
                               mi_row + hbs, mi_col,
                               mi_row_top, mi_col_top,
                               dst_buf1, dst_stride1,
                               top_bsize, subsize, 0, 0);
          if (bsize < top_bsize)
            dec_extend_all(pbi, xd, tile, 0, subsize, top_bsize,
                           mi_row + hbs, mi_col,
                           mi_row_top, mi_col_top,
                           dst_buf1, dst_stride1);
          else
            dec_extend_dir(pbi, xd, tile, 0, subsize, top_bsize,
                           mi_row + hbs, mi_col,
                           mi_row_top, mi_col_top,
                           dst_buf1, dst_stride1, 1);

          // weighted average to smooth the boundary
          for (i = 0; i < MAX_MB_PLANE; i++) {
            xd->plane[i].dst.buf = dst_buf[i];
            xd->plane[i].dst.stride = dst_stride[i];
            vp10_build_masked_inter_predictor_complex(
                xd, dst_buf[i], dst_stride[i], dst_buf1[i], dst_stride1[i],
                mi_row, mi_col, mi_row_top, mi_col_top,
                bsize, top_bsize, PARTITION_HORZ, i);
          }
        }
      }
      break;
    case PARTITION_VERT:
      if (bsize == BLOCK_8X8) {
        // First half
        dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col, mi_row, mi_col,
                             mi_row_top, mi_col_top, dst_buf, dst_stride,
                             top_bsize, BLOCK_8X8, 1, 0);
        if (bsize < top_bsize)
          dec_extend_all(pbi, xd, tile, 0, subsize, top_bsize, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf, dst_stride);

        // Second half
        dec_predict_b_extend(pbi, xd, tile, 1, mi_row, mi_col, mi_row, mi_col,
                             mi_row_top, mi_col_top, dst_buf1, dst_stride1,
                             top_bsize, BLOCK_8X8, 1, 1);
        if (bsize < top_bsize)
          dec_extend_all(pbi, xd, tile, 1, subsize, top_bsize, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf1, dst_stride1);

        // Smooth
        xd->plane[0].dst.buf = dst_buf[0];
        xd->plane[0].dst.stride = dst_stride[0];
        vp10_build_masked_inter_predictor_complex(xd,
                                                  dst_buf[0], dst_stride[0],
                                                  dst_buf1[0], dst_stride1[0],
                                                  mi_row, mi_col,
                                                  mi_row_top, mi_col_top,
                                                  bsize, top_bsize,
                                                  PARTITION_VERT, 0);
      } else {
        // First half
        dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col, mi_row, mi_col,
                             mi_row_top, mi_col_top, dst_buf, dst_stride,
                             top_bsize, subsize, 0, 0);
        if (bsize < top_bsize)
          dec_extend_all(pbi, xd, tile, 0, subsize, top_bsize, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf, dst_stride);
        else
          dec_extend_dir(pbi, xd, tile, 0, subsize, top_bsize, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf, dst_stride, 3);

        // Second half
        if (mi_col + hbs < cm->mi_cols) {
          dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col + hbs,
                               mi_row, mi_col + hbs, mi_row_top, mi_col_top,
                               dst_buf1, dst_stride1, top_bsize, subsize, 0, 0);
          if (bsize < top_bsize)
            dec_extend_all(pbi, xd, tile, 0, subsize, top_bsize,
                           mi_row, mi_col + hbs, mi_row_top, mi_col_top,
                           dst_buf1, dst_stride1);
          else
            dec_extend_dir(pbi, xd, tile, 0, subsize, top_bsize,
                           mi_row, mi_col + hbs, mi_row_top, mi_col_top,
                           dst_buf1, dst_stride1, 2);

          // Smooth
          for (i = 0; i < MAX_MB_PLANE; i++) {
            xd->plane[i].dst.buf = dst_buf[i];
            xd->plane[i].dst.stride = dst_stride[i];
            vp10_build_masked_inter_predictor_complex(
                xd, dst_buf[i], dst_stride[i], dst_buf1[i], dst_stride1[i],
                mi_row, mi_col, mi_row_top, mi_col_top,
                bsize, top_bsize, PARTITION_VERT, i);
          }
        }
      }
      break;
    case PARTITION_SPLIT:
      if (bsize == BLOCK_8X8) {
        dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col, mi_row, mi_col,
                             mi_row_top, mi_col_top, dst_buf, dst_stride,
                             top_bsize, BLOCK_8X8, 1, 0);
        dec_predict_b_extend(pbi, xd, tile, 1, mi_row, mi_col, mi_row, mi_col,
                             mi_row_top, mi_col_top, dst_buf1, dst_stride1,
                             top_bsize, BLOCK_8X8, 1, 1);
        dec_predict_b_extend(pbi, xd, tile, 2, mi_row, mi_col, mi_row, mi_col,
                             mi_row_top, mi_col_top, dst_buf2, dst_stride2,
                             top_bsize, BLOCK_8X8, 1, 1);
        dec_predict_b_extend(pbi, xd, tile, 3, mi_row, mi_col, mi_row, mi_col,
                             mi_row_top, mi_col_top, dst_buf3, dst_stride3,
                             top_bsize, BLOCK_8X8, 1, 1);
        if (bsize < top_bsize) {
          dec_extend_all(pbi, xd, tile, 0, subsize, top_bsize, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf, dst_stride);
          dec_extend_all(pbi, xd, tile, 1, subsize, top_bsize, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf1, dst_stride1);
          dec_extend_all(pbi, xd, tile, 2, subsize, top_bsize, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf2, dst_stride2);
          dec_extend_all(pbi, xd, tile, 3, subsize, top_bsize, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf3, dst_stride3);
        }
      } else {
        dec_predict_sb_complex(pbi, xd, tile, mi_row, mi_col,
                               mi_row_top, mi_col_top, subsize, top_bsize,
                               dst_buf, dst_stride);
        if (mi_row < cm->mi_rows && mi_col + hbs < cm->mi_cols)
          dec_predict_sb_complex(pbi, xd, tile, mi_row, mi_col + hbs,
                                 mi_row_top, mi_col_top, subsize, top_bsize,
                                 dst_buf1, dst_stride1);
        if (mi_row + hbs < cm->mi_rows && mi_col < cm->mi_cols)
          dec_predict_sb_complex(pbi, xd, tile, mi_row + hbs, mi_col,
                                 mi_row_top, mi_col_top, subsize, top_bsize,
                                 dst_buf2, dst_stride2);
        if (mi_row + hbs < cm->mi_rows && mi_col + hbs < cm->mi_cols)
          dec_predict_sb_complex(pbi, xd, tile, mi_row + hbs, mi_col + hbs,
                                 mi_row_top, mi_col_top, subsize, top_bsize,
                                 dst_buf3, dst_stride3);
      }
        for (i = 0; i < MAX_MB_PLANE; i++) {
          if (bsize == BLOCK_8X8 && i != 0)
            continue;  // Skip <4x4 chroma smoothing
          if (mi_row < cm->mi_rows && mi_col + hbs < cm->mi_cols) {
            vp10_build_masked_inter_predictor_complex(xd,
                                                      dst_buf[i], dst_stride[i],
                                                      dst_buf1[i],
                                                      dst_stride1[i],
                                                      mi_row, mi_col,
                                                      mi_row_top, mi_col_top,
                                                      bsize, top_bsize,
                                                      PARTITION_VERT, i);
            if (mi_row + hbs < cm->mi_rows) {
              vp10_build_masked_inter_predictor_complex(xd,
                                                        dst_buf2[i],
                                                        dst_stride2[i],
                                                        dst_buf3[i],
                                                        dst_stride3[i],
                                                        mi_row, mi_col,
                                                        mi_row_top, mi_col_top,
                                                        bsize, top_bsize,
                                                        PARTITION_VERT, i);
              vp10_build_masked_inter_predictor_complex(xd,
                                                        dst_buf[i],
                                                        dst_stride[i],
                                                        dst_buf2[i],
                                                        dst_stride2[i],
                                                        mi_row, mi_col,
                                                        mi_row_top, mi_col_top,
                                                        bsize, top_bsize,
                                                        PARTITION_HORZ, i);
            }
          } else if (mi_row + hbs < cm->mi_rows && mi_col < cm->mi_cols) {
            vp10_build_masked_inter_predictor_complex(xd,
                                                      dst_buf[i],
                                                      dst_stride[i],
                                                      dst_buf2[i],
                                                      dst_stride2[i],
                                                      mi_row, mi_col,
                                                      mi_row_top, mi_col_top,
                                                      bsize, top_bsize,
                                                      PARTITION_HORZ, i);
          }
        }
      break;
#if CONFIG_EXT_PARTITION_TYPES
    case PARTITION_HORZ_A:
      dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col, mi_row, mi_col,
                           mi_row_top, mi_col_top, dst_buf, dst_stride,
                           top_bsize, bsize2, 0, 0);
      dec_extend_all(pbi, xd, tile, 0, bsize2, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, dst_buf, dst_stride);

      dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col + hbs, mi_row,
                           mi_col + hbs, mi_row_top, mi_col_top,
                           dst_buf1, dst_stride1, top_bsize, bsize2, 0, 0);
      dec_extend_all(pbi, xd, tile, 0, bsize2, top_bsize, mi_row, mi_col + hbs,
                     mi_row_top, mi_col_top, dst_buf1, dst_stride1);

      dec_predict_b_extend(pbi, xd, tile, 0, mi_row + hbs, mi_col,
                           mi_row + hbs, mi_col, mi_row_top, mi_col_top,
                           dst_buf2, dst_stride2, top_bsize, subsize, 0, 0);
      if (bsize < top_bsize)
        dec_extend_all(pbi, xd, tile, 0, subsize, top_bsize,
                       mi_row + hbs, mi_col, mi_row_top, mi_col_top,
                       dst_buf2, dst_stride2);
      else
        dec_extend_dir(pbi, xd, tile, 0, subsize, top_bsize,
                       mi_row + hbs, mi_col, mi_row_top, mi_col_top,
                       dst_buf2, dst_stride2, 1);

      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = dst_buf[i];
        xd->plane[i].dst.stride = dst_stride[i];
        vp10_build_masked_inter_predictor_complex(xd,
                                                  dst_buf[i], dst_stride[i],
                                                  dst_buf1[i], dst_stride1[i],
                                                  mi_row, mi_col,
                                                  mi_row_top, mi_col_top,
                                                  bsize, top_bsize,
                                                  PARTITION_VERT, i);
      }
      for (i = 0; i < MAX_MB_PLANE; i++) {
        vp10_build_masked_inter_predictor_complex(xd,
                                                  dst_buf[i], dst_stride[i],
                                                  dst_buf2[i], dst_stride2[i],
                                                  mi_row, mi_col,
                                                  mi_row_top, mi_col_top,
                                                  bsize, top_bsize,
                                                  PARTITION_HORZ, i);
      }
      break;
    case PARTITION_VERT_A:

      dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col, mi_row, mi_col,
                           mi_row_top, mi_col_top, dst_buf, dst_stride,
                           top_bsize, bsize2, 0, 0);
      dec_extend_all(pbi, xd, tile, 0, bsize2, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, dst_buf, dst_stride);

      dec_predict_b_extend(pbi, xd, tile, 0, mi_row + hbs, mi_col,
                           mi_row + hbs, mi_col, mi_row_top, mi_col_top,
                           dst_buf1, dst_stride1, top_bsize, bsize2, 0, 0);
      dec_extend_all(pbi, xd, tile, 0, bsize2, top_bsize, mi_row + hbs, mi_col,
                     mi_row_top, mi_col_top, dst_buf1, dst_stride1);

      dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col + hbs,
                           mi_row, mi_col + hbs, mi_row_top, mi_col_top,
                           dst_buf2, dst_stride2,
                           top_bsize, subsize, 0, 0);
      if (bsize < top_bsize)
        dec_extend_all(pbi, xd, tile, 0, subsize, top_bsize,
                       mi_row, mi_col + hbs, mi_row_top, mi_col_top,
                       dst_buf2, dst_stride2);
      else
        dec_extend_dir(pbi, xd, tile, 0, subsize, top_bsize,
                       mi_row, mi_col + hbs, mi_row_top, mi_col_top,
                       dst_buf2, dst_stride2, 2);

      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = dst_buf[i];
        xd->plane[i].dst.stride = dst_stride[i];
        vp10_build_masked_inter_predictor_complex(xd,
                                                  dst_buf[i], dst_stride[i],
                                                  dst_buf1[i], dst_stride1[i],
                                                  mi_row, mi_col,
                                                  mi_row_top, mi_col_top,
                                                  bsize, top_bsize,
                                                  PARTITION_HORZ, i);
      }
      for (i = 0; i < MAX_MB_PLANE; i++) {
        vp10_build_masked_inter_predictor_complex(xd,
                                                  dst_buf[i], dst_stride[i],
                                                  dst_buf2[i], dst_stride2[i],
                                                  mi_row, mi_col,
                                                  mi_row_top, mi_col_top,
                                                  bsize, top_bsize,
                                                  PARTITION_VERT, i);
      }
      break;
    case PARTITION_HORZ_B:
      dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col, mi_row, mi_col,
                           mi_row_top, mi_col_top, dst_buf, dst_stride,
                           top_bsize, subsize, 0, 0);
      if (bsize < top_bsize)
        dec_extend_all(pbi, xd, tile, 0, subsize, top_bsize, mi_row, mi_col,
                       mi_row_top, mi_col_top, dst_buf, dst_stride);
      else
        dec_extend_dir(pbi, xd, tile, 0, subsize, top_bsize, mi_row, mi_col,
                       mi_row_top, mi_col_top, dst_buf, dst_stride, 0);

      dec_predict_b_extend(pbi, xd, tile, 0, mi_row + hbs, mi_col, mi_row + hbs,
                           mi_col, mi_row_top, mi_col_top,
                           dst_buf1, dst_stride1, top_bsize, bsize2, 0, 0);
      dec_extend_all(pbi, xd, tile, 0, bsize2, top_bsize, mi_row + hbs, mi_col,
                     mi_row_top, mi_col_top, dst_buf1, dst_stride1);

      dec_predict_b_extend(pbi, xd, tile, 0, mi_row + hbs, mi_col + hbs,
                           mi_row + hbs, mi_col + hbs, mi_row_top, mi_col_top,
                           dst_buf2, dst_stride2, top_bsize, bsize2, 0, 0);
      dec_extend_all(pbi, xd, tile, 0, bsize2, top_bsize,
                     mi_row + hbs, mi_col + hbs,
                     mi_row_top, mi_col_top, dst_buf2, dst_stride2);

      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = dst_buf1[i];
        xd->plane[i].dst.stride = dst_stride1[i];
        vp10_build_masked_inter_predictor_complex(xd,
                                                  dst_buf1[i], dst_stride1[i],
                                                  dst_buf2[i], dst_stride2[i],
                                                  mi_row, mi_col,
                                                  mi_row_top, mi_col_top,
                                                  bsize, top_bsize,
                                                  PARTITION_VERT, i);
      }
      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = dst_buf[i];
        xd->plane[i].dst.stride = dst_stride[i];
        vp10_build_masked_inter_predictor_complex(xd,
                                                  dst_buf[i], dst_stride[i],
                                                  dst_buf1[i], dst_stride1[i],
                                                  mi_row, mi_col,
                                                  mi_row_top, mi_col_top,
                                                  bsize, top_bsize,
                                                  PARTITION_HORZ, i);
      }
      break;
    case PARTITION_VERT_B:
      dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col, mi_row, mi_col,
                           mi_row_top, mi_col_top, dst_buf, dst_stride,
                           top_bsize, subsize, 0, 0);
      if (bsize < top_bsize)
        dec_extend_all(pbi, xd, tile, 0, subsize, top_bsize, mi_row, mi_col,
                       mi_row_top, mi_col_top, dst_buf, dst_stride);
      else
        dec_extend_dir(pbi, xd, tile, 0, subsize, top_bsize, mi_row, mi_col,
                       mi_row_top, mi_col_top, dst_buf, dst_stride, 3);

      dec_predict_b_extend(pbi, xd, tile, 0, mi_row, mi_col + hbs,
                           mi_row, mi_col + hbs, mi_row_top, mi_col_top,
                           dst_buf1, dst_stride1, top_bsize, bsize2, 0, 0);
      dec_extend_all(pbi, xd, tile, 0, bsize2, top_bsize, mi_row, mi_col + hbs,
                     mi_row_top, mi_col_top, dst_buf1, dst_stride1);

      dec_predict_b_extend(pbi, xd, tile, 0, mi_row + hbs, mi_col + hbs,
                           mi_row + hbs, mi_col + hbs, mi_row_top, mi_col_top,
                           dst_buf2, dst_stride2, top_bsize, bsize2, 0, 0);
      dec_extend_all(pbi, xd, tile, 0, bsize2, top_bsize,
                     mi_row + hbs, mi_col + hbs,
                     mi_row_top, mi_col_top, dst_buf2, dst_stride2);

      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = dst_buf1[i];
        xd->plane[i].dst.stride = dst_stride1[i];
        vp10_build_masked_inter_predictor_complex(xd,
                                                  dst_buf1[i], dst_stride1[i],
                                                  dst_buf2[i], dst_stride2[i],
                                                  mi_row, mi_col,
                                                  mi_row_top, mi_col_top,
                                                  bsize, top_bsize,
                                                  PARTITION_HORZ, i);
      }
      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = dst_buf[i];
        xd->plane[i].dst.stride = dst_stride[i];
        vp10_build_masked_inter_predictor_complex(xd,
                                                  dst_buf[i], dst_stride[i],
                                                  dst_buf1[i], dst_stride1[i],
                                                  mi_row, mi_col,
                                                  mi_row_top, mi_col_top,
                                                  bsize, top_bsize,
                                                  PARTITION_VERT, i);
      }
      break;
#endif  // CONFIG_EXT_PARTITION_TYPES
    default:
      assert(0);
  }
}

static void set_segment_id_supertx(const VP10_COMMON *const cm,
                                   const int mi_row, const int mi_col,
                                   const BLOCK_SIZE bsize) {
  const struct segmentation *seg = &cm->seg;
  const int miw =
      VPXMIN(num_8x8_blocks_wide_lookup[bsize], cm->mi_cols - mi_col);
  const int mih =
      VPXMIN(num_8x8_blocks_high_lookup[bsize], cm->mi_rows - mi_row);
  const int mi_offset = mi_row * cm->mi_stride + mi_col;
  MODE_INFO **const mip = cm->mi_grid_visible + mi_offset;
  int r, c;
  int seg_id_supertx = MAX_SEGMENTS;

  if (!seg->enabled) {
    seg_id_supertx = 0;
  } else {
    // Find the minimum segment_id
    for (r = 0 ; r < mih ; r++)
      for (c = 0 ; c < miw ; c++)
        seg_id_supertx = VPXMIN(mip[r * cm->mi_stride + c]->mbmi.segment_id,
                                seg_id_supertx);
    assert(0 <= seg_id_supertx && seg_id_supertx < MAX_SEGMENTS);
  }

  // Assign the the segment_id back to segment_id_supertx
  for (r = 0 ; r < mih ; r++)
    for (c = 0 ; c < miw ; c++)
      mip[r * cm->mi_stride + c]->mbmi.segment_id_supertx = seg_id_supertx;
}
#endif  // CONFIG_SUPERTX

static void decode_block(VP10Decoder *const pbi, MACROBLOCKD *const xd,
#if CONFIG_SUPERTX
                         int supertx_enabled,
#endif  // CONFIG_SUPERTX
                         int mi_row, int mi_col,
                         vp10_reader *r,
#if CONFIG_EXT_PARTITION_TYPES
                         PARTITION_TYPE partition,
#endif  // CONFIG_EXT_PARTITION_TYPES
                         BLOCK_SIZE bsize,
                         int bwl, int bhl) {
  VP10_COMMON *const cm = &pbi->common;
  const int less8x8 = bsize < BLOCK_8X8;
  const int bw = 1 << (bwl - 1);
  const int bh = 1 << (bhl - 1);
  const int x_mis = VPXMIN(bw, cm->mi_cols - mi_col);
  const int y_mis = VPXMIN(bh, cm->mi_rows - mi_row);

#if CONFIG_SUPERTX
  MB_MODE_INFO *mbmi;
  if (supertx_enabled) {
    mbmi = set_mb_offsets(cm, xd, bsize, mi_row, mi_col,
                          bw, bh, x_mis, y_mis);
  } else {
    mbmi = set_offsets(cm, xd, bsize, mi_row, mi_col,
                       bw, bh, x_mis, y_mis, bwl, bhl);
  }
#if CONFIG_EXT_PARTITION_TYPES
  xd->mi[0]->mbmi.partition = partition;
#endif
  vp10_read_mode_info(pbi, xd, supertx_enabled,
                      mi_row, mi_col, r, x_mis, y_mis);
#else
  MB_MODE_INFO *mbmi = set_offsets(cm, xd, bsize, mi_row, mi_col,
                                   bw, bh, x_mis, y_mis, bwl, bhl);
#if CONFIG_EXT_PARTITION_TYPES
  xd->mi[0]->mbmi.partition = partition;
#endif
  vp10_read_mode_info(pbi, xd, mi_row, mi_col, r, x_mis, y_mis);
#endif  // CONFIG_SUPERTX

  if (bsize >= BLOCK_8X8 && (cm->subsampling_x || cm->subsampling_y)) {
    const BLOCK_SIZE uv_subsize =
        ss_size_lookup[bsize][cm->subsampling_x][cm->subsampling_y];
    if (uv_subsize == BLOCK_INVALID)
      vpx_internal_error(xd->error_info,
                         VPX_CODEC_CORRUPT_FRAME, "Invalid block size.");
  }

#if CONFIG_SUPERTX
  mbmi->segment_id_supertx = MAX_SEGMENTS;

  if (supertx_enabled) {
    xd->corrupted |= vp10_reader_has_error(r);
    return;
  }
#endif  // CONFIG_SUPERTX

  if (mbmi->skip) {
    dec_reset_skip_context(xd);
  }
  if (!is_inter_block(mbmi)) {
    int plane;
    for (plane = 0; plane <= 1; ++plane) {
      if (mbmi->palette_mode_info.palette_size[plane])
        vp10_decode_palette_tokens(xd, plane, r);
    }
    for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
      const struct macroblockd_plane *const pd = &xd->plane[plane];
      const TX_SIZE tx_size =
          plane ? dec_get_uv_tx_size(mbmi, pd->n4_wl, pd->n4_hl)
          : mbmi->tx_size;
      const int num_4x4_w = pd->n4_w;
      const int num_4x4_h = pd->n4_h;
      const int step = (1 << tx_size);
      int row, col;
      const int max_blocks_wide = num_4x4_w +
          (xd->mb_to_right_edge >= 0 ?
           0 : xd->mb_to_right_edge >> (5 + pd->subsampling_x));
      const int max_blocks_high = num_4x4_h +
          (xd->mb_to_bottom_edge >= 0 ?
           0 : xd->mb_to_bottom_edge >> (5 + pd->subsampling_y));

      for (row = 0; row < max_blocks_high; row += step)
        for (col = 0; col < max_blocks_wide; col += step)
          predict_and_reconstruct_intra_block(xd,
                                              r,
                                              mbmi, plane,
                                              row, col, tx_size);
    }
  } else {
    // Prediction
    vp10_build_inter_predictors_sb(xd, mi_row, mi_col,
                                   VPXMAX(bsize, BLOCK_8X8));
#if CONFIG_OBMC
    if (mbmi->motion_variation == OBMC_CAUSAL) {
#if CONFIG_VP9_HIGHBITDEPTH
      DECLARE_ALIGNED(16, uint8_t,
                      tmp_buf1[2 * MAX_MB_PLANE * MAX_SB_SQUARE]);
      DECLARE_ALIGNED(16, uint8_t,
                      tmp_buf2[2 * MAX_MB_PLANE * MAX_SB_SQUARE]);
#else
      DECLARE_ALIGNED(16, uint8_t,
                      tmp_buf1[MAX_MB_PLANE * MAX_SB_SQUARE]);
      DECLARE_ALIGNED(16, uint8_t,
                      tmp_buf2[MAX_MB_PLANE * MAX_SB_SQUARE]);
#endif  // CONFIG_VP9_HIGHBITDEPTH
      uint8_t *dst_buf1[MAX_MB_PLANE], *dst_buf2[MAX_MB_PLANE];
      int dst_stride1[MAX_MB_PLANE] = {MAX_SB_SIZE, MAX_SB_SIZE, MAX_SB_SIZE};
      int dst_stride2[MAX_MB_PLANE] = {MAX_SB_SIZE, MAX_SB_SIZE, MAX_SB_SIZE};

      assert(mbmi->sb_type >= BLOCK_8X8);
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        int len = sizeof(uint16_t);
        dst_buf1[0] = CONVERT_TO_BYTEPTR(tmp_buf1);
        dst_buf1[1] = CONVERT_TO_BYTEPTR(tmp_buf1 + MAX_SB_SQUARE * len);
        dst_buf1[2] = CONVERT_TO_BYTEPTR(tmp_buf1 + MAX_SB_SQUARE * 2 * len);
        dst_buf2[0] = CONVERT_TO_BYTEPTR(tmp_buf2);
        dst_buf2[1] = CONVERT_TO_BYTEPTR(tmp_buf2 + MAX_SB_SQUARE * len);
        dst_buf2[2] = CONVERT_TO_BYTEPTR(tmp_buf2 + MAX_SB_SQUARE * 2 * len);
      } else {
#endif  // CONFIG_VP9_HIGHBITDEPTH
        dst_buf1[0] = tmp_buf1;
        dst_buf1[1] = tmp_buf1 + MAX_SB_SQUARE;
        dst_buf1[2] = tmp_buf1 + MAX_SB_SQUARE * 2;
        dst_buf2[0] = tmp_buf2;
        dst_buf2[1] = tmp_buf2 + MAX_SB_SQUARE;
        dst_buf2[2] = tmp_buf2 + MAX_SB_SQUARE * 2;
#if CONFIG_VP9_HIGHBITDEPTH
      }
#endif  // CONFIG_VP9_HIGHBITDEPTH
      vp10_build_prediction_by_above_preds(cm, xd, mi_row, mi_col,
                                           dst_buf1, dst_stride1);
      vp10_build_prediction_by_left_preds(cm, xd, mi_row, mi_col,
                                          dst_buf2, dst_stride2);
      vp10_setup_dst_planes(xd->plane, get_frame_new_buffer(cm),
                            mi_row, mi_col);
      vp10_build_obmc_inter_prediction(cm, xd, mi_row, mi_col, 0, NULL, NULL,
                                       dst_buf1, dst_stride1,
                                       dst_buf2, dst_stride2);
    }
#endif  // CONFIG_OBMC

    // Reconstruction
    if (!mbmi->skip) {
      int eobtotal = 0;
      int plane;

      for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
        const struct macroblockd_plane *const pd = &xd->plane[plane];
        const int num_4x4_w = pd->n4_w;
        const int num_4x4_h = pd->n4_h;
        int row, col;
#if CONFIG_VAR_TX
        // TODO(jingning): This can be simplified for decoder performance.
        const BLOCK_SIZE plane_bsize =
            get_plane_block_size(VPXMAX(bsize, BLOCK_8X8), pd);
        const TX_SIZE max_tx_size = max_txsize_lookup[plane_bsize];
        const BLOCK_SIZE txb_size = txsize_to_bsize[max_tx_size];
        int bw = num_4x4_blocks_wide_lookup[txb_size];
        int block = 0;
        const int step = 1 << (max_tx_size << 1);

        for (row = 0; row < num_4x4_h; row += bw) {
          for (col = 0; col < num_4x4_w; col += bw) {
            decode_reconstruct_tx(xd, r, mbmi, plane, plane_bsize,
                                  block, row, col, max_tx_size, &eobtotal);
            block += step;
          }
        }
#else
        const TX_SIZE tx_size =
            plane ? dec_get_uv_tx_size(mbmi, pd->n4_wl, pd->n4_hl)
            : mbmi->tx_size;
        const int step = (1 << tx_size);
        const int max_blocks_wide = num_4x4_w +
            (xd->mb_to_right_edge >= 0 ?
             0 : xd->mb_to_right_edge >> (5 + pd->subsampling_x));
        const int max_blocks_high = num_4x4_h +
            (xd->mb_to_bottom_edge >= 0 ?
             0 : xd->mb_to_bottom_edge >> (5 + pd->subsampling_y));

        for (row = 0; row < max_blocks_high; row += step)
          for (col = 0; col < max_blocks_wide; col += step)
            eobtotal += reconstruct_inter_block(xd,
                                                r,
                                                mbmi->segment_id,
                                                plane, row, col,
                                                tx_size);
#endif
      }

      if (!less8x8 && eobtotal == 0)
        mbmi->has_no_coeffs = 1;  // skip loopfilter
    }
  }

  xd->corrupted |= vp10_reader_has_error(r);
}

static INLINE int dec_partition_plane_context(const MACROBLOCKD *xd,
                                              int mi_row, int mi_col,
                                              int bsl) {
  const PARTITION_CONTEXT *above_ctx = xd->above_seg_context + mi_col;
  const PARTITION_CONTEXT *left_ctx =
    xd->left_seg_context + (mi_row & MAX_MIB_MASK);
  int above = (*above_ctx >> bsl) & 1 , left = (*left_ctx >> bsl) & 1;

//  assert(bsl >= 0);

  return (left * 2 + above) + bsl * PARTITION_PLOFFSET;
}

#if !CONFIG_EXT_PARTITION_TYPES
static INLINE void dec_update_partition_context(MACROBLOCKD *xd,
                                                int mi_row, int mi_col,
                                                BLOCK_SIZE subsize,
                                                int bw) {
  PARTITION_CONTEXT *const above_ctx = xd->above_seg_context + mi_col;
  PARTITION_CONTEXT *const left_ctx =
    xd->left_seg_context + (mi_row & MAX_MIB_MASK);

  // update the partition context at the end notes. set partition bits
  // of block sizes larger than the current one to be one, and partition
  // bits of smaller block sizes to be zero.
  memset(above_ctx, partition_context_lookup[subsize].above, bw);
  memset(left_ctx, partition_context_lookup[subsize].left, bw);
}
#endif  // !CONFIG_EXT_PARTITION_TYPES

static PARTITION_TYPE read_partition(VP10_COMMON *cm, MACROBLOCKD *xd,
                                     int mi_row, int mi_col, vp10_reader *r,
                                     int has_rows, int has_cols,
#if CONFIG_EXT_PARTITION_TYPES
                                     BLOCK_SIZE bsize,
#endif
                                     int bsl) {
  const int ctx = dec_partition_plane_context(xd, mi_row, mi_col, bsl);
  const vpx_prob *const probs = cm->fc->partition_prob[ctx];
  FRAME_COUNTS *counts = xd->counts;
  PARTITION_TYPE p;

  if (has_rows && has_cols)
#if CONFIG_EXT_PARTITION_TYPES
    if (bsize <= BLOCK_8X8)
      p = (PARTITION_TYPE)vp10_read_tree(r, vp10_partition_tree, probs);
    else
      p = (PARTITION_TYPE)vp10_read_tree(r, vp10_ext_partition_tree, probs);
#else
    p = (PARTITION_TYPE)vp10_read_tree(r, vp10_partition_tree, probs);
#endif  // CONFIG_EXT_PARTITION_TYPES
  else if (!has_rows && has_cols)
    p = vp10_read(r, probs[1]) ? PARTITION_SPLIT : PARTITION_HORZ;
  else if (has_rows && !has_cols)
    p = vp10_read(r, probs[2]) ? PARTITION_SPLIT : PARTITION_VERT;
  else
    p = PARTITION_SPLIT;

  if (counts)
    ++counts->partition[ctx][p];

  return p;
}

#if CONFIG_SUPERTX
static int read_skip(VP10_COMMON *cm, const MACROBLOCKD *xd,
                     int segment_id, vp10_reader *r) {
  if (segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP)) {
    return 1;
  } else {
    const int ctx = vp10_get_skip_context(xd);
    const int skip = vp10_read(r, cm->fc->skip_probs[ctx]);
    FRAME_COUNTS *counts = xd->counts;
    if (counts)
      ++counts->skip[ctx][skip];
    return skip;
  }
}
#endif  // CONFIG_SUPERTX

// TODO(slavarnway): eliminate bsize and subsize in future commits
static void decode_partition(VP10Decoder *const pbi, MACROBLOCKD *const xd,
#if CONFIG_SUPERTX
                             int supertx_enabled,
#endif
                             int mi_row, int mi_col,
                             vp10_reader* r,
                             BLOCK_SIZE bsize, int n4x4_l2) {
  VP10_COMMON *const cm = &pbi->common;
  const int n8x8_l2 = n4x4_l2 - 1;
  const int num_8x8_wh = 1 << n8x8_l2;
  const int hbs = num_8x8_wh >> 1;
  PARTITION_TYPE partition;
  BLOCK_SIZE subsize;
#if CONFIG_EXT_PARTITION_TYPES
  BLOCK_SIZE bsize2 = get_subsize(bsize, PARTITION_SPLIT);
#endif
  const int has_rows = (mi_row + hbs) < cm->mi_rows;
  const int has_cols = (mi_col + hbs) < cm->mi_cols;
#if CONFIG_SUPERTX
  const int read_token = !supertx_enabled;
  int skip = 0;
  TX_SIZE supertx_size = b_width_log2_lookup[bsize];
  const TileInfo *const tile = &xd->tile;
  int txfm = DCT_DCT;
#endif  // CONFIG_SUPERTX

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  partition = read_partition(cm, xd, mi_row, mi_col, r, has_rows, has_cols,
#if CONFIG_EXT_PARTITION_TYPES
                             bsize,
#endif
                             n8x8_l2);
  subsize = subsize_lookup[partition][bsize];  // get_subsize(bsize, partition);
#if CONFIG_SUPERTX
  if (!frame_is_intra_only(cm) &&
      partition != PARTITION_NONE &&
      bsize <= MAX_SUPERTX_BLOCK_SIZE &&
      !supertx_enabled &&
      !xd->lossless[0]) {
    const int supertx_context =
        partition_supertx_context_lookup[partition];
    supertx_enabled = vp10_read(
        r, cm->fc->supertx_prob[supertx_context][supertx_size]);
    if (xd->counts)
      xd->counts->supertx[supertx_context][supertx_size][supertx_enabled]++;
#if CONFIG_VAR_TX
    if (supertx_enabled)
      xd->supertx_size = supertx_size;
#endif
  }
#endif  // CONFIG_SUPERTX
  if (!hbs) {
    // calculate bmode block dimensions (log 2)
    xd->bmode_blocks_wl = 1 >> !!(partition & PARTITION_VERT);
    xd->bmode_blocks_hl = 1 >> !!(partition & PARTITION_HORZ);
    decode_block(pbi, xd,
#if CONFIG_SUPERTX
                 supertx_enabled,
#endif  // CONFIG_SUPERTX
                 mi_row, mi_col, r,
#if CONFIG_EXT_PARTITION_TYPES
                 partition,
#endif  // CONFIG_EXT_PARTITION_TYPES
                 subsize, 1, 1);
  } else {
    switch (partition) {
      case PARTITION_NONE:
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif  // CONFIG_SUPERTX
                     mi_row, mi_col, r,
#if CONFIG_EXT_PARTITION_TYPES
                 partition,
#endif  // CONFIG_EXT_PARTITION_TYPES
                     subsize, n4x4_l2, n4x4_l2);
        break;
      case PARTITION_HORZ:
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif  // CONFIG_SUPERTX
                     mi_row, mi_col, r,
#if CONFIG_EXT_PARTITION_TYPES
                 partition,
#endif  // CONFIG_EXT_PARTITION_TYPES
                     subsize, n4x4_l2, n8x8_l2);
        if (has_rows)
          decode_block(pbi, xd,
#if CONFIG_SUPERTX
                       supertx_enabled,
#endif  // CONFIG_SUPERTX
                       mi_row + hbs, mi_col, r,
#if CONFIG_EXT_PARTITION_TYPES
                 partition,
#endif  // CONFIG_EXT_PARTITION_TYPES
                       subsize, n4x4_l2, n8x8_l2);
        break;
      case PARTITION_VERT:
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif  // CONFIG_SUPERTX
                     mi_row, mi_col, r,
#if CONFIG_EXT_PARTITION_TYPES
                 partition,
#endif  // CONFIG_EXT_PARTITION_TYPES
                     subsize, n8x8_l2, n4x4_l2);
        if (has_cols)
          decode_block(pbi, xd,
#if CONFIG_SUPERTX
                       supertx_enabled,
#endif  // CONFIG_SUPERTX
                       mi_row, mi_col + hbs, r,
#if CONFIG_EXT_PARTITION_TYPES
                 partition,
#endif  // CONFIG_EXT_PARTITION_TYPES
                       subsize, n8x8_l2, n4x4_l2);
        break;
      case PARTITION_SPLIT:
        decode_partition(pbi, xd,
#if CONFIG_SUPERTX
                         supertx_enabled,
#endif  // CONFIG_SUPERTX
                         mi_row, mi_col, r,
                         subsize, n8x8_l2);
        decode_partition(pbi, xd,
#if CONFIG_SUPERTX
                         supertx_enabled,
#endif  // CONFIG_SUPERTX
                         mi_row, mi_col + hbs, r,
                         subsize, n8x8_l2);
        decode_partition(pbi, xd,
#if CONFIG_SUPERTX
                         supertx_enabled,
#endif  // CONFIG_SUPERTX
                         mi_row + hbs, mi_col, r,
                         subsize, n8x8_l2);
        decode_partition(pbi, xd,
#if CONFIG_SUPERTX
                         supertx_enabled,
#endif  // CONFIG_SUPERTX
                         mi_row + hbs, mi_col + hbs, r,
                         subsize, n8x8_l2);
        break;
#if CONFIG_EXT_PARTITION_TYPES
      case PARTITION_HORZ_A:
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row,       mi_col,       r,
                     partition, bsize2, n8x8_l2, n8x8_l2);
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row,       mi_col + hbs, r,
                     partition, bsize2, n8x8_l2, n8x8_l2);
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row + hbs, mi_col, r,
                     partition, subsize, n4x4_l2, n8x8_l2);
        break;
      case PARTITION_HORZ_B:
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row, mi_col, r,
                     partition, subsize, n4x4_l2, n8x8_l2);
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row + hbs, mi_col,       r,
                     partition, bsize2, n8x8_l2, n8x8_l2);
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row + hbs, mi_col + hbs, r,
                     partition, bsize2, n8x8_l2, n8x8_l2);
        break;
      case PARTITION_VERT_A:
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row,       mi_col,       r,
                     partition, bsize2, n8x8_l2, n8x8_l2);
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row + hbs, mi_col,       r,
                     partition, bsize2, n8x8_l2, n8x8_l2);
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row, mi_col + hbs, r,
                     partition, subsize, n8x8_l2, n4x4_l2);
        break;
      case PARTITION_VERT_B:
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row, mi_col, r,
                     partition, subsize, n8x8_l2, n4x4_l2);
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row,       mi_col + hbs, r,
                     partition, bsize2, n8x8_l2, n8x8_l2);
        decode_block(pbi, xd,
#if CONFIG_SUPERTX
                     supertx_enabled,
#endif
                     mi_row + hbs, mi_col + hbs, r,
                     partition, bsize2, n8x8_l2, n8x8_l2);
        break;
#endif
      default:
        assert(0 && "Invalid partition type");
    }
  }

#if CONFIG_SUPERTX
  if (supertx_enabled && read_token) {
    uint8_t *dst_buf[3];
    int dst_stride[3], i;
    int offset = mi_row * cm->mi_stride + mi_col;

    set_segment_id_supertx(cm, mi_row, mi_col, bsize);

    xd->mi = cm->mi_grid_visible + offset;
    xd->mi[0] = cm->mi + offset;
    set_mi_row_col(xd, tile, mi_row, num_8x8_blocks_high_lookup[bsize],
                   mi_col, num_8x8_blocks_wide_lookup[bsize],
                   cm->mi_rows, cm->mi_cols);
    set_skip_context(xd, mi_row, mi_col);
    skip = read_skip(cm, xd, xd->mi[0]->mbmi.segment_id_supertx, r);
    if (skip) {
      reset_skip_context(xd, bsize);
    } else {
#if CONFIG_EXT_TX
      if (get_ext_tx_types(supertx_size, bsize, 1) > 1) {
        int eset = get_ext_tx_set(supertx_size, bsize, 1);
        if (eset > 0) {
          txfm = vp10_read_tree(r, vp10_ext_tx_inter_tree[eset],
                               cm->fc->inter_ext_tx_prob[eset][supertx_size]);
          if (xd->counts)
            ++xd->counts->inter_ext_tx[eset][supertx_size][txfm];
        }
      }
#else
      if (supertx_size < TX_32X32) {
        txfm = vp10_read_tree(r, vp10_ext_tx_tree,
                             cm->fc->inter_ext_tx_prob[supertx_size]);
        if (xd->counts)
          ++xd->counts->inter_ext_tx[supertx_size][txfm];
      }
#endif  // CONFIG_EXT_TX
    }


    vp10_setup_dst_planes(xd->plane, get_frame_new_buffer(cm), mi_row, mi_col);
    for (i = 0; i < MAX_MB_PLANE; i++) {
      dst_buf[i] = xd->plane[i].dst.buf;
      dst_stride[i] = xd->plane[i].dst.stride;
    }
    dec_predict_sb_complex(pbi, xd, tile, mi_row, mi_col, mi_row, mi_col,
                           bsize, bsize, dst_buf, dst_stride);

    if (!skip) {
      int eobtotal = 0;
      MB_MODE_INFO *mbmi;
      set_offsets_topblock(cm, xd, tile, bsize, mi_row, mi_col);
      mbmi = &xd->mi[0]->mbmi;
      mbmi->tx_type = txfm;
      assert(mbmi->segment_id_supertx != MAX_SEGMENTS);
      for (i = 0; i < MAX_MB_PLANE; ++i) {
        const struct macroblockd_plane *const pd = &xd->plane[i];
        const int num_4x4_w = pd->n4_w;
        const int num_4x4_h = pd->n4_h;
        int row, col;
        const TX_SIZE tx_size =
            i ? dec_get_uv_tx_size(mbmi, pd->n4_wl, pd->n4_hl)
            : mbmi->tx_size;
        const int step = (1 << tx_size);
        const int max_blocks_wide = num_4x4_w +
            (xd->mb_to_right_edge >= 0 ?
             0 : xd->mb_to_right_edge >> (5 + pd->subsampling_x));
        const int max_blocks_high = num_4x4_h +
            (xd->mb_to_bottom_edge >= 0 ?
             0 : xd->mb_to_bottom_edge >> (5 + pd->subsampling_y));

        for (row = 0; row < max_blocks_high; row += step)
          for (col = 0; col < max_blocks_wide; col += step)
            eobtotal += reconstruct_inter_block(xd,
                                                r,
                                                mbmi->segment_id_supertx,
                                                i, row, col,
                                                tx_size);
      }
      if (!(subsize < BLOCK_8X8) && eobtotal == 0)
        skip = 1;
    }
    set_param_topblock(cm, xd, bsize, mi_row, mi_col, txfm, skip);
  }
#endif  // CONFIG_SUPERTX

#if CONFIG_EXT_PARTITION_TYPES
  if (bsize >= BLOCK_8X8) {
    switch (partition) {
      case PARTITION_SPLIT:
        if (bsize > BLOCK_8X8)
          break;
      case PARTITION_NONE:
      case PARTITION_HORZ:
      case PARTITION_VERT:
        update_partition_context(xd, mi_row, mi_col, subsize, bsize);
        break;
      case PARTITION_HORZ_A:
        update_partition_context(xd, mi_row, mi_col, bsize2, subsize);
        update_partition_context(xd, mi_row + hbs, mi_col, subsize, subsize);
        break;
      case PARTITION_HORZ_B:
        update_partition_context(xd, mi_row, mi_col, subsize, subsize);
        update_partition_context(xd, mi_row + hbs, mi_col, bsize2, subsize);
        break;
      case PARTITION_VERT_A:
        update_partition_context(xd, mi_row, mi_col, bsize2, subsize);
        update_partition_context(xd, mi_row, mi_col + hbs, subsize, subsize);
        break;
      case PARTITION_VERT_B:
        update_partition_context(xd, mi_row, mi_col, subsize, subsize);
        update_partition_context(xd, mi_row, mi_col + hbs, bsize2, subsize);
        break;
      default:
        assert(0 && "Invalid partition type");
    }
  }
#else
  // update partition context
  if (bsize >= BLOCK_8X8 &&
      (bsize == BLOCK_8X8 || partition != PARTITION_SPLIT))
    dec_update_partition_context(xd, mi_row, mi_col, subsize, num_8x8_wh);
#endif  // CONFIG_EXT_PARTITION_TYPES
}

#if !CONFIG_ANS
static void setup_bool_decoder(const uint8_t *data,
                               const uint8_t *data_end,
                               const size_t read_size,
                               struct vpx_internal_error_info *error_info,
                               vp10_reader *r,
                               vpx_decrypt_cb decrypt_cb,
                               void *decrypt_state) {
  // Validate the calculated partition length. If the buffer
  // described by the partition can't be fully read, then restrict
  // it to the portion that can be (for EC mode) or throw an error.
  if (!read_is_valid(data, read_size, data_end))
    vpx_internal_error(error_info, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt tile length");

  if (vpx_reader_init(r, data, read_size, decrypt_cb, decrypt_state))
    vpx_internal_error(error_info, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder %d", 1);
}
#else
static void setup_token_decoder(const uint8_t *data,
                                const uint8_t *data_end,
                                const size_t read_size,
                                struct vpx_internal_error_info *error_info,
                                struct AnsDecoder *const ans,
                                vpx_decrypt_cb decrypt_cb,
                                void *decrypt_state) {
  (void) decrypt_cb;
  (void) decrypt_state;
  // Validate the calculated partition length. If the buffer
  // described by the partition can't be fully read, then restrict
  // it to the portion that can be (for EC mode) or throw an error.
  if (!read_is_valid(data, read_size, data_end))
    vpx_internal_error(error_info, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt tile length");

  if (read_size > INT_MAX || ans_read_init(ans, data, (int)read_size))
    vpx_internal_error(error_info, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate token decoder %d", 1);
}
#endif

static void read_coef_probs_common(vp10_coeff_probs_model *coef_probs,
                                   vp10_reader *r) {
  int i, j, k, l, m;

  if (vp10_read_bit(r))
    for (i = 0; i < PLANE_TYPES; ++i)
      for (j = 0; j < REF_TYPES; ++j)
        for (k = 0; k < COEF_BANDS; ++k)
          for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l)
            for (m = 0; m < UNCONSTRAINED_NODES; ++m)
              vp10_diff_update_prob(r, &coef_probs[i][j][k][l][m]);
}

static void read_coef_probs(FRAME_CONTEXT *fc, TX_MODE tx_mode,
                            vp10_reader *r) {
    const TX_SIZE max_tx_size = tx_mode_to_biggest_tx_size[tx_mode];
    TX_SIZE tx_size;
    for (tx_size = TX_4X4; tx_size <= max_tx_size; ++tx_size)
      read_coef_probs_common(fc->coef_probs[tx_size], r);
#if CONFIG_ANS
    vp10_coef_pareto_cdfs(fc);
#endif  // CONFIG_ANS
}

static void setup_segmentation(VP10_COMMON *const cm,
                               struct vpx_read_bit_buffer *rb) {
  struct segmentation *const seg = &cm->seg;
  int i, j;

  seg->update_map = 0;
  seg->update_data = 0;

  seg->enabled = vpx_rb_read_bit(rb);
  if (!seg->enabled)
    return;

  // Segmentation map update
  if (frame_is_intra_only(cm) || cm->error_resilient_mode) {
    seg->update_map = 1;
  } else {
    seg->update_map = vpx_rb_read_bit(rb);
  }
  if (seg->update_map) {
    if (frame_is_intra_only(cm) || cm->error_resilient_mode) {
      seg->temporal_update = 0;
    } else {
      seg->temporal_update = vpx_rb_read_bit(rb);
    }
  }

  // Segmentation data update
  seg->update_data = vpx_rb_read_bit(rb);
  if (seg->update_data) {
    seg->abs_delta = vpx_rb_read_bit(rb);

    vp10_clearall_segfeatures(seg);

    for (i = 0; i < MAX_SEGMENTS; i++) {
      for (j = 0; j < SEG_LVL_MAX; j++) {
        int data = 0;
        const int feature_enabled = vpx_rb_read_bit(rb);
        if (feature_enabled) {
          vp10_enable_segfeature(seg, i, j);
          data = decode_unsigned_max(rb, vp10_seg_feature_data_max(j));
          if (vp10_is_segfeature_signed(j))
            data = vpx_rb_read_bit(rb) ? -data : data;
        }
        vp10_set_segdata(seg, i, j, data);
      }
    }
  }
}

#if CONFIG_LOOP_RESTORATION
static void setup_restoration(VP10_COMMON *cm,
                              struct vpx_read_bit_buffer *rb) {
  RestorationInfo *rst = &cm->rst_info;
  if (vpx_rb_read_bit(rb)) {
    if (vpx_rb_read_bit(rb)) {
      rst->restoration_type = RESTORE_BILATERAL;
      rst->restoration_level =
          vpx_rb_read_literal(rb, vp10_restoration_level_bits(cm));
    } else {
      rst->restoration_type = RESTORE_WIENER;
      rst->vfilter[0] = vpx_rb_read_literal(rb, WIENER_FILT_TAP0_BITS) +
          WIENER_FILT_TAP0_MINV;
      rst->vfilter[1] = vpx_rb_read_literal(rb, WIENER_FILT_TAP1_BITS) +
          WIENER_FILT_TAP1_MINV;
      rst->vfilter[2] = vpx_rb_read_literal(rb, WIENER_FILT_TAP2_BITS) +
          WIENER_FILT_TAP2_MINV;
      rst->hfilter[0] = vpx_rb_read_literal(rb, WIENER_FILT_TAP0_BITS) +
          WIENER_FILT_TAP0_MINV;
      rst->hfilter[1] = vpx_rb_read_literal(rb, WIENER_FILT_TAP1_BITS) +
          WIENER_FILT_TAP1_MINV;
      rst->hfilter[2] = vpx_rb_read_literal(rb, WIENER_FILT_TAP2_BITS) +
          WIENER_FILT_TAP2_MINV;
    }
  } else {
    rst->restoration_type = RESTORE_NONE;
  }
}
#endif  // CONFIG_LOOP_RESTORATION

static void setup_loopfilter(VP10_COMMON *cm,
                             struct vpx_read_bit_buffer *rb) {
  struct loopfilter *lf = &cm->lf;
  lf->filter_level = vpx_rb_read_literal(rb, 6);
  lf->sharpness_level = vpx_rb_read_literal(rb, 3);

  // Read in loop filter deltas applied at the MB level based on mode or ref
  // frame.
  lf->mode_ref_delta_update = 0;

  lf->mode_ref_delta_enabled = vpx_rb_read_bit(rb);
  if (lf->mode_ref_delta_enabled) {
    lf->mode_ref_delta_update = vpx_rb_read_bit(rb);
    if (lf->mode_ref_delta_update) {
      int i;

      for (i = 0; i < MAX_REF_FRAMES; i++)
        if (vpx_rb_read_bit(rb))
          lf->ref_deltas[i] = vpx_rb_read_inv_signed_literal(rb, 6);

      for (i = 0; i < MAX_MODE_LF_DELTAS; i++)
        if (vpx_rb_read_bit(rb))
          lf->mode_deltas[i] = vpx_rb_read_inv_signed_literal(rb, 6);
    }
  }
}

static INLINE int read_delta_q(struct vpx_read_bit_buffer *rb) {
  return vpx_rb_read_bit(rb) ?
      vpx_rb_read_inv_signed_literal(rb, 6) : 0;
}

static void setup_quantization(VP10_COMMON *const cm,
                               struct vpx_read_bit_buffer *rb) {
  cm->base_qindex = vpx_rb_read_literal(rb, QINDEX_BITS);
  cm->y_dc_delta_q = read_delta_q(rb);
  cm->uv_dc_delta_q = read_delta_q(rb);
  cm->uv_ac_delta_q = read_delta_q(rb);
  cm->dequant_bit_depth = cm->bit_depth;
}

static void setup_segmentation_dequant(VP10_COMMON *const cm) {
  // Build y/uv dequant values based on segmentation.
#if CONFIG_NEW_QUANT
  int b;
  int dq;
#endif  //  CONFIG_NEW_QUANT
  if (cm->seg.enabled) {
    int i;
    for (i = 0; i < MAX_SEGMENTS; ++i) {
      const int qindex = vp10_get_qindex(&cm->seg, i, cm->base_qindex);
      cm->y_dequant[i][0] = vp10_dc_quant(qindex, cm->y_dc_delta_q,
                                         cm->bit_depth);
      cm->y_dequant[i][1] = vp10_ac_quant(qindex, 0, cm->bit_depth);
      cm->uv_dequant[i][0] = vp10_dc_quant(qindex, cm->uv_dc_delta_q,
                                          cm->bit_depth);
      cm->uv_dequant[i][1] = vp10_ac_quant(qindex, cm->uv_ac_delta_q,
                                          cm->bit_depth);
#if CONFIG_NEW_QUANT
      for (dq = 0; dq < QUANT_PROFILES; dq ++) {
        for (b = 0; b < COEF_BANDS; ++b) {
          vp10_get_dequant_val_nuq(
              cm->y_dequant[i][b != 0], qindex, b,
              cm->y_dequant_nuq[i][dq][b], NULL, dq);
          vp10_get_dequant_val_nuq(
              cm->uv_dequant[i][b != 0], qindex, b,
              cm->uv_dequant_nuq[i][dq][b], NULL, dq);
        }
      }
#endif  //  CONFIG_NEW_QUANT
    }
  } else {
    const int qindex = cm->base_qindex;
    // When segmentation is disabled, only the first value is used.  The
    // remaining are don't cares.
    cm->y_dequant[0][0] = vp10_dc_quant(qindex, cm->y_dc_delta_q, cm->bit_depth);
    cm->y_dequant[0][1] = vp10_ac_quant(qindex, 0, cm->bit_depth);
    cm->uv_dequant[0][0] = vp10_dc_quant(qindex, cm->uv_dc_delta_q,
                                        cm->bit_depth);
    cm->uv_dequant[0][1] = vp10_ac_quant(qindex, cm->uv_ac_delta_q,
                                        cm->bit_depth);
#if CONFIG_NEW_QUANT
    for (dq = 0; dq < QUANT_PROFILES; dq ++) {
      for (b = 0; b < COEF_BANDS; ++b) {
        vp10_get_dequant_val_nuq(
            cm->y_dequant[0][b != 0], qindex, b,
            cm->y_dequant_nuq[0][dq][b], NULL, dq);
        vp10_get_dequant_val_nuq(
            cm->uv_dequant[0][b != 0], qindex, b,
            cm->uv_dequant_nuq[0][dq][b], NULL, dq);
      }
    }
#endif  //  CONFIG_NEW_QUANT
  }
}

static INTERP_FILTER read_interp_filter(struct vpx_read_bit_buffer *rb) {
  return vpx_rb_read_bit(rb) ?
      SWITCHABLE : vpx_rb_read_literal(rb, 2 + CONFIG_EXT_INTERP);
}

static void setup_render_size(VP10_COMMON *cm,
                              struct vpx_read_bit_buffer *rb) {
  cm->render_width = cm->width;
  cm->render_height = cm->height;
  if (vpx_rb_read_bit(rb))
    vp10_read_frame_size(rb, &cm->render_width, &cm->render_height);
}

static void resize_mv_buffer(VP10_COMMON *cm) {
  vpx_free(cm->cur_frame->mvs);
  cm->cur_frame->mi_rows = cm->mi_rows;
  cm->cur_frame->mi_cols = cm->mi_cols;
  CHECK_MEM_ERROR(cm, cm->cur_frame->mvs,
                  (MV_REF *)vpx_calloc(cm->mi_rows * cm->mi_cols,
                                       sizeof(*cm->cur_frame->mvs)));
}

static void resize_context_buffers(VP10_COMMON *cm, int width, int height) {
#if CONFIG_SIZE_LIMIT
  if (width > DECODE_WIDTH_LIMIT || height > DECODE_HEIGHT_LIMIT)
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Dimensions of %dx%d beyond allowed size of %dx%d.",
                       width, height, DECODE_WIDTH_LIMIT, DECODE_HEIGHT_LIMIT);
#endif
  if (cm->width != width || cm->height != height) {
    const int new_mi_rows =
        ALIGN_POWER_OF_TWO(height, MI_SIZE_LOG2) >> MI_SIZE_LOG2;
    const int new_mi_cols =
        ALIGN_POWER_OF_TWO(width,  MI_SIZE_LOG2) >> MI_SIZE_LOG2;

    // Allocations in vp10_alloc_context_buffers() depend on individual
    // dimensions as well as the overall size.
    if (new_mi_cols > cm->mi_cols || new_mi_rows > cm->mi_rows) {
      if (vp10_alloc_context_buffers(cm, width, height))
        vpx_internal_error(&cm->error, VPX_CODEC_MEM_ERROR,
                           "Failed to allocate context buffers");
    } else {
      vp10_set_mb_mi(cm, width, height);
    }
    vp10_init_context_buffers(cm);
    cm->width = width;
    cm->height = height;
  }
  if (cm->cur_frame->mvs == NULL || cm->mi_rows > cm->cur_frame->mi_rows ||
      cm->mi_cols > cm->cur_frame->mi_cols) {
    resize_mv_buffer(cm);
  }
}

static void setup_frame_size(VP10_COMMON *cm, struct vpx_read_bit_buffer *rb) {
  int width, height;
  BufferPool *const pool = cm->buffer_pool;
  vp10_read_frame_size(rb, &width, &height);
  resize_context_buffers(cm, width, height);
  setup_render_size(cm, rb);

  lock_buffer_pool(pool);
  if (vpx_realloc_frame_buffer(
          get_frame_new_buffer(cm), cm->width, cm->height,
          cm->subsampling_x, cm->subsampling_y,
#if CONFIG_VP9_HIGHBITDEPTH
          cm->use_highbitdepth,
#endif
          VP9_DEC_BORDER_IN_PIXELS,
          cm->byte_alignment,
          &pool->frame_bufs[cm->new_fb_idx].raw_frame_buffer, pool->get_fb_cb,
          pool->cb_priv)) {
    unlock_buffer_pool(pool);
    vpx_internal_error(&cm->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate frame buffer");
  }
  unlock_buffer_pool(pool);

  pool->frame_bufs[cm->new_fb_idx].buf.subsampling_x = cm->subsampling_x;
  pool->frame_bufs[cm->new_fb_idx].buf.subsampling_y = cm->subsampling_y;
  pool->frame_bufs[cm->new_fb_idx].buf.bit_depth = (unsigned int)cm->bit_depth;
  pool->frame_bufs[cm->new_fb_idx].buf.color_space = cm->color_space;
  pool->frame_bufs[cm->new_fb_idx].buf.color_range = cm->color_range;
  pool->frame_bufs[cm->new_fb_idx].buf.render_width  = cm->render_width;
  pool->frame_bufs[cm->new_fb_idx].buf.render_height = cm->render_height;
}

static INLINE int valid_ref_frame_img_fmt(vpx_bit_depth_t ref_bit_depth,
                                          int ref_xss, int ref_yss,
                                          vpx_bit_depth_t this_bit_depth,
                                          int this_xss, int this_yss) {
  return ref_bit_depth == this_bit_depth && ref_xss == this_xss &&
         ref_yss == this_yss;
}

static void setup_frame_size_with_refs(VP10_COMMON *cm,
                                       struct vpx_read_bit_buffer *rb) {
  int width, height;
  int found = 0, i;
  int has_valid_ref_frame = 0;
  BufferPool *const pool = cm->buffer_pool;
  for (i = 0; i < REFS_PER_FRAME; ++i) {
    if (vpx_rb_read_bit(rb)) {
      YV12_BUFFER_CONFIG *const buf = cm->frame_refs[i].buf;
      width = buf->y_crop_width;
      height = buf->y_crop_height;
      cm->render_width = buf->render_width;
      cm->render_height = buf->render_height;
      found = 1;
      break;
    }
  }

  if (!found) {
    vp10_read_frame_size(rb, &width, &height);
    setup_render_size(cm, rb);
  }

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
                         "Referenced frame has incompatible color format");
  }

  resize_context_buffers(cm, width, height);

  lock_buffer_pool(pool);
  if (vpx_realloc_frame_buffer(
          get_frame_new_buffer(cm), cm->width, cm->height,
          cm->subsampling_x, cm->subsampling_y,
#if CONFIG_VP9_HIGHBITDEPTH
          cm->use_highbitdepth,
#endif
          VP9_DEC_BORDER_IN_PIXELS,
          cm->byte_alignment,
          &pool->frame_bufs[cm->new_fb_idx].raw_frame_buffer, pool->get_fb_cb,
          pool->cb_priv)) {
    unlock_buffer_pool(pool);
    vpx_internal_error(&cm->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate frame buffer");
  }
  unlock_buffer_pool(pool);

  pool->frame_bufs[cm->new_fb_idx].buf.subsampling_x = cm->subsampling_x;
  pool->frame_bufs[cm->new_fb_idx].buf.subsampling_y = cm->subsampling_y;
  pool->frame_bufs[cm->new_fb_idx].buf.bit_depth = (unsigned int)cm->bit_depth;
  pool->frame_bufs[cm->new_fb_idx].buf.color_space = cm->color_space;
  pool->frame_bufs[cm->new_fb_idx].buf.color_range = cm->color_range;
  pool->frame_bufs[cm->new_fb_idx].buf.render_width  = cm->render_width;
  pool->frame_bufs[cm->new_fb_idx].buf.render_height = cm->render_height;
}

static void read_tile_info(VP10Decoder *const pbi,
                            struct vpx_read_bit_buffer *const rb) {
  VP10_COMMON *const cm = &pbi->common;
#if CONFIG_EXT_TILE
  // Read the tile width/height
#if CONFIG_EXT_PARTITION
  if (cm->sb_size == BLOCK_128X128) {
    cm->tile_width  = vpx_rb_read_literal(rb, 5) + 1;
    cm->tile_height = vpx_rb_read_literal(rb, 5) + 1;
  } else
#endif  // CONFIG_EXT_PARTITION
  {
    cm->tile_width  = vpx_rb_read_literal(rb, 6) + 1;
    cm->tile_height = vpx_rb_read_literal(rb, 6) + 1;
  }

  cm->tile_width  <<= cm->mib_size_log2;
  cm->tile_height <<= cm->mib_size_log2;

  cm->tile_width  = VPXMIN(cm->tile_width, cm->mi_cols);
  cm->tile_height = VPXMIN(cm->tile_height, cm->mi_rows);

  // Get the number of tiles
  cm->tile_cols = 1;
  while (cm->tile_cols * cm->tile_width < cm->mi_cols)
    ++cm->tile_cols;

  cm->tile_rows = 1;
  while (cm->tile_rows * cm->tile_height < cm->mi_rows)
    ++cm->tile_rows;

  if (cm->tile_cols * cm->tile_rows > 1) {
    // Read the number of bytes used to store tile size
    pbi->tile_col_size_bytes  = vpx_rb_read_literal(rb, 2) + 1;
    pbi->tile_size_bytes = vpx_rb_read_literal(rb, 2) + 1;
  }
#else
  int min_log2_tile_cols, max_log2_tile_cols, max_ones;
  vp10_get_tile_n_bits(cm->mi_cols, &min_log2_tile_cols, &max_log2_tile_cols);

  // columns
  max_ones = max_log2_tile_cols - min_log2_tile_cols;
  cm->log2_tile_cols = min_log2_tile_cols;
  while (max_ones-- && vpx_rb_read_bit(rb))
    cm->log2_tile_cols++;

  if (cm->log2_tile_cols > 6)
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Invalid number of tile columns");

  // rows
  cm->log2_tile_rows = vpx_rb_read_bit(rb);
  if (cm->log2_tile_rows)
    cm->log2_tile_rows += vpx_rb_read_bit(rb);

  cm->tile_cols = 1 << cm->log2_tile_cols;
  cm->tile_rows = 1 << cm->log2_tile_rows;

  cm->tile_width = ALIGN_POWER_OF_TWO(cm->mi_cols, MAX_MIB_SIZE_LOG2);
  cm->tile_width >>= cm->log2_tile_cols;
  cm->tile_height = ALIGN_POWER_OF_TWO(cm->mi_rows, MAX_MIB_SIZE_LOG2);
  cm->tile_height >>= cm->log2_tile_rows;

  // round to integer multiples of superblock size
  cm->tile_width  = ALIGN_POWER_OF_TWO(cm->tile_width, MAX_MIB_SIZE_LOG2);
  cm->tile_height = ALIGN_POWER_OF_TWO(cm->tile_height, MAX_MIB_SIZE_LOG2);

  // tile size magnitude
  if (cm->tile_rows > 1 || cm->tile_cols > 1) {
    pbi->tile_size_bytes = vpx_rb_read_literal(rb, 2) + 1;
  }
#endif  // CONFIG_EXT_TILE
}

static int mem_get_varsize(const uint8_t *src, const int sz) {
  switch (sz) {
    case 1:
      return src[0];
    case 2:
      return mem_get_le16(src);
    case 3:
      return mem_get_le24(src);
    case 4:
      return mem_get_le32(src);
    default:
      assert("Invalid size" && 0);
      return -1;
  }
}

#if CONFIG_EXT_TILE
// Reads the next tile returning its size and adjusting '*data' accordingly
// based on 'is_last'.
static void get_tile_buffer(const uint8_t *const data_end,
                            struct vpx_internal_error_info *error_info,
                            const uint8_t **data,
                            vpx_decrypt_cb decrypt_cb, void *decrypt_state,
                            TileBufferDec (*const tile_buffers)[MAX_TILE_COLS],
                            int tile_size_bytes, int col, int row) {
  size_t size;

  size_t copy_size  = 0;
  const uint8_t *copy_data = NULL;

  if (!read_is_valid(*data, tile_size_bytes, data_end))
    vpx_internal_error(error_info, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt tile length");
  if (decrypt_cb) {
    uint8_t be_data[4];
    decrypt_cb(decrypt_state, *data, be_data, tile_size_bytes);

    // Only read number of bytes in cm->tile_size_bytes.
    size = mem_get_varsize(be_data, tile_size_bytes);
  } else {
    size = mem_get_varsize(*data, tile_size_bytes);
  }

  // The top bit indicates copy mode
  if ((size >> (tile_size_bytes * 8 - 1)) == 1) {
    // The remaining bits in the top byte signal the row offset
    int offset = (size >> (tile_size_bytes - 1) * 8) & 0x7f;

    // Currently, only use tiles in same column as reference tiles.
    copy_data = tile_buffers[row - offset][col].data;
    copy_size = tile_buffers[row - offset][col].size;
    size = 0;
  }

  *data += tile_size_bytes;

  if (size > (size_t)(data_end - *data))
    vpx_internal_error(error_info, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt tile size");

  if (size > 0) {
    tile_buffers[row][col].data = *data;
    tile_buffers[row][col].size = size;
  } else {
    tile_buffers[row][col].data = copy_data;
    tile_buffers[row][col].size = copy_size;
  }

  *data += size;

  tile_buffers[row][col].raw_data_end = *data;
}

static void get_tile_buffers(
    VP10Decoder *pbi,
    const uint8_t *data, const uint8_t *data_end,
    TileBufferDec (*const tile_buffers)[MAX_TILE_COLS]) {
  VP10_COMMON *const cm = &pbi->common;
  const int tile_cols = cm->tile_cols;
  const int tile_rows = cm->tile_rows;
  const int have_tiles = tile_cols * tile_rows > 1;

  if (!have_tiles)  {
    const uint32_t tile_size = data_end - data;
    tile_buffers[0][0].data = data;
    tile_buffers[0][0].size = tile_size;
    tile_buffers[0][0].raw_data_end = NULL;
  } else {
    // We locate only the tile buffers that are required, which are the ones
    // specified by pbi->dec_tile_col and pbi->dec_tile_row. Also, we always
    // need the last (bottom right) tile buffer, as we need to know where the
    // end of the compressed frame buffer is for proper superframe decoding.

    const uint8_t *tile_col_data_end[MAX_TILE_COLS];
    const uint8_t *const data_start = data;

    const int dec_tile_row = VPXMIN(pbi->dec_tile_row, tile_rows);
    const int single_row = pbi->dec_tile_row >= 0;
    const int tile_rows_start = single_row ? dec_tile_row : 0;
    const int tile_rows_end = single_row ? tile_rows_start + 1 : tile_rows;
    const int dec_tile_col = VPXMIN(pbi->dec_tile_col, tile_cols);
    const int single_col = pbi->dec_tile_col >= 0;
    const int tile_cols_start = single_col ? dec_tile_col : 0;
    const int tile_cols_end = single_col ? tile_cols_start + 1 : tile_cols;

    const int tile_col_size_bytes = pbi->tile_col_size_bytes;
    const int tile_size_bytes = pbi->tile_size_bytes;

    size_t tile_col_size;
    int r, c;

    // Read tile column sizes for all columns (we need the last tile buffer)
    for (c = 0; c < tile_cols; ++c) {
      const int is_last = c == tile_cols - 1;
      if (!is_last) {
        tile_col_size = mem_get_varsize(data, tile_col_size_bytes);
        data += tile_col_size_bytes;
        tile_col_data_end[c] = data + tile_col_size;
      } else {
        tile_col_size = data_end - data;
        tile_col_data_end[c] = data_end;
      }
      data += tile_col_size;
    }

    data = data_start;

    // Read the required tile sizes.
    for (c = tile_cols_start; c < tile_cols_end; ++c) {
      const int is_last = c == tile_cols - 1;

      if (c > 0)
        data = tile_col_data_end[c - 1];

      if (!is_last)
        data += tile_col_size_bytes;

      // Get the whole of the last column, otherwise stop at the required tile.
      for (r = 0; r < (is_last ? tile_rows : tile_rows_end); ++r) {
        tile_buffers[r][c].col = c;

        get_tile_buffer(tile_col_data_end[c],
                        &pbi->common.error, &data,
                        pbi->decrypt_cb, pbi->decrypt_state,
                        tile_buffers, tile_size_bytes, c, r);
      }
    }

    // If we have not read the last column, then read it to get the last tile.
    if (tile_cols_end != tile_cols) {
      c = tile_cols - 1;

      data = tile_col_data_end[c - 1];

      for (r = 0; r < tile_rows; ++r) {
        tile_buffers[r][c].col = c;

        get_tile_buffer(tile_col_data_end[c],
                        &pbi->common.error, &data,
                        pbi->decrypt_cb, pbi->decrypt_state,
                        tile_buffers, tile_size_bytes, c, r);
      }
    }
  }
}
#else
// Reads the next tile returning its size and adjusting '*data' accordingly
// based on 'is_last'.
static void get_tile_buffer(const uint8_t *const data_end,
                            const int tile_size_bytes, int is_last,
                            struct vpx_internal_error_info *error_info,
                            const uint8_t **data,
                            vpx_decrypt_cb decrypt_cb, void *decrypt_state,
                            TileBufferDec *const buf) {
  size_t size;

  if (!is_last) {
    if (!read_is_valid(*data, 4, data_end))
      vpx_internal_error(error_info, VPX_CODEC_CORRUPT_FRAME,
                         "Truncated packet or corrupt tile length");

    if (decrypt_cb) {
      uint8_t be_data[4];
      decrypt_cb(decrypt_state, *data, be_data, tile_size_bytes);
      size = mem_get_varsize(be_data, tile_size_bytes);
    } else {
      size = mem_get_varsize(*data, tile_size_bytes);
    }
    *data += tile_size_bytes;

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

static void get_tile_buffers(
    VP10Decoder *pbi,
    const uint8_t *data, const uint8_t *data_end,
    TileBufferDec (*const tile_buffers)[MAX_TILE_COLS]) {
  VP10_COMMON *const cm = &pbi->common;
  int r, c;
  const int tile_cols = cm->tile_cols;
  const int tile_rows = cm->tile_rows;

  for (r = 0; r < tile_rows; ++r) {
    for (c = 0; c < tile_cols; ++c) {
      const int is_last = (r == tile_rows - 1) && (c == tile_cols - 1);
      TileBufferDec *const buf = &tile_buffers[r][c];
      buf->col = c;
      get_tile_buffer(data_end, pbi->tile_size_bytes,
                      is_last, &cm->error, &data,
                      pbi->decrypt_cb, pbi->decrypt_state, buf);
    }
  }
}
#endif  // CONFIG_EXT_TILE

static const uint8_t *decode_tiles(VP10Decoder *pbi,
                                   const uint8_t *data,
                                   const uint8_t *data_end) {
  VP10_COMMON *const cm = &pbi->common;
  const VPxWorkerInterface *const winterface = vpx_get_worker_interface();
  const int tile_cols = cm->tile_cols;
  const int tile_rows = cm->tile_rows;
  const int n_tiles = tile_cols * tile_rows;
  TileBufferDec (*const tile_buffers)[MAX_TILE_COLS] = pbi->tile_buffers;
#if CONFIG_EXT_TILE
  const int dec_tile_row = VPXMIN(pbi->dec_tile_row, tile_rows);
  const int single_row = pbi->dec_tile_row >= 0;
  const int tile_rows_start = single_row ? dec_tile_row : 0;
  const int tile_rows_end = single_row ? dec_tile_row + 1 : tile_rows;
  const int dec_tile_col = VPXMIN(pbi->dec_tile_col, tile_cols);
  const int single_col = pbi->dec_tile_col >= 0;
  const int tile_cols_start = single_col ? dec_tile_col : 0;
  const int tile_cols_end = single_col ? tile_cols_start + 1 : tile_cols;
  const int inv_col_order = pbi->inv_tile_order && !single_col;
  const int inv_row_order = pbi->inv_tile_order && !single_row;
#else
  const int tile_rows_start = 0;
  const int tile_rows_end = tile_rows;
  const int tile_cols_start = 0;
  const int tile_cols_end = tile_cols;
  const int inv_col_order = pbi->inv_tile_order;
  const int inv_row_order = pbi->inv_tile_order;
#endif  // CONFIG_EXT_TILE
  int tile_row, tile_col;

#if CONFIG_ENTROPY
  cm->do_subframe_update = n_tiles == 1;
#endif  // CONFIG_ENTROPY

  if (cm->lf.filter_level && !cm->skip_loop_filter &&
      pbi->lf_worker.data1 == NULL) {
    CHECK_MEM_ERROR(cm, pbi->lf_worker.data1,
                    vpx_memalign(32, sizeof(LFWorkerData)));
    pbi->lf_worker.hook = (VPxWorkerHook)vp10_loop_filter_worker;
    if (pbi->max_threads > 1 && !winterface->reset(&pbi->lf_worker)) {
      vpx_internal_error(&cm->error, VPX_CODEC_ERROR,
                         "Loop filter thread creation failed");
    }
  }

  if (cm->lf.filter_level && !cm->skip_loop_filter) {
    LFWorkerData *const lf_data = (LFWorkerData*)pbi->lf_worker.data1;
    // Be sure to sync as we might be resuming after a failed frame decode.
    winterface->sync(&pbi->lf_worker);
    vp10_loop_filter_data_reset(lf_data, get_frame_new_buffer(cm), cm,
                               pbi->mb.plane);
  }

  assert(tile_rows <= MAX_TILE_ROWS);
  assert(tile_cols <= MAX_TILE_COLS);

  get_tile_buffers(pbi, data, data_end, tile_buffers);

  if (pbi->tile_data == NULL || n_tiles != pbi->allocated_tiles) {
    vpx_free(pbi->tile_data);
    CHECK_MEM_ERROR(
        cm,
        pbi->tile_data,
        vpx_memalign(32, n_tiles * (sizeof(*pbi->tile_data))));
    pbi->allocated_tiles = n_tiles;
  }

  // Load all tile information into tile_data.
  for (tile_row = tile_rows_start; tile_row < tile_rows_end; ++tile_row) {
    for (tile_col = tile_cols_start; tile_col < tile_cols_end; ++tile_col) {
      const TileBufferDec *const buf = &tile_buffers[tile_row][tile_col];
      TileData *const td = pbi->tile_data + tile_cols * tile_row + tile_col;

      td->cm = cm;
      td->xd = pbi->mb;
      td->xd.corrupted = 0;
      td->xd.counts =
          cm->refresh_frame_context == REFRESH_FRAME_CONTEXT_BACKWARD ?
              &cm->counts : NULL;
      vp10_zero(td->dqcoeff);
      vp10_tile_init(&td->xd.tile, td->cm, tile_row, tile_col);
#if !CONFIG_ANS
      setup_bool_decoder(buf->data, data_end, buf->size, &cm->error,
                         &td->bit_reader, pbi->decrypt_cb,
                         pbi->decrypt_state);
#else
      setup_token_decoder(buf->data, data_end, buf->size, &cm->error,
                          &td->bit_reader, pbi->decrypt_cb,
                          pbi->decrypt_state);
#endif
      vp10_init_macroblockd(cm, &td->xd, td->dqcoeff);
      td->xd.plane[0].color_index_map = td->color_index_map[0];
      td->xd.plane[1].color_index_map = td->color_index_map[1];
    }
  }

  for (tile_row = tile_rows_start; tile_row < tile_rows_end; ++tile_row) {
    const int row = inv_row_order ? tile_rows - 1 - tile_row : tile_row;
    int mi_row = 0;
    TileInfo tile_info;

    vp10_tile_set_row(&tile_info, cm, row);

    for (tile_col = tile_cols_start; tile_col < tile_cols_end; ++tile_col) {
      const int col = inv_col_order ? tile_cols - 1 - tile_col : tile_col;
      TileData *const td = pbi->tile_data + tile_cols * row + col;

      vp10_tile_set_col(&tile_info, cm, col);

      vp10_zero_above_context(cm, tile_info.mi_col_start, tile_info.mi_col_end);

      for (mi_row = tile_info.mi_row_start; mi_row < tile_info.mi_row_end;
           mi_row += cm->mib_size) {
        int mi_col;

        vp10_zero_left_context(&td->xd);

        for (mi_col = tile_info.mi_col_start; mi_col < tile_info.mi_col_end;
             mi_col += cm->mib_size) {
          decode_partition(pbi, &td->xd,
#if CONFIG_SUPERTX
                           0,
#endif  // CONFIG_SUPERTX
                           mi_row, mi_col, &td->bit_reader,
                           cm->sb_size, b_width_log2_lookup[cm->sb_size]);
        }
        pbi->mb.corrupted |= td->xd.corrupted;
        if (pbi->mb.corrupted)
            vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                               "Failed to decode tile data");
#if CONFIG_ENTROPY
        if (cm->do_subframe_update &&
            cm->refresh_frame_context == REFRESH_FRAME_CONTEXT_BACKWARD) {
          if ((mi_row + MI_SIZE) % (MI_SIZE *
              VPXMAX(cm->mi_rows / MI_SIZE / COEF_PROBS_BUFS, 1)) == 0 &&
              mi_row + MI_SIZE < cm->mi_rows &&
              cm->coef_probs_update_idx < COEF_PROBS_BUFS - 1) {
            vp10_partial_adapt_probs(cm, mi_row, mi_col);
            ++cm->coef_probs_update_idx;
          }
        }
#endif  // CONFIG_ENTROPY
      }
    }

    assert(mi_row > 0);

#if !CONFIG_VAR_TX
    // Loopfilter one tile row.
    if (cm->lf.filter_level && !cm->skip_loop_filter) {
      LFWorkerData *const lf_data = (LFWorkerData*)pbi->lf_worker.data1;
      const int lf_start = VPXMAX(0, tile_info.mi_row_start - cm->mib_size);
      const int lf_end = tile_info.mi_row_end - cm->mib_size;

      // Delay the loopfilter if the first tile row is only
      // a single superblock high.
      if (lf_end <= 0)
        continue;

      // Decoding has completed. Finish up the loop filter in this thread.
      if (tile_info.mi_row_end >= cm->mi_rows)
        continue;

      winterface->sync(&pbi->lf_worker);
      lf_data->start = lf_start;
      lf_data->stop = lf_end;
      if (pbi->max_threads > 1) {
        winterface->launch(&pbi->lf_worker);
      } else {
        winterface->execute(&pbi->lf_worker);
      }
    }

    // After loopfiltering, the last 7 row pixels in each superblock row may
    // still be changed by the longest loopfilter of the next superblock row.
    if (cm->frame_parallel_decode)
      vp10_frameworker_broadcast(pbi->cur_buf, mi_row << cm->mib_size_log2);
#endif  // !CONFIG_VAR_TX
  }

#if CONFIG_VAR_TX
  // Loopfilter the whole frame.
  vp10_loop_filter_frame(get_frame_new_buffer(cm), cm, &pbi->mb,
                         cm->lf.filter_level, 0, 0);
#else
  // Loopfilter remaining rows in the frame.
  if (cm->lf.filter_level && !cm->skip_loop_filter) {
    LFWorkerData *const lf_data = (LFWorkerData*)pbi->lf_worker.data1;
    winterface->sync(&pbi->lf_worker);
    lf_data->start = lf_data->stop;
    lf_data->stop = cm->mi_rows;
    winterface->execute(&pbi->lf_worker);
  }
#endif  // CONFIG_VAR_TX

  if (cm->frame_parallel_decode)
    vp10_frameworker_broadcast(pbi->cur_buf, INT_MAX);

#if CONFIG_EXT_TILE
  if (n_tiles == 1) {
#if CONFIG_ANS
    return data_end;
#else
    // Find the end of the single tile buffer
    return vpx_reader_find_end(&pbi->tile_data->bit_reader);
#endif  // CONFIG_ANS
  } else {
    // Return the end of the last tile buffer
    return tile_buffers[tile_rows - 1][tile_cols - 1].raw_data_end;
  }
#else
#if CONFIG_ANS
  return data_end;
#else
  {
    // Get last tile data.
    TileData *const td = pbi->tile_data + tile_cols * tile_rows - 1;
    return vpx_reader_find_end(&td->bit_reader);
  }
#endif  // CONFIG_ANS
#endif  // CONFIG_EXT_TILE
}

static int tile_worker_hook(TileWorkerData *const tile_data,
                            const TileInfo *const tile) {
  VP10Decoder *const pbi = tile_data->pbi;
  const VP10_COMMON *const cm = &pbi->common;
  int mi_row, mi_col;

  if (setjmp(tile_data->error_info.jmp)) {
    tile_data->error_info.setjmp = 0;
    tile_data->xd.corrupted = 1;
    return 0;
  }

  tile_data->error_info.setjmp = 1;
  tile_data->xd.error_info = &tile_data->error_info;

  vp10_zero_above_context(&pbi->common, tile->mi_col_start, tile->mi_col_end);

  for (mi_row = tile->mi_row_start; mi_row < tile->mi_row_end;
       mi_row += cm->mib_size) {
    vp10_zero_left_context(&tile_data->xd);

    for (mi_col = tile->mi_col_start; mi_col < tile->mi_col_end;
         mi_col += cm->mib_size) {
      decode_partition(pbi, &tile_data->xd,
#if CONFIG_SUPERTX
                       0,
#endif
                       mi_row, mi_col, &tile_data->bit_reader,
                       cm->sb_size, b_width_log2_lookup[cm->sb_size]);
    }
  }
  return !tile_data->xd.corrupted;
}

// sorts in descending order
static int compare_tile_buffers(const void *a, const void *b) {
  const TileBufferDec *const buf1 = (const TileBufferDec*)a;
  const TileBufferDec *const buf2 = (const TileBufferDec*)b;
  return (int)(buf2->size - buf1->size);
}

static const uint8_t *decode_tiles_mt(VP10Decoder *pbi,
                                      const uint8_t *data,
                                      const uint8_t *data_end) {
  VP10_COMMON *const cm = &pbi->common;
  const VPxWorkerInterface *const winterface = vpx_get_worker_interface();
  const int tile_cols = cm->tile_cols;
  const int tile_rows = cm->tile_rows;
  const int num_workers = VPXMIN(pbi->max_threads & ~1, tile_cols);
  TileBufferDec (*const tile_buffers)[MAX_TILE_COLS] = pbi->tile_buffers;
#if CONFIG_EXT_TILE
  const int dec_tile_row = VPXMIN(pbi->dec_tile_row, tile_rows);
  const int single_row = pbi->dec_tile_row >= 0;
  const int tile_rows_start = single_row ? dec_tile_row : 0;
  const int tile_rows_end = single_row ? dec_tile_row + 1 : tile_rows;
  const int dec_tile_col = VPXMIN(pbi->dec_tile_col, tile_cols);
  const int single_col = pbi->dec_tile_col >= 0;
  const int tile_cols_start = single_col ? dec_tile_col : 0;
  const int tile_cols_end = single_col ? tile_cols_start + 1 : tile_cols;
#else
  const int tile_rows_start = 0;
  const int tile_rows_end = tile_rows;
  const int tile_cols_start = 0;
  const int tile_cols_end = tile_cols;
#endif  // CONFIG_EXT_TILE
  int tile_row, tile_col;
  int i;

#if !(CONFIG_ANS || CONFIG_EXT_TILE)
  int final_worker = -1;
#endif  // !(CONFIG_ANS || CONFIG_EXT_TILE)

  assert(tile_rows <= MAX_TILE_ROWS);
  assert(tile_cols <= MAX_TILE_COLS);

  assert(tile_cols * tile_rows > 1);

#if CONFIG_ANS
  // TODO(any): This might just work now. Needs to be tested.
  abort();  // FIXME: Tile parsing broken
#endif  // CONFIG_ANS

  // TODO(jzern): See if we can remove the restriction of passing in max
  // threads to the decoder.
  if (pbi->num_tile_workers == 0) {
    const int num_threads = pbi->max_threads & ~1;
    CHECK_MEM_ERROR(cm, pbi->tile_workers,
                    vpx_malloc(num_threads * sizeof(*pbi->tile_workers)));
    // Ensure tile data offsets will be properly aligned. This may fail on
    // platforms without DECLARE_ALIGNED().
    assert((sizeof(*pbi->tile_worker_data) % 16) == 0);
    CHECK_MEM_ERROR(cm, pbi->tile_worker_data,
                    vpx_memalign(32, num_threads *
                                 sizeof(*pbi->tile_worker_data)));
    CHECK_MEM_ERROR(cm, pbi->tile_worker_info,
                    vpx_malloc(num_threads * sizeof(*pbi->tile_worker_info)));
    for (i = 0; i < num_threads; ++i) {
      VPxWorker *const worker = &pbi->tile_workers[i];
      ++pbi->num_tile_workers;

      winterface->init(worker);
      if (i < num_threads - 1 && !winterface->reset(worker)) {
        vpx_internal_error(&cm->error, VPX_CODEC_ERROR,
                           "Tile decoder thread creation failed");
      }
    }
  }

  // Reset tile decoding hook
  for (i = 0; i < num_workers; ++i) {
    VPxWorker *const worker = &pbi->tile_workers[i];
    winterface->sync(worker);
    worker->hook = (VPxWorkerHook)tile_worker_hook;
    worker->data1 = &pbi->tile_worker_data[i];
    worker->data2 = &pbi->tile_worker_info[i];
  }

  // Initialize thread frame counts.
  if (cm->refresh_frame_context == REFRESH_FRAME_CONTEXT_BACKWARD) {
    for (i = 0; i < num_workers; ++i) {
      TileWorkerData *const twd = (TileWorkerData*)pbi->tile_workers[i].data1;
      vp10_zero(twd->counts);
    }
  }

  // Load tile data into tile_buffers
  get_tile_buffers(pbi, data, data_end, tile_buffers);

  for (tile_row = tile_rows_start ; tile_row < tile_rows_end ; ++tile_row) {
    // Sort the buffers in this tile row based on size in descending order.
    qsort(&tile_buffers[tile_row][tile_cols_start],
          tile_cols_end - tile_cols_start, sizeof(tile_buffers[0][0]),
          compare_tile_buffers);

    // Rearrange the tile buffers in this tile row such that per-tile group
    // the largest, and presumably the most difficult tile will be decoded in
    // the main thread. This should help minimize the number of instances
    // where the main thread is waiting for a worker to complete.
    {
      int group_start;
      for (group_start = tile_cols_start ; group_start < tile_cols_end ;
           group_start += num_workers) {
        const int group_end = VPXMIN(group_start + num_workers, tile_cols);
        const TileBufferDec largest = tile_buffers[tile_row][group_start];
        memmove(&tile_buffers[tile_row][group_start],
                &tile_buffers[tile_row][group_start + 1],
                (group_end - group_start - 1) * sizeof(tile_buffers[0][0]));
        tile_buffers[tile_row][group_end - 1] = largest;
      }
    }

    for (tile_col = tile_cols_start ; tile_col < tile_cols_end ; ) {
      // Launch workers for individual columns
      for (i = 0; i < num_workers && tile_col < tile_cols_end;
           ++i, ++tile_col) {
        TileBufferDec *const buf = &tile_buffers[tile_row][tile_col];
        VPxWorker *const worker = &pbi->tile_workers[i];
        TileWorkerData *const twd = (TileWorkerData*)worker->data1;
        TileInfo *const tile_info = (TileInfo*)worker->data2;

        twd->pbi = pbi;
        twd->xd = pbi->mb;
        twd->xd.corrupted = 0;
        twd->xd.counts =
            cm->refresh_frame_context == REFRESH_FRAME_CONTEXT_BACKWARD ?
                &twd->counts : NULL;
        vp10_zero(twd->dqcoeff);
        vp10_tile_init(tile_info, cm, tile_row, buf->col);
        vp10_tile_init(&twd->xd.tile, cm, tile_row, buf->col);
#if !CONFIG_ANS
        setup_bool_decoder(buf->data, data_end, buf->size, &cm->error,
                           &twd->bit_reader,
                           pbi->decrypt_cb, pbi->decrypt_state);
#else
        setup_token_decoder(buf->data, data_end, buf->size, &cm->error,
                            &twd->bit_reader, pbi->decrypt_cb,
                            pbi->decrypt_state);
#endif  // CONFIG_ANS
        vp10_init_macroblockd(cm, &twd->xd, twd->dqcoeff);
        twd->xd.plane[0].color_index_map = twd->color_index_map[0];
        twd->xd.plane[1].color_index_map = twd->color_index_map[1];

        worker->had_error = 0;
        if (i == num_workers - 1 || tile_col == tile_cols_end - 1) {
          winterface->execute(worker);
        } else {
          winterface->launch(worker);
        }

#if !(CONFIG_ANS || CONFIG_EXT_TILE)
        if (tile_row == tile_rows - 1 && buf->col == tile_cols - 1) {
          final_worker = i;
        }
#endif  // !(CONFIG_ANS || CONFIG_EXT_TILE)
      }

      // Sync all workers
      for (; i > 0; --i) {
        VPxWorker *const worker = &pbi->tile_workers[i - 1];
        // TODO(jzern): The tile may have specific error data associated with
        // its vpx_internal_error_info which could be propagated to the main
        // info in cm. Additionally once the threads have been synced and an
        // error is detected, there's no point in continuing to decode tiles.
        pbi->mb.corrupted |= !winterface->sync(worker);
      }
    }
  }

  // Accumulate thread frame counts.
  if (cm->refresh_frame_context == REFRESH_FRAME_CONTEXT_BACKWARD) {
    for (i = 0; i < num_workers; ++i) {
      TileWorkerData *const twd = (TileWorkerData*)pbi->tile_workers[i].data1;
      vp10_accumulate_frame_counts(cm, &twd->counts);
    }
  }

#if CONFIG_EXT_TILE
  // Return the end of the last tile buffer
  return tile_buffers[tile_rows - 1][tile_cols - 1].raw_data_end;
#else
#if CONFIG_ANS
  return data_end;
#else
  assert(final_worker != -1);
  {
    TileWorkerData *const twd =
        (TileWorkerData*)pbi->tile_workers[final_worker].data1;
    return vpx_reader_find_end(&twd->bit_reader);
  }
#endif  // CONFIG_ANS
#endif  // CONFIG_EXT_TILE
}

static void error_handler(void *data) {
  VP10_COMMON *const cm = (VP10_COMMON *)data;
  vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME, "Truncated packet");
}

static void read_bitdepth_colorspace_sampling(
    VP10_COMMON *cm, struct vpx_read_bit_buffer *rb) {
  if (cm->profile >= PROFILE_2) {
    cm->bit_depth = vpx_rb_read_bit(rb) ? VPX_BITS_12 : VPX_BITS_10;
#if CONFIG_VP9_HIGHBITDEPTH
    cm->use_highbitdepth = 1;
#endif
  } else {
    cm->bit_depth = VPX_BITS_8;
#if CONFIG_VP9_HIGHBITDEPTH
    cm->use_highbitdepth = 0;
#endif
  }
  cm->color_space = vpx_rb_read_literal(rb, 3);
  if (cm->color_space != VPX_CS_SRGB) {
    // [16,235] (including xvycc) vs [0,255] range
    cm->color_range = vpx_rb_read_bit(rb);
    if (cm->profile == PROFILE_1 || cm->profile == PROFILE_3) {
      cm->subsampling_x = vpx_rb_read_bit(rb);
      cm->subsampling_y = vpx_rb_read_bit(rb);
      if (cm->subsampling_x == 1 && cm->subsampling_y == 1)
        vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                           "4:2:0 color not supported in profile 1 or 3");
      if (vpx_rb_read_bit(rb))
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
      if (vpx_rb_read_bit(rb))
        vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                           "Reserved bit set");
    } else {
      vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                         "4:4:4 color not supported in profile 0 or 2");
    }
  }
}

static size_t read_uncompressed_header(VP10Decoder *pbi,
                                       struct vpx_read_bit_buffer *rb) {
  VP10_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  BufferPool *const pool = cm->buffer_pool;
  RefCntBuffer *const frame_bufs = pool->frame_bufs;
  int i, mask, ref_index = 0;
  size_t sz;
#if CONFIG_EXT_REFS
  cm->last3_frame_type = cm->last2_frame_type;
  cm->last2_frame_type = cm->last_frame_type;
#endif  // CONFIG_EXT_REFS
  cm->last_frame_type = cm->frame_type;
  cm->last_intra_only = cm->intra_only;

#if CONFIG_EXT_REFS
  // NOTE: By default all coded frames to be used as a reference
  cm->is_reference_frame = 1;
#endif  // CONFIG_EXT_REFS

  if (vpx_rb_read_literal(rb, 2) != VP9_FRAME_MARKER)
      vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                         "Invalid frame marker");

  cm->profile = vp10_read_profile(rb);
#if CONFIG_VP9_HIGHBITDEPTH
  if (cm->profile >= MAX_PROFILES)
    vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                       "Unsupported bitstream profile");
#else
  if (cm->profile >= PROFILE_2)
    vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                       "Unsupported bitstream profile");
#endif

  cm->show_existing_frame = vpx_rb_read_bit(rb);

  if (cm->show_existing_frame) {
    // Show an existing frame directly.
    const int frame_to_show = cm->ref_frame_map[vpx_rb_read_literal(rb, 3)];

    lock_buffer_pool(pool);
    if (frame_to_show < 0 || frame_bufs[frame_to_show].ref_count < 1) {
      unlock_buffer_pool(pool);
      vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                         "Buffer %d does not contain a decoded frame",
                         frame_to_show);
    }
    ref_cnt_fb(frame_bufs, &cm->new_fb_idx, frame_to_show);
    unlock_buffer_pool(pool);

    cm->lf.filter_level = 0;
    cm->show_frame = 1;
    pbi->refresh_frame_flags = 0;

    if (cm->frame_parallel_decode) {
      for (i = 0; i < REF_FRAMES; ++i)
        cm->next_ref_frame_map[i] = cm->ref_frame_map[i];
    }

    return 0;
  }

  cm->frame_type = (FRAME_TYPE) vpx_rb_read_bit(rb);
  cm->show_frame = vpx_rb_read_bit(rb);
  cm->error_resilient_mode = vpx_rb_read_bit(rb);

  if (cm->frame_type == KEY_FRAME) {
    if (!vp10_read_sync_code(rb))
      vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                         "Invalid frame sync code");

    read_bitdepth_colorspace_sampling(cm, rb);
    pbi->refresh_frame_flags = (1 << REF_FRAMES) - 1;

    for (i = 0; i < REFS_PER_FRAME; ++i) {
      cm->frame_refs[i].idx = INVALID_IDX;
      cm->frame_refs[i].buf = NULL;
    }

    setup_frame_size(cm, rb);
    if (pbi->need_resync) {
      memset(&cm->ref_frame_map, -1, sizeof(cm->ref_frame_map));
      pbi->need_resync = 0;
    }
    if (frame_is_intra_only(cm))
      cm->allow_screen_content_tools = vpx_rb_read_bit(rb);
  } else {
    cm->intra_only = cm->show_frame ? 0 : vpx_rb_read_bit(rb);

    if (cm->error_resilient_mode) {
        cm->reset_frame_context = RESET_FRAME_CONTEXT_ALL;
    } else {
      if (cm->intra_only) {
          cm->reset_frame_context =
              vpx_rb_read_bit(rb) ? RESET_FRAME_CONTEXT_ALL
                                  : RESET_FRAME_CONTEXT_CURRENT;
      } else {
          cm->reset_frame_context =
              vpx_rb_read_bit(rb) ? RESET_FRAME_CONTEXT_CURRENT
                                  : RESET_FRAME_CONTEXT_NONE;
          if (cm->reset_frame_context == RESET_FRAME_CONTEXT_CURRENT)
            cm->reset_frame_context =
                  vpx_rb_read_bit(rb) ? RESET_FRAME_CONTEXT_ALL
                                      : RESET_FRAME_CONTEXT_CURRENT;
      }
    }

    if (cm->intra_only) {
      if (!vp10_read_sync_code(rb))
        vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                           "Invalid frame sync code");

      read_bitdepth_colorspace_sampling(cm, rb);

      pbi->refresh_frame_flags = vpx_rb_read_literal(rb, REF_FRAMES);
      setup_frame_size(cm, rb);
      if (pbi->need_resync) {
        memset(&cm->ref_frame_map, -1, sizeof(cm->ref_frame_map));
        pbi->need_resync = 0;
      }
    } else if (pbi->need_resync != 1) {  /* Skip if need resync */
      pbi->refresh_frame_flags = vpx_rb_read_literal(rb, REF_FRAMES);

#if CONFIG_EXT_REFS
      if (!pbi->refresh_frame_flags) {
        // NOTE: "pbi->refresh_frame_flags == 0" indicates that the coded frame
        //       will not be used as a reference
        cm->is_reference_frame = 0;
      }
#endif  // CONFIG_EXT_REFS

      for (i = 0; i < REFS_PER_FRAME; ++i) {
        const int ref = vpx_rb_read_literal(rb, REF_FRAMES_LOG2);
        const int idx = cm->ref_frame_map[ref];
        RefBuffer *const ref_frame = &cm->frame_refs[i];
        ref_frame->idx = idx;
        ref_frame->buf = &frame_bufs[idx].buf;
        cm->ref_frame_sign_bias[LAST_FRAME + i] = vpx_rb_read_bit(rb);
      }

      setup_frame_size_with_refs(cm, rb);

      cm->allow_high_precision_mv = vpx_rb_read_bit(rb);
      cm->interp_filter = read_interp_filter(rb);

      for (i = 0; i < REFS_PER_FRAME; ++i) {
        RefBuffer *const ref_buf = &cm->frame_refs[i];
#if CONFIG_VP9_HIGHBITDEPTH
        vp10_setup_scale_factors_for_frame(&ref_buf->sf,
                                          ref_buf->buf->y_crop_width,
                                          ref_buf->buf->y_crop_height,
                                          cm->width, cm->height,
                                          cm->use_highbitdepth);
#else
        vp10_setup_scale_factors_for_frame(&ref_buf->sf,
                                          ref_buf->buf->y_crop_width,
                                          ref_buf->buf->y_crop_height,
                                          cm->width, cm->height);
#endif
      }
    }
  }
#if CONFIG_VP9_HIGHBITDEPTH
  get_frame_new_buffer(cm)->bit_depth = cm->bit_depth;
#endif
  get_frame_new_buffer(cm)->color_space = cm->color_space;
  get_frame_new_buffer(cm)->color_range = cm->color_range;
  get_frame_new_buffer(cm)->render_width  = cm->render_width;
  get_frame_new_buffer(cm)->render_height = cm->render_height;

  if (pbi->need_resync) {
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Keyframe / intra-only frame required to reset decoder"
                       " state");
  }

  if (!cm->error_resilient_mode) {
    cm->refresh_frame_context =
        vpx_rb_read_bit(rb) ? REFRESH_FRAME_CONTEXT_FORWARD
                            : REFRESH_FRAME_CONTEXT_BACKWARD;
  } else {
    cm->refresh_frame_context = REFRESH_FRAME_CONTEXT_FORWARD;
  }

  // This flag will be overridden by the call to vp10_setup_past_independence
  // below, forcing the use of context 0 for those frame types.
  cm->frame_context_idx = vpx_rb_read_literal(rb, FRAME_CONTEXTS_LOG2);

  // Generate next_ref_frame_map.
  lock_buffer_pool(pool);
  for (mask = pbi->refresh_frame_flags; mask; mask >>= 1) {
    if (mask & 1) {
      cm->next_ref_frame_map[ref_index] = cm->new_fb_idx;
      ++frame_bufs[cm->new_fb_idx].ref_count;
    } else {
      cm->next_ref_frame_map[ref_index] = cm->ref_frame_map[ref_index];
    }
    // Current thread holds the reference frame.
    if (cm->ref_frame_map[ref_index] >= 0)
      ++frame_bufs[cm->ref_frame_map[ref_index]].ref_count;
    ++ref_index;
  }

  for (; ref_index < REF_FRAMES; ++ref_index) {
    cm->next_ref_frame_map[ref_index] = cm->ref_frame_map[ref_index];

    // Current thread holds the reference frame.
    if (cm->ref_frame_map[ref_index] >= 0)
      ++frame_bufs[cm->ref_frame_map[ref_index]].ref_count;
  }
  unlock_buffer_pool(pool);
  pbi->hold_ref_buf = 1;

  if (frame_is_intra_only(cm) || cm->error_resilient_mode)
    vp10_setup_past_independence(cm);

#if CONFIG_EXT_PARTITION
  set_sb_size(cm, vpx_rb_read_bit(rb) ? BLOCK_128X128 : BLOCK_64X64);
#else
  set_sb_size(cm, BLOCK_64X64);
#endif  // CONFIG_EXT_PARTITION

  setup_loopfilter(cm, rb);
#if CONFIG_LOOP_RESTORATION
  setup_restoration(cm, rb);
#endif  // CONFIG_LOOP_RESTORATION
  setup_quantization(cm, rb);
#if CONFIG_VP9_HIGHBITDEPTH
  xd->bd = (int)cm->bit_depth;
#endif

#if CONFIG_ENTROPY
  vp10_default_coef_probs(cm);
  if (cm->frame_type == KEY_FRAME || cm->error_resilient_mode ||
      cm->reset_frame_context == RESET_FRAME_CONTEXT_ALL) {
    for (i = 0; i < FRAME_CONTEXTS; ++i)
      cm->frame_contexts[i] = *cm->fc;
  } else if (cm->reset_frame_context == RESET_FRAME_CONTEXT_CURRENT) {
    cm->frame_contexts[cm->frame_context_idx] = *cm->fc;
  }
#endif  // CONFIG_ENTROPY

  setup_segmentation(cm, rb);

  {
    int i;
    for (i = 0; i < MAX_SEGMENTS; ++i) {
      const int qindex = cm->seg.enabled ?
          vp10_get_qindex(&cm->seg, i, cm->base_qindex) : cm->base_qindex;
      xd->lossless[i] = qindex == 0 &&
          cm->y_dc_delta_q == 0 &&
          cm->uv_dc_delta_q == 0 &&
          cm->uv_ac_delta_q == 0;
    }
  }

  setup_segmentation_dequant(cm);
  cm->tx_mode = (!cm->seg.enabled && xd->lossless[0]) ? ONLY_4X4
                                                      : read_tx_mode(rb);
  cm->reference_mode = read_frame_reference_mode(cm, rb);

  read_tile_info(pbi, rb);
  sz = vpx_rb_read_literal(rb, 16);

  if (sz == 0)
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Invalid header size");

  return sz;
}

#if CONFIG_EXT_TX
static void read_ext_tx_probs(FRAME_CONTEXT *fc, vp10_reader *r) {
  int i, j, k;
  int s;
  for (s = 1; s < EXT_TX_SETS_INTER; ++s) {
    if (vp10_read(r, GROUP_DIFF_UPDATE_PROB)) {
      for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
        if (!use_inter_ext_tx_for_txsize[s][i]) continue;
        for (j = 0; j < num_ext_tx_set_inter[s] - 1; ++j)
          vp10_diff_update_prob(r, &fc->inter_ext_tx_prob[s][i][j]);
      }
    }
  }

  for (s = 1; s < EXT_TX_SETS_INTRA; ++s) {
    if (vp10_read(r, GROUP_DIFF_UPDATE_PROB)) {
      for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
        if (!use_intra_ext_tx_for_txsize[s][i]) continue;
        for (j = 0; j < INTRA_MODES; ++j)
          for (k = 0; k < num_ext_tx_set_intra[s] - 1; ++k)
            vp10_diff_update_prob(r, &fc->intra_ext_tx_prob[s][i][j][k]);
      }
    }
  }
}

#else

static void read_ext_tx_probs(FRAME_CONTEXT *fc, vp10_reader *r) {
  int i, j, k;
  if (vp10_read(r, GROUP_DIFF_UPDATE_PROB)) {
    for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
      for (j = 0; j < TX_TYPES; ++j)
        for (k = 0; k < TX_TYPES - 1; ++k)
          vp10_diff_update_prob(r, &fc->intra_ext_tx_prob[i][j][k]);
    }
  }
  if (vp10_read(r, GROUP_DIFF_UPDATE_PROB)) {
    for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
      for (k = 0; k < TX_TYPES - 1; ++k)
        vp10_diff_update_prob(r, &fc->inter_ext_tx_prob[i][k]);
    }
  }
}
#endif  // CONFIG_EXT_TX

#if CONFIG_SUPERTX
static void read_supertx_probs(FRAME_CONTEXT *fc, vp10_reader *r) {
  int i, j;
  if (vp10_read(r, GROUP_DIFF_UPDATE_PROB)) {
    for (i = 0; i < PARTITION_SUPERTX_CONTEXTS; ++i) {
      for (j = 1; j < TX_SIZES; ++j) {
        vp10_diff_update_prob(r, &fc->supertx_prob[i][j]);
      }
    }
  }
}
#endif  // CONFIG_SUPERTX

static int read_compressed_header(VP10Decoder *pbi, const uint8_t *data,
                                  size_t partition_size) {
  VP10_COMMON *const cm = &pbi->common;
#if CONFIG_SUPERTX
  MACROBLOCKD *const xd = &pbi->mb;
#endif
  FRAME_CONTEXT *const fc = cm->fc;
  vp10_reader r;
  int k, i, j;

#if !CONFIG_ANS
  if (vpx_reader_init(&r, data, partition_size, pbi->decrypt_cb,
                      pbi->decrypt_state))
    vpx_internal_error(&cm->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder 0");
#else
  if (ans_read_init(&r, data, partition_size))
    vpx_internal_error(&cm->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate compressed header ANS decoder");
#endif  // !CONFIG_ANS

  if (cm->tx_mode == TX_MODE_SELECT) {
    for (i = 0; i < TX_SIZES - 1; ++i)
      for (j = 0; j < TX_SIZE_CONTEXTS; ++j)
        for (k = 0; k < i + 1; ++k)
          vp10_diff_update_prob(&r, &fc->tx_size_probs[i][j][k]);
  }

  read_coef_probs(fc, cm->tx_mode, &r);

#if CONFIG_VAR_TX
  for (k = 0; k < TXFM_PARTITION_CONTEXTS; ++k)
    vp10_diff_update_prob(&r, &fc->txfm_partition_prob[k]);
#endif

  for (k = 0; k < SKIP_CONTEXTS; ++k)
    vp10_diff_update_prob(&r, &fc->skip_probs[k]);

  if (cm->seg.enabled) {
    if (cm->seg.temporal_update) {
      for (k = 0; k < PREDICTION_PROBS; k++)
        vp10_diff_update_prob(&r, &cm->fc->seg.pred_probs[k]);
    }
    for (k = 0; k < MAX_SEGMENTS - 1; k++)
      vp10_diff_update_prob(&r, &cm->fc->seg.tree_probs[k]);
  }

  for (j = 0; j < INTRA_MODES; j++)
    for (i = 0; i < INTRA_MODES - 1; ++i)
      vp10_diff_update_prob(&r, &fc->uv_mode_prob[j][i]);

#if CONFIG_EXT_PARTITION_TYPES
    for (i = 0; i < PARTITION_TYPES - 1; ++i)
      vp10_diff_update_prob(&r, &fc->partition_prob[0][i]);
    for (j = 1; j < PARTITION_CONTEXTS; ++j)
      for (i = 0; i < EXT_PARTITION_TYPES - 1; ++i)
        vp10_diff_update_prob(&r, &fc->partition_prob[j][i]);
#else
  for (j = 0; j < PARTITION_CONTEXTS; ++j)
    for (i = 0; i < PARTITION_TYPES - 1; ++i)
      vp10_diff_update_prob(&r, &fc->partition_prob[j][i]);
#endif  // CONFIG_EXT_PARTITION_TYPES

#if CONFIG_EXT_INTRA
  for (i = 0; i < INTRA_FILTERS + 1; ++i)
    for (j = 0; j < INTRA_FILTERS - 1; ++j)
      vp10_diff_update_prob(&r, &fc->intra_filter_probs[i][j]);
#endif  // CONFIG_EXT_INTRA

  if (frame_is_intra_only(cm)) {
    vp10_copy(cm->kf_y_prob, vp10_kf_y_mode_prob);
    for (k = 0; k < INTRA_MODES; k++)
      for (j = 0; j < INTRA_MODES; j++)
        for (i = 0; i < INTRA_MODES - 1; ++i)
          vp10_diff_update_prob(&r, &cm->kf_y_prob[k][j][i]);
  } else {
#if !CONFIG_REF_MV
    nmv_context *const nmvc = &fc->nmvc;
#endif

    read_inter_mode_probs(fc, &r);

#if CONFIG_EXT_INTER
    read_inter_compound_mode_probs(fc, &r);
    if (cm->reference_mode != COMPOUND_REFERENCE) {
      for (i = 0; i < BLOCK_SIZE_GROUPS; i++) {
        if (is_interintra_allowed_bsize_group(i)) {
          vp10_diff_update_prob(&r, &fc->interintra_prob[i]);
        }
      }
      for (i = 0; i < BLOCK_SIZE_GROUPS; i++) {
        for (j = 0; j < INTERINTRA_MODES - 1; j++)
          vp10_diff_update_prob(&r, &fc->interintra_mode_prob[i][j]);
      }
      for (i = 0; i < BLOCK_SIZES; i++) {
        if (is_interintra_allowed_bsize(i) && is_interintra_wedge_used(i)) {
          vp10_diff_update_prob(&r, &fc->wedge_interintra_prob[i]);
        }
      }
    }
    if (cm->reference_mode != SINGLE_REFERENCE) {
      for (i = 0; i < BLOCK_SIZES; i++) {
        if (is_interinter_wedge_used(i)) {
          vp10_diff_update_prob(&r, &fc->wedge_interinter_prob[i]);
        }
      }
    }
#endif  // CONFIG_EXT_INTER

#if CONFIG_OBMC || CONFIG_WARPED_MOTION
    for (i = BLOCK_8X8; i < BLOCK_SIZES; ++i) {
      for (j = 0; j < MOTION_VARIATIONS - 1; ++j)
        vp10_diff_update_prob(&r, &fc->motvar_prob[i][j]);
    }
#endif  // CONFIG_OBMC || CONFIG_WARPED_MOTION

    if (cm->interp_filter == SWITCHABLE)
      read_switchable_interp_probs(fc, &r);

    for (i = 0; i < INTRA_INTER_CONTEXTS; i++)
      vp10_diff_update_prob(&r, &fc->intra_inter_prob[i]);

    if (cm->reference_mode != SINGLE_REFERENCE)
      setup_compound_reference_mode(cm);

    read_frame_reference_mode_probs(cm, &r);

    for (j = 0; j < BLOCK_SIZE_GROUPS; j++)
      for (i = 0; i < INTRA_MODES - 1; ++i)
        vp10_diff_update_prob(&r, &fc->y_mode_prob[j][i]);

#if CONFIG_REF_MV
    for (i = 0; i < NMV_CONTEXTS; ++i)
      read_mv_probs(&fc->nmvc[i], cm->allow_high_precision_mv, &r);
#else
    read_mv_probs(nmvc, cm->allow_high_precision_mv, &r);
#endif
    read_ext_tx_probs(fc, &r);
#if CONFIG_SUPERTX
    if (!xd->lossless[0])
      read_supertx_probs(fc, &r);
#endif
  }

  return vp10_reader_has_error(&r);
}

#ifdef NDEBUG
#define debug_check_frame_counts(cm) (void)0
#else  // !NDEBUG
// Counts should only be incremented when frame_parallel_decoding_mode and
// error_resilient_mode are disabled.
static void debug_check_frame_counts(const VP10_COMMON *const cm) {
  FRAME_COUNTS zero_counts;
  vp10_zero(zero_counts);
  assert(cm->refresh_frame_context != REFRESH_FRAME_CONTEXT_BACKWARD ||
         cm->error_resilient_mode);
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
#if CONFIG_EXT_INTER
  assert(!memcmp(cm->counts.inter_compound_mode,
                 zero_counts.inter_compound_mode,
                 sizeof(cm->counts.inter_compound_mode)));
  assert(!memcmp(cm->counts.interintra, zero_counts.interintra,
                 sizeof(cm->counts.interintra)));
  assert(!memcmp(cm->counts.wedge_interintra, zero_counts.wedge_interintra,
                 sizeof(cm->counts.wedge_interintra)));
  assert(!memcmp(cm->counts.wedge_interinter, zero_counts.wedge_interinter,
                 sizeof(cm->counts.wedge_interinter)));
#endif  // CONFIG_EXT_INTER
#if CONFIG_OBMC || CONFIG_WARPED_MOTION
  assert(!memcmp(cm->counts.motvar, zero_counts.motvar,
                 sizeof(cm->counts.motvar)));
#endif  // CONFIG_OBMC || CONFIG_WARPED_MOTION
  assert(!memcmp(cm->counts.intra_inter, zero_counts.intra_inter,
                 sizeof(cm->counts.intra_inter)));
  assert(!memcmp(cm->counts.comp_inter, zero_counts.comp_inter,
                 sizeof(cm->counts.comp_inter)));
  assert(!memcmp(cm->counts.single_ref, zero_counts.single_ref,
                 sizeof(cm->counts.single_ref)));
  assert(!memcmp(cm->counts.comp_ref, zero_counts.comp_ref,
                 sizeof(cm->counts.comp_ref)));
#if CONFIG_EXT_REFS
  assert(!memcmp(cm->counts.comp_bwdref, zero_counts.comp_bwdref,
                 sizeof(cm->counts.comp_bwdref)));
#endif  // CONFIG_EXT_REFS
  assert(!memcmp(&cm->counts.tx_size, &zero_counts.tx_size,
                 sizeof(cm->counts.tx_size)));
  assert(!memcmp(cm->counts.skip, zero_counts.skip, sizeof(cm->counts.skip)));
#if CONFIG_REF_MV
  assert(!memcmp(&cm->counts.mv[0], &zero_counts.mv[0],
                 sizeof(cm->counts.mv[0])));
  assert(!memcmp(&cm->counts.mv[1], &zero_counts.mv[1],
                 sizeof(cm->counts.mv[0])));
#else
  assert(!memcmp(&cm->counts.mv, &zero_counts.mv, sizeof(cm->counts.mv)));
#endif
  assert(!memcmp(cm->counts.inter_ext_tx, zero_counts.inter_ext_tx,
                 sizeof(cm->counts.inter_ext_tx)));
  assert(!memcmp(cm->counts.intra_ext_tx, zero_counts.intra_ext_tx,
                 sizeof(cm->counts.intra_ext_tx)));
}
#endif  // NDEBUG

static struct vpx_read_bit_buffer *init_read_bit_buffer(
    VP10Decoder *pbi,
    struct vpx_read_bit_buffer *rb,
    const uint8_t *data,
    const uint8_t *data_end,
    uint8_t clear_data[MAX_VP9_HEADER_SIZE]) {
  rb->bit_offset = 0;
  rb->error_handler = error_handler;
  rb->error_handler_data = &pbi->common;
  if (pbi->decrypt_cb) {
    const int n = (int)VPXMIN(MAX_VP9_HEADER_SIZE, data_end - data);
    pbi->decrypt_cb(pbi->decrypt_state, data, clear_data, n);
    rb->bit_buffer = clear_data;
    rb->bit_buffer_end = clear_data + n;
  } else {
    rb->bit_buffer = data;
    rb->bit_buffer_end = data_end;
  }
  return rb;
}

//------------------------------------------------------------------------------

int vp10_read_sync_code(struct vpx_read_bit_buffer *const rb) {
  return vpx_rb_read_literal(rb, 8) == VP10_SYNC_CODE_0 &&
         vpx_rb_read_literal(rb, 8) == VP10_SYNC_CODE_1 &&
         vpx_rb_read_literal(rb, 8) == VP10_SYNC_CODE_2;
}

void vp10_read_frame_size(struct vpx_read_bit_buffer *rb,
                         int *width, int *height) {
  *width = vpx_rb_read_literal(rb, 16) + 1;
  *height = vpx_rb_read_literal(rb, 16) + 1;
}

BITSTREAM_PROFILE vp10_read_profile(struct vpx_read_bit_buffer *rb) {
  int profile = vpx_rb_read_bit(rb);
  profile |= vpx_rb_read_bit(rb) << 1;
  if (profile > 2)
    profile += vpx_rb_read_bit(rb);
  return (BITSTREAM_PROFILE) profile;
}

void vp10_decode_frame(VP10Decoder *pbi,
                      const uint8_t *data, const uint8_t *data_end,
                      const uint8_t **p_data_end) {
  VP10_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  struct vpx_read_bit_buffer rb;
  int context_updated = 0;
  uint8_t clear_data[MAX_VP9_HEADER_SIZE];
  const size_t first_partition_size = read_uncompressed_header(pbi,
      init_read_bit_buffer(pbi, &rb, data, data_end, clear_data));
  YV12_BUFFER_CONFIG *const new_fb = get_frame_new_buffer(cm);
  xd->cur_buf = new_fb;

  if (!first_partition_size) {
    // showing a frame directly
#if CONFIG_EXT_REFS
    if (cm->show_existing_frame)
      *p_data_end = data + vpx_rb_bytes_read(&rb);
    else
#endif  // CONFIG_EXT_REFS
      *p_data_end = data + (cm->profile <= PROFILE_2 ? 1 : 2);

    return;
  }

  data += vpx_rb_bytes_read(&rb);
  if (!read_is_valid(data, first_partition_size, data_end))
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt header length");

  cm->use_prev_frame_mvs = !cm->error_resilient_mode &&
                           cm->width == cm->last_width &&
                           cm->height == cm->last_height &&
                           !cm->last_intra_only &&
                           cm->last_show_frame &&
                           (cm->last_frame_type != KEY_FRAME);
#if CONFIG_EXT_REFS
  // NOTE(zoeliu): As cm->prev_frame can take neither a frame of
  //               show_exisiting_frame=1, nor can it take a frame not used as
  //               a reference, it is probable that by the time it is being
  //               referred to, the frame buffer it originally points to may
  //               already get expired and have been reassigned to the current
  //               newly coded frame. Hence, we need to check whether this is
  //               the case, and if yes, we have 2 choices:
  //               (1) Simply disable the use of previous frame mvs; or
  //               (2) Have cm->prev_frame point to one reference frame buffer,
  //                   e.g. LAST_FRAME.
  if (cm->use_prev_frame_mvs && !dec_is_ref_frame_buf(pbi, cm->prev_frame)) {
    // Reassign the LAST_FRAME buffer to cm->prev_frame.
    RefBuffer *last_fb_ref_buf = &cm->frame_refs[LAST_FRAME - LAST_FRAME];
    cm->prev_frame = &cm->buffer_pool->frame_bufs[last_fb_ref_buf->idx];
  }
#endif  // CONFIG_EXT_REFS

  vp10_setup_block_planes(xd, cm->subsampling_x, cm->subsampling_y);

  *cm->fc = cm->frame_contexts[cm->frame_context_idx];
  if (!cm->fc->initialized)
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Uninitialized entropy context.");

  vp10_zero(cm->counts);

  xd->corrupted = 0;
  new_fb->corrupted = read_compressed_header(pbi, data, first_partition_size);
  if (new_fb->corrupted)
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Decode failed. Frame data header is corrupted.");

  if (cm->lf.filter_level && !cm->skip_loop_filter) {
    vp10_loop_filter_frame_init(cm, cm->lf.filter_level);
  }

  // If encoded in frame parallel mode, frame context is ready after decoding
  // the frame header.
  if (cm->frame_parallel_decode &&
      cm->refresh_frame_context != REFRESH_FRAME_CONTEXT_BACKWARD) {
    VPxWorker *const worker = pbi->frame_worker_owner;
    FrameWorkerData *const frame_worker_data = worker->data1;
    if (cm->refresh_frame_context == REFRESH_FRAME_CONTEXT_FORWARD) {
      context_updated = 1;
      cm->frame_contexts[cm->frame_context_idx] = *cm->fc;
    }
    vp10_frameworker_lock_stats(worker);
    pbi->cur_buf->row = -1;
    pbi->cur_buf->col = -1;
    frame_worker_data->frame_context_ready = 1;
    // Signal the main thread that context is ready.
    vp10_frameworker_signal_stats(worker);
    vp10_frameworker_unlock_stats(worker);
  }

#if CONFIG_ENTROPY
  vp10_copy(cm->starting_coef_probs, cm->fc->coef_probs);
  cm->coef_probs_update_idx = 0;
#endif  // CONFIG_ENTROPY

  if (pbi->max_threads > 1
#if CONFIG_EXT_TILE
      && pbi->dec_tile_col < 0  // Decoding all columns
#endif  // CONFIG_EXT_TILE
      && cm->tile_cols > 1) {
    // Multi-threaded tile decoder
    *p_data_end = decode_tiles_mt(pbi, data + first_partition_size, data_end);
    if (!xd->corrupted) {
      if (!cm->skip_loop_filter) {
        // If multiple threads are used to decode tiles, then we use those
        // threads to do parallel loopfiltering.
        vp10_loop_filter_frame_mt(new_fb, cm, pbi->mb.plane,
                                 cm->lf.filter_level, 0, 0, pbi->tile_workers,
                                 pbi->num_tile_workers, &pbi->lf_row_sync);
      }
    } else {
      vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                         "Decode failed. Frame data is corrupted.");

    }
  } else {
    *p_data_end = decode_tiles(pbi, data + first_partition_size, data_end);
  }
#if CONFIG_LOOP_RESTORATION
  if (cm->rst_info.restoration_type != RESTORE_NONE) {
    vp10_loop_restoration_init(&cm->rst_internal,
                               &cm->rst_info,
                               cm->frame_type == KEY_FRAME);
    vp10_loop_restoration_rows(new_fb, cm, 0, cm->mi_rows, 0);
  }
#endif  // CONFIG_LOOP_RESTORATION

  if (!xd->corrupted) {
    if (cm->refresh_frame_context == REFRESH_FRAME_CONTEXT_BACKWARD) {
#if CONFIG_ENTROPY
      cm->partial_prob_update = 0;
#endif  // CONFIG_ENTROPY
      vp10_adapt_coef_probs(cm);
      vp10_adapt_intra_frame_probs(cm);

      if (!frame_is_intra_only(cm)) {
        vp10_adapt_inter_frame_probs(cm);
        vp10_adapt_mv_probs(cm, cm->allow_high_precision_mv);
      }
    } else {
      debug_check_frame_counts(cm);
    }
  } else {
    vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                       "Decode failed. Frame data is corrupted.");
  }

  // Non frame parallel update frame context here.
  if (!cm->error_resilient_mode && !context_updated)
    cm->frame_contexts[cm->frame_context_idx] = *cm->fc;
}
