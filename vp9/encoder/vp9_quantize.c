/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>
#include "vpx_mem/vpx_mem.h"

#include "vp9/encoder/vp9_onyx_int.h"
#include "vp9/encoder/vp9_quantize.h"
#include "vp9/common/vp9_quant_common.h"

#include "vp9/common/vp9_seg_common.h"

#ifdef ENC_DEBUG
extern int enc_debug;
#endif

static INLINE int plane_idx(int plane) {
  return plane == 0 ? 0 :
         plane == 1 ? 16 : 20;
}

static void quantize(int16_t *zbin_boost_orig_ptr,
                     int16_t *coeff_ptr, int n_coeffs, int skip_block,
                     int16_t *zbin_ptr, int16_t *round_ptr, int16_t *quant_ptr,
                     uint8_t *quant_shift_ptr,
                     int16_t *qcoeff_ptr, int16_t *dqcoeff_ptr,
                     int16_t *dequant_ptr, int zbin_oq_value,
                     uint16_t *eob_ptr,
                     const int *scan, int mul) {
  int i, rc, eob;
  int zbin;
  int x, y, z, sz;
  int zero_run = 0;
  int16_t *zbin_boost_ptr = zbin_boost_orig_ptr;

  vpx_memset(qcoeff_ptr, 0, n_coeffs*sizeof(int16_t));
  vpx_memset(dqcoeff_ptr, 0, n_coeffs*sizeof(int16_t));

  eob = -1;

  if (!skip_block) {
    for (i = 0; i < n_coeffs; i++) {
      rc   = scan[i];
      z    = coeff_ptr[rc] * mul;

      zbin = (zbin_ptr[rc != 0] + zbin_boost_ptr[zero_run] + zbin_oq_value);
      zero_run += (zero_run < 15);

      sz = (z >> 31);                               // sign of z
      x  = (z ^ sz) - sz;                           // x = abs(z)

      if (x >= zbin) {
        x += (round_ptr[rc != 0]);
        y  = ((int)(((int)(x * quant_ptr[rc != 0]) >> 16) + x))
            >> quant_shift_ptr[rc != 0];            // quantize (x)
        x  = (y ^ sz) - sz;                         // get the sign back
        qcoeff_ptr[rc]  = x;                        // write to destination
        dqcoeff_ptr[rc] = x * dequant_ptr[rc != 0] / mul;  // dequantized value

        if (y) {
          eob = i;                                  // last nonzero coeffs
          zero_run = 0;
        }
      }
    }
  }

  *eob_ptr = eob + 1;
}

void vp9_regular_quantize_b_4x4(MACROBLOCK *mb, int b_idx, TX_TYPE tx_type,
                                int y_blocks) {
  MACROBLOCKD *const xd = &mb->e_mbd;
  const struct plane_block_idx pb_idx = plane_block_idx(y_blocks, b_idx);
  const int *pt_scan = get_scan_4x4(tx_type);

  quantize(mb->plane[pb_idx.plane].zrun_zbin_boost,
           BLOCK_OFFSET(mb->plane[pb_idx.plane].coeff, pb_idx.block, 16),
           16, mb->skip_block,
           mb->plane[pb_idx.plane].zbin,
           mb->plane[pb_idx.plane].round,
           mb->plane[pb_idx.plane].quant,
           mb->plane[pb_idx.plane].quant_shift,
           BLOCK_OFFSET(xd->plane[pb_idx.plane].qcoeff, pb_idx.block, 16),
           BLOCK_OFFSET(xd->plane[pb_idx.plane].dqcoeff, pb_idx.block, 16),
           xd->plane[pb_idx.plane].dequant,
           mb->plane[pb_idx.plane].zbin_extra,
           &xd->plane[pb_idx.plane].eobs[pb_idx.block],
           pt_scan, 1);
}

void vp9_regular_quantize_b_8x8(MACROBLOCK *mb, int b_idx, TX_TYPE tx_type,
                                int y_blocks) {
  MACROBLOCKD *const xd = &mb->e_mbd;
  const struct plane_block_idx pb_idx = plane_block_idx(y_blocks, b_idx);
  const int *pt_scan = get_scan_8x8(tx_type);

  quantize(mb->plane[pb_idx.plane].zrun_zbin_boost,
           BLOCK_OFFSET(mb->plane[pb_idx.plane].coeff, pb_idx.block, 16),
           64, mb->skip_block,
           mb->plane[pb_idx.plane].zbin,
           mb->plane[pb_idx.plane].round,
           mb->plane[pb_idx.plane].quant,
           mb->plane[pb_idx.plane].quant_shift,
           BLOCK_OFFSET(xd->plane[pb_idx.plane].qcoeff, pb_idx.block, 16),
           BLOCK_OFFSET(xd->plane[pb_idx.plane].dqcoeff, pb_idx.block, 16),
           xd->plane[pb_idx.plane].dequant,
           mb->plane[pb_idx.plane].zbin_extra,
           &xd->plane[pb_idx.plane].eobs[pb_idx.block],
           pt_scan, 1);
}

void vp9_regular_quantize_b_16x16(MACROBLOCK *mb, int b_idx, TX_TYPE tx_type,
                                  int y_blocks) {
  MACROBLOCKD *const xd = &mb->e_mbd;
  const struct plane_block_idx pb_idx = plane_block_idx(y_blocks, b_idx);
  const int *pt_scan = get_scan_16x16(tx_type);

  quantize(mb->plane[pb_idx.plane].zrun_zbin_boost,
           BLOCK_OFFSET(mb->plane[pb_idx.plane].coeff, pb_idx.block, 16),
           256, mb->skip_block,
           mb->plane[pb_idx.plane].zbin,
           mb->plane[pb_idx.plane].round,
           mb->plane[pb_idx.plane].quant,
           mb->plane[pb_idx.plane].quant_shift,
           BLOCK_OFFSET(xd->plane[pb_idx.plane].qcoeff, pb_idx.block, 16),
           BLOCK_OFFSET(xd->plane[pb_idx.plane].dqcoeff, pb_idx.block, 16),
           xd->plane[pb_idx.plane].dequant,
           mb->plane[pb_idx.plane].zbin_extra,
           &xd->plane[pb_idx.plane].eobs[pb_idx.block],
           pt_scan, 1);
}

void vp9_regular_quantize_b_32x32(MACROBLOCK *mb, int b_idx, int y_blocks) {
  MACROBLOCKD *const xd = &mb->e_mbd;
  const struct plane_block_idx pb_idx = plane_block_idx(y_blocks, b_idx);

  quantize(mb->plane[pb_idx.plane].zrun_zbin_boost,
           BLOCK_OFFSET(mb->plane[pb_idx.plane].coeff, pb_idx.block, 16),
           1024, mb->skip_block,
           mb->plane[pb_idx.plane].zbin,
           mb->plane[pb_idx.plane].round,
           mb->plane[pb_idx.plane].quant,
           mb->plane[pb_idx.plane].quant_shift,
           BLOCK_OFFSET(xd->plane[pb_idx.plane].qcoeff, pb_idx.block, 16),
           BLOCK_OFFSET(xd->plane[pb_idx.plane].dqcoeff, pb_idx.block, 16),
           xd->plane[pb_idx.plane].dequant,
           mb->plane[pb_idx.plane].zbin_extra,
           &xd->plane[pb_idx.plane].eobs[pb_idx.block],
           vp9_default_zig_zag1d_32x32, 2);
}

void vp9_quantize_sby_32x32(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bw = 1 << (b_width_log2(bsize) - 3);
  const int bh = 1 << (b_height_log2(bsize) - 3);
  int n;

  for (n = 0; n < bw * bh; n++)
    vp9_regular_quantize_b_32x32(x, n * 64, bw * bh * 64);
}

void vp9_quantize_sby_16x16(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = b_width_log2(bsize) - 2, bw = 1 << bwl;
  const int bh = 1 << (b_height_log2(bsize) - 2);
  const int bstride = 16 << bwl;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;
    TX_TYPE tx_type = get_tx_type_16x16(&x->e_mbd,
                                        4 * x_idx + y_idx * bstride);
    x->quantize_b_16x16(x, n * 16, tx_type, 16 * bw * bh);
  }
}

void vp9_quantize_sby_8x8(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = b_width_log2(bsize) - 1, bw = 1 << bwl;
  const int bh = 1 << (b_height_log2(bsize) - 1);
  const int bstride = 4 << bwl;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;
    TX_TYPE tx_type = get_tx_type_8x8(&x->e_mbd,
                                      2 * x_idx + y_idx * bstride);
    x->quantize_b_8x8(x, n * 4, tx_type, 4 * bw * bh);
  }
}

void vp9_quantize_sby_4x4(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = b_width_log2(bsize), bw = 1 << bwl;
  const int bh = 1 << b_height_log2(bsize);
  MACROBLOCKD *const xd = &x->e_mbd;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const TX_TYPE tx_type = get_tx_type_4x4(xd, n);
    x->quantize_b_4x4(x, n, tx_type, bw * bh);
  }
}

void vp9_quantize_sbuv_32x32(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  assert(bsize == BLOCK_SIZE_SB64X64);
  vp9_regular_quantize_b_32x32(x, 256, 256);
  vp9_regular_quantize_b_32x32(x, 320, 256);
}

void vp9_quantize_sbuv_16x16(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = b_width_log2(bsize) - 2;
  const int bhl = b_height_log2(bsize) - 2;
  const int uoff = 16 << (bhl + bwl);
  int i;

  for (i = uoff; i < ((uoff * 3) >> 1); i += 16)
    x->quantize_b_16x16(x, i, DCT_DCT, uoff);
}

void vp9_quantize_sbuv_8x8(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = b_width_log2(bsize) - 2;
  const int bhl = b_height_log2(bsize) - 2;
  const int uoff = 16 << (bhl + bwl);
  int i;

  for (i = uoff; i < ((uoff * 3) >> 1); i += 4)
    x->quantize_b_8x8(x, i, DCT_DCT, uoff);
}

void vp9_quantize_sbuv_4x4(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = b_width_log2(bsize) - 2;
  const int bhl = b_height_log2(bsize) - 2;
  const int uoff = 16 << (bhl + bwl);
  int i;

  for (i = uoff; i < ((uoff * 3) >> 1); i++)
    x->quantize_b_4x4(x, i, DCT_DCT, uoff);
}

/* quantize_b_pair function pointer in MACROBLOCK structure is set to one of
 * these two C functions if corresponding optimized routine is not available.
 * NEON optimized version implements currently the fast quantization for pair
 * of blocks. */
void vp9_regular_quantize_b_4x4_pair(MACROBLOCK *x, int b_idx1, int b_idx2,
                                     int y_blocks) {
  vp9_regular_quantize_b_4x4(x, b_idx1, DCT_DCT, y_blocks);
  vp9_regular_quantize_b_4x4(x, b_idx2, DCT_DCT, y_blocks);
}

static void invert_quant(int16_t *quant, uint8_t *shift, int d) {
  unsigned t;
  int l;
  t = d;
  for (l = 0; t > 1; l++)
    t >>= 1;
  t = 1 + (1 << (16 + l)) / d;
  *quant = (int16_t)(t - (1 << 16));
  *shift = l;
}

void vp9_init_quantizer(VP9_COMP *cpi) {
  int i;
  int quant_val;
  int q;

  static const int zbin_boost[16] = { 0,  0,  0,  8,  8,  8, 10, 12,
                                     14, 16, 20, 24, 28, 32, 36, 40 };

  for (q = 0; q < QINDEX_RANGE; q++) {
    int qzbin_factor = (vp9_dc_quant(q, 0) < 148) ? 84 : 80;
    int qrounding_factor = 48;
    if (q == 0) {
      qzbin_factor = 64;
      qrounding_factor = 64;
    }
    // dc values
    quant_val = vp9_dc_quant(q, cpi->common.y_dc_delta_q);
    invert_quant(cpi->Y1quant[q] + 0, cpi->Y1quant_shift[q] + 0, quant_val);
    cpi->Y1zbin[q][0] = ROUND_POWER_OF_TWO(qzbin_factor * quant_val, 7);
    cpi->Y1round[q][0] = (qrounding_factor * quant_val) >> 7;
    cpi->common.y_dequant[q][0] = quant_val;
    cpi->zrun_zbin_boost_y1[q][0] = (quant_val * zbin_boost[0]) >> 7;

    quant_val = vp9_dc_uv_quant(q, cpi->common.uv_dc_delta_q);
    invert_quant(cpi->UVquant[q] + 0, cpi->UVquant_shift[q] + 0, quant_val);
    cpi->UVzbin[q][0] = ROUND_POWER_OF_TWO(qzbin_factor * quant_val, 7);
    cpi->UVround[q][0] = (qrounding_factor * quant_val) >> 7;
    cpi->common.uv_dequant[q][0] = quant_val;
    cpi->zrun_zbin_boost_uv[q][0] = (quant_val * zbin_boost[0]) >> 7;

    // all the 4x4 ac values =;
    for (i = 1; i < 16; i++) {
      int rc = vp9_default_zig_zag1d_4x4[i];

      quant_val = vp9_ac_yquant(q);
      invert_quant(cpi->Y1quant[q] + rc, cpi->Y1quant_shift[q] + rc, quant_val);
      cpi->Y1zbin[q][rc] = ROUND_POWER_OF_TWO(qzbin_factor * quant_val, 7);
      cpi->Y1round[q][rc] = (qrounding_factor * quant_val) >> 7;
      cpi->common.y_dequant[q][rc] = quant_val;
      cpi->zrun_zbin_boost_y1[q][i] =
          ROUND_POWER_OF_TWO(quant_val * zbin_boost[i], 7);

      quant_val = vp9_ac_uv_quant(q, cpi->common.uv_ac_delta_q);
      invert_quant(cpi->UVquant[q] + rc, cpi->UVquant_shift[q] + rc, quant_val);
      cpi->UVzbin[q][rc] = ROUND_POWER_OF_TWO(qzbin_factor * quant_val, 7);
      cpi->UVround[q][rc] = (qrounding_factor * quant_val) >> 7;
      cpi->common.uv_dequant[q][rc] = quant_val;
      cpi->zrun_zbin_boost_uv[q][i] =
          ROUND_POWER_OF_TWO(quant_val * zbin_boost[i], 7);
    }
  }
}

void vp9_mb_init_quantizer(VP9_COMP *cpi, MACROBLOCK *x) {
  int i;
  int qindex;
  MACROBLOCKD *xd = &x->e_mbd;
  int zbin_extra;
  int segment_id = xd->mode_info_context->mbmi.segment_id;

  // Select the baseline MB Q index allowing for any segment level change.
  if (vp9_segfeature_active(xd, segment_id, SEG_LVL_ALT_Q)) {
    if (xd->mb_segment_abs_delta == SEGMENT_ABSDATA) {
      // Abs Value
      qindex = vp9_get_segdata(xd, segment_id, SEG_LVL_ALT_Q);
    } else {
      // Delta Value
      qindex = cpi->common.base_qindex +
                 vp9_get_segdata(xd, segment_id, SEG_LVL_ALT_Q);

      // Clamp to valid range
      qindex = clamp(qindex, 0, MAXQ);
    }
  } else {
    qindex = cpi->common.base_qindex;
  }

  // Y
  zbin_extra = (cpi->common.y_dequant[qindex][1] *
                 (cpi->zbin_mode_boost + x->act_zbin_adj)) >> 7;

  x->plane[0].quant = cpi->Y1quant[qindex];
  x->plane[0].quant_shift = cpi->Y1quant_shift[qindex];
  x->plane[0].zbin = cpi->Y1zbin[qindex];
  x->plane[0].round = cpi->Y1round[qindex];
  x->plane[0].zrun_zbin_boost = cpi->zrun_zbin_boost_y1[qindex];
  x->plane[0].zbin_extra = (int16_t)zbin_extra;
  x->e_mbd.plane[0].dequant = cpi->common.y_dequant[qindex];

  // UV
  zbin_extra = (cpi->common.uv_dequant[qindex][1] *
                (cpi->zbin_mode_boost + x->act_zbin_adj)) >> 7;

  for (i = 1; i < 3; i++) {
    x->plane[i].quant = cpi->UVquant[qindex];
    x->plane[i].quant_shift = cpi->UVquant_shift[qindex];
    x->plane[i].zbin = cpi->UVzbin[qindex];
    x->plane[i].round = cpi->UVround[qindex];
    x->plane[i].zrun_zbin_boost = cpi->zrun_zbin_boost_uv[qindex];
    x->plane[i].zbin_extra = (int16_t)zbin_extra;
    x->e_mbd.plane[i].dequant = cpi->common.uv_dequant[qindex];
  }

  x->skip_block = vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP);

  /* save this macroblock QIndex for vp9_update_zbin_extra() */
  x->e_mbd.q_index = qindex;
}

void vp9_update_zbin_extra(VP9_COMP *cpi, MACROBLOCK *x) {
  const int qindex = x->e_mbd.q_index;
  const int y_zbin_extra = (cpi->common.y_dequant[qindex][1] *
                (cpi->zbin_mode_boost + x->act_zbin_adj)) >> 7;
  const int uv_zbin_extra = (cpi->common.uv_dequant[qindex][1] *
                  (cpi->zbin_mode_boost + x->act_zbin_adj)) >> 7;

  x->plane[0].zbin_extra = (int16_t)y_zbin_extra;
  x->plane[1].zbin_extra = (int16_t)uv_zbin_extra;
  x->plane[2].zbin_extra = (int16_t)uv_zbin_extra;
}

void vp9_frame_init_quantizer(VP9_COMP *cpi) {
  // Clear Zbin mode boost for default case
  cpi->zbin_mode_boost = 0;

  // MB level quantizer setup
  vp9_mb_init_quantizer(cpi, &cpi->mb);
}

void vp9_set_quantizer(struct VP9_COMP *cpi, int Q) {
  VP9_COMMON *cm = &cpi->common;

  cm->base_qindex = Q;

  // Set lossless mode
  if (cm->base_qindex <= 4)
    cm->base_qindex = 0;

  // if any of the delta_q values are changing update flag will
  // have to be set.
  cm->y_dc_delta_q = 0;
  cm->uv_dc_delta_q = 0;
  cm->uv_ac_delta_q = 0;

  // quantizer has to be reinitialized if any delta_q changes.
  // As there are not any here for now this is inactive code.
  // if(update)
  //    vp9_init_quantizer(cpi);
}
