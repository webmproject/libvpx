/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <string.h>
#include <math.h>

#include "./aom_scale_rtcd.h"
#include "aom/aom_integer.h"
#include "av1/common/dering.h"
#include "av1/common/onyxc_int.h"
#include "av1/common/reconinter.h"
#include "av1/common/od_dering.h"

int compute_level_from_index(int global_level, int gi) {
  static const int dering_gains[DERING_REFINEMENT_LEVELS] = { 0, 11, 16, 22 };
  int level;
  if (global_level == 0) return 0;
  level = (global_level * dering_gains[gi] + 8) >> 4;
  return clamp(level, gi, MAX_DERING_LEVEL - 1);
}

int sb_all_skip(const AV1_COMMON *const cm, int mi_row, int mi_col) {
  int r, c;
  int maxc, maxr;
  int skip = 1;
  maxc = cm->mi_cols - mi_col;
  maxr = cm->mi_rows - mi_row;
  if (maxr > MAX_MIB_SIZE) maxr = MAX_MIB_SIZE;
  if (maxc > MAX_MIB_SIZE) maxc = MAX_MIB_SIZE;
  for (r = 0; r < maxr; r++) {
    for (c = 0; c < maxc; c++) {
      skip = skip &&
             cm->mi_grid_visible[(mi_row + r) * cm->mi_stride + mi_col + c]
                 ->mbmi.skip;
    }
  }
  return skip;
}

void av1_dering_frame(YV12_BUFFER_CONFIG *frame, AV1_COMMON *cm,
                      MACROBLOCKD *xd, int global_level) {
  int r, c;
  int sbr, sbc;
  int nhsb, nvsb;
  od_dering_in *src[3];
  unsigned char *bskip;
  int dir[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS] = { { 0 } };
  int stride;
  int bsize[3];
  int dec[3];
  int pli;
  int coeff_shift = AOMMAX(cm->bit_depth - 8, 0);
  nvsb = (cm->mi_rows + MAX_MIB_SIZE - 1) / MAX_MIB_SIZE;
  nhsb = (cm->mi_cols + MAX_MIB_SIZE - 1) / MAX_MIB_SIZE;
  bskip = aom_malloc(sizeof(*bskip) * cm->mi_rows * cm->mi_cols);
  av1_setup_dst_planes(xd->plane, frame, 0, 0);
  for (pli = 0; pli < 3; pli++) {
    dec[pli] = xd->plane[pli].subsampling_x;
    bsize[pli] = 8 >> dec[pli];
  }
  stride = bsize[0] * cm->mi_cols;
  for (pli = 0; pli < 3; pli++) {
    src[pli] = aom_malloc(sizeof(*src) * cm->mi_rows * cm->mi_cols * 64);
    for (r = 0; r < bsize[pli] * cm->mi_rows; ++r) {
      for (c = 0; c < bsize[pli] * cm->mi_cols; ++c) {
#if CONFIG_AOM_HIGHBITDEPTH
        if (cm->use_highbitdepth) {
          src[pli][r * stride + c] = CONVERT_TO_SHORTPTR(
              xd->plane[pli].dst.buf)[r * xd->plane[pli].dst.stride + c];
        } else {
#endif
          src[pli][r * stride + c] =
              xd->plane[pli].dst.buf[r * xd->plane[pli].dst.stride + c];
#if CONFIG_AOM_HIGHBITDEPTH
        }
#endif
      }
    }
  }
  for (r = 0; r < cm->mi_rows; ++r) {
    for (c = 0; c < cm->mi_cols; ++c) {
      const MB_MODE_INFO *mbmi =
          &cm->mi_grid_visible[r * cm->mi_stride + c]->mbmi;
      bskip[r * cm->mi_cols + c] = mbmi->skip;
    }
  }
  for (sbr = 0; sbr < nvsb; sbr++) {
    for (sbc = 0; sbc < nhsb; sbc++) {
      int level;
      int nhb, nvb;
      nhb = AOMMIN(MAX_MIB_SIZE, cm->mi_cols - MAX_MIB_SIZE * sbc);
      nvb = AOMMIN(MAX_MIB_SIZE, cm->mi_rows - MAX_MIB_SIZE * sbr);
      for (pli = 0; pli < 3; pli++) {
        int16_t dst[MAX_MIB_SIZE * MAX_MIB_SIZE * 8 * 8];
        int threshold;
        level = compute_level_from_index(
            global_level,
            cm->mi_grid_visible[MAX_MIB_SIZE * sbr * cm->mi_stride +
                                MAX_MIB_SIZE * sbc]
                ->mbmi.dering_gain);
        /* FIXME: This is a temporary hack that uses more conservative
           deringing for chroma. */
        if (pli) level = (level * 5 + 4) >> 3;
        if (sb_all_skip(cm, sbr * MAX_MIB_SIZE, sbc * MAX_MIB_SIZE)) continue;
        threshold = level << coeff_shift;
        od_dering(dst, MAX_MIB_SIZE * bsize[pli],
                  &src[pli][sbr * stride * bsize[pli] * MAX_MIB_SIZE +
                            sbc * bsize[pli] * MAX_MIB_SIZE],
                  stride, nhb, nvb, sbc, sbr, nhsb, nvsb, dec[pli], dir, pli,
                  &bskip[MAX_MIB_SIZE * sbr * cm->mi_cols + MAX_MIB_SIZE * sbc],
                  cm->mi_cols, threshold, coeff_shift);
        for (r = 0; r < bsize[pli] * nvb; ++r) {
          for (c = 0; c < bsize[pli] * nhb; ++c) {
#if CONFIG_AOM_HIGHBITDEPTH
            if (cm->use_highbitdepth) {
              CONVERT_TO_SHORTPTR(xd->plane[pli].dst.buf)
              [xd->plane[pli].dst.stride *
                   (bsize[pli] * MAX_MIB_SIZE * sbr + r) +
               sbc * bsize[pli] * MAX_MIB_SIZE + c] =
                  dst[r * MAX_MIB_SIZE * bsize[pli] + c];
            } else {
#endif
              xd->plane[pli].dst.buf[xd->plane[pli].dst.stride *
                                         (bsize[pli] * MAX_MIB_SIZE * sbr + r) +
                                     sbc * bsize[pli] * MAX_MIB_SIZE + c] =
                  dst[r * MAX_MIB_SIZE * bsize[pli] + c];
#if CONFIG_AOM_HIGHBITDEPTH
            }
#endif
          }
        }
      }
    }
  }
  for (pli = 0; pli < 3; pli++) {
    aom_free(src[pli]);
  }
  aom_free(bskip);
}
