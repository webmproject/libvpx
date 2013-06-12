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

#include "vp9/common/vp9_blockd.h"

void vp9_print_modes_and_motion_vectors(MODE_INFO *mi, int rows, int cols,
                                        int frame, char *file) {
  int mi_row;
  int mi_col;
  int mi_index = 0;
  FILE *mvs = fopen(file, "a");

  // Print out the macroblock Y modes
  fprintf(mvs, "SB Types for Frame %d\n", frame);

  for (mi_row = 0; mi_row < rows; mi_row++) {
    for (mi_col = 0; mi_col < cols; mi_col++) {
      fprintf(mvs, "%2d ", mi[mi_index].mbmi.sb_type);

      mi_index++;
    }

    fprintf(mvs, "\n");
    mi_index += 8;
  }

  // Print out the macroblock Y modes
  fprintf(mvs, "Mb Modes for Frame %d\n", frame);
  mi_index = 0;
  for (mi_row = 0; mi_row < rows; mi_row++) {
    for (mi_col = 0; mi_col < cols; mi_col++) {
      fprintf(mvs, "%2d ", mi[mi_index].mbmi.mode);

      mi_index++;
    }

    fprintf(mvs, "\n");
    mi_index += 8;
  }

  fprintf(mvs, "\n");

  mi_index = 0;
  fprintf(mvs, "Mb mv ref for Frame %d\n", frame);

  for (mi_row = 0; mi_row < rows; mi_row++) {
    for (mi_col = 0; mi_col < cols; mi_col++) {
      fprintf(mvs, "%2d ", mi[mi_index].mbmi.ref_frame[0]);

      mi_index++;
    }

    fprintf(mvs, "\n");
    mi_index += 8;
  }
  fprintf(mvs, "\n");

  mi_index = 0;
  fprintf(mvs, "Mb mv ref for Frame %d\n", frame);

  for (mi_row = 0; mi_row < rows; mi_row++) {
    for (mi_col = 0; mi_col < cols; mi_col++) {
      fprintf(mvs, "%4d:%4d ", mi[mi_index].mbmi.mv[0].as_mv.row,
              mi[mi_index].mbmi.mv[0].as_mv.col);

      mi_index++;
    }

    fprintf(mvs, "\n");
    mi_index += 8;
  }

  fprintf(mvs, "\n");

  /* print out the macroblock txform sizes */
  mi_index = 0;
  fprintf(mvs, "TXFM size for Frame %d\n", frame);

  for (mi_row = 0; mi_row < rows; mi_row++) {
    for (mi_col = 0; mi_col < cols; mi_col++) {
      fprintf(mvs, "%2d ", mi[mi_index].mbmi.txfm_size);

      mi_index++;
    }

    mi_index += 8;
    fprintf(mvs, "\n");
  }

  fprintf(mvs, "\n");

  /* print out the macroblock UV modes */
  mi_index = 0;
  fprintf(mvs, "UV Modes for Frame %d\n", frame);

  for (mi_row = 0; mi_row < rows; mi_row++) {
    for (mi_col = 0; mi_col < cols; mi_col++) {
      fprintf(mvs, "%2d ", mi[mi_index].mbmi.uv_mode);

      mi_index++;
    }

    mi_index += 8;
    fprintf(mvs, "\n");
  }

  fprintf(mvs, "\n");

  /* print out the macroblock mvs */
  mi_index = 0;
  fprintf(mvs, "MVs for Frame %d\n", frame);

  for (mi_row = 0; mi_row < rows; mi_row++) {
    for (mi_col = 0; mi_col < cols; mi_col++) {
      fprintf(mvs, "%5d:%-5d", mi[mi_index].mbmi.mv[0].as_mv.row / 2,
              mi[mi_index].mbmi.mv[0].as_mv.col / 2);

      mi_index++;
    }

    mi_index += 8;
    fprintf(mvs, "\n");
  }

  fprintf(mvs, "\n");

  fclose(mvs);
}
