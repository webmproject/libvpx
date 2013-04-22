/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "./vpx_config.h"
#include "vpx_mem/vpx_mem.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9_rtcd.h"

#if CONFIG_NEWBINTRAMODES
static int find_grad_measure(uint8_t *x, int stride, int n, int tx, int ty,
                             int dx, int dy) {
  int i, j;
  int count = 0, gsum = 0, gdiv;
  /* TODO: Make this code more efficient by breaking up into two loops */
  for (i = -ty; i < n; ++i)
    for (j = -tx; j < n; ++j) {
      int g;
      if (i >= 0 && j >= 0) continue;
      if (i + dy >= 0 && j + dx >= 0) continue;
      if (i + dy < -ty || i + dy >= n || j + dx < -tx || j + dx >= n) continue;
      g = abs(x[(i + dy) * stride + j + dx] - x[i * stride + j]);
      gsum += g * g;
      count++;
    }
  gdiv = (dx * dx + dy * dy) * count;
  return ((gsum << 8) + (gdiv >> 1)) / gdiv;
}

#if CONTEXT_PRED_REPLACEMENTS == 6
B_PREDICTION_MODE vp9_find_dominant_direction(uint8_t *ptr,
                                              int stride, int n,
                                              int tx, int ty) {
  int g[8], i, imin, imax;
  g[1] = find_grad_measure(ptr, stride, n, tx, ty,  2, 1);
  g[2] = find_grad_measure(ptr, stride, n, tx, ty,  1, 1);
  g[3] = find_grad_measure(ptr, stride, n, tx, ty,  1, 2);
  g[5] = find_grad_measure(ptr, stride, n, tx, ty, -1, 2);
  g[6] = find_grad_measure(ptr, stride, n, tx, ty, -1, 1);
  g[7] = find_grad_measure(ptr, stride, n, tx, ty, -2, 1);
  imin = 1;
  for (i = 2; i < 8; i += 1 + (i == 3))
    imin = (g[i] < g[imin] ? i : imin);
  imax = 1;
  for (i = 2; i < 8; i += 1 + (i == 3))
    imax = (g[i] > g[imax] ? i : imax);
  /*
  printf("%d %d %d %d %d %d = %d %d\n",
         g[1], g[2], g[3], g[5], g[6], g[7], imin, imax);
         */
  switch (imin) {
    case 1:
      return B_D153_PRED;
    case 2:
      return B_D135_PRED;
    case 3:
      return B_D117_PRED;
    case 5:
      return B_D63_PRED;
    case 6:
      return B_D45_PRED;
    case 7:
      return B_D27_PRED;
    default:
      assert(0);
  }
}
#elif CONTEXT_PRED_REPLACEMENTS == 4
B_PREDICTION_MODE vp9_find_dominant_direction(uint8_t *ptr,
                                              int stride, int n,
                                              int tx, int ty) {
  int g[8], i, imin, imax;
  g[1] = find_grad_measure(ptr, stride, n, tx, ty,  2, 1);
  g[3] = find_grad_measure(ptr, stride, n, tx, ty,  1, 2);
  g[5] = find_grad_measure(ptr, stride, n, tx, ty, -1, 2);
  g[7] = find_grad_measure(ptr, stride, n, tx, ty, -2, 1);
  imin = 1;
  for (i = 3; i < 8; i+=2)
    imin = (g[i] < g[imin] ? i : imin);
  imax = 1;
  for (i = 3; i < 8; i+=2)
    imax = (g[i] > g[imax] ? i : imax);
  /*
  printf("%d %d %d %d = %d %d\n",
         g[1], g[3], g[5], g[7], imin, imax);
         */
  switch (imin) {
    case 1:
      return B_D153_PRED;
    case 3:
      return B_D117_PRED;
    case 5:
      return B_D63_PRED;
    case 7:
      return B_D27_PRED;
    default:
      assert(0);
  }
}
#elif CONTEXT_PRED_REPLACEMENTS == 0
B_PREDICTION_MODE vp9_find_dominant_direction(uint8_t *ptr,
                                              int stride, int n,
                                              int tx, int ty) {
  int g[8], i, imin, imax;
  g[0] = find_grad_measure(ptr, stride, n, tx, ty,  1, 0);
  g[1] = find_grad_measure(ptr, stride, n, tx, ty,  2, 1);
  g[2] = find_grad_measure(ptr, stride, n, tx, ty,  1, 1);
  g[3] = find_grad_measure(ptr, stride, n, tx, ty,  1, 2);
  g[4] = find_grad_measure(ptr, stride, n, tx, ty,  0, 1);
  g[5] = find_grad_measure(ptr, stride, n, tx, ty, -1, 2);
  g[6] = find_grad_measure(ptr, stride, n, tx, ty, -1, 1);
  g[7] = find_grad_measure(ptr, stride, n, tx, ty, -2, 1);
  imax = 0;
  for (i = 1; i < 8; i++)
    imax = (g[i] > g[imax] ? i : imax);
  imin = 0;
  for (i = 1; i < 8; i++)
    imin = (g[i] < g[imin] ? i : imin);

  switch (imin) {
    case 0:
      return B_H_PRED;
    case 1:
      return B_D153_PRED;
    case 2:
      return B_D135_PRED;
    case 3:
      return B_D117_PRED;
    case 4:
      return B_V_PRED;
    case 5:
      return B_D63_PRED;
    case 6:
      return B_D45_PRED;
    case 7:
      return B_D27_PRED;
    default:
      assert(0);
  }
}
#endif

B_PREDICTION_MODE vp9_find_bpred_context(MACROBLOCKD *xd, BLOCKD *x) {
  const int block_idx = x - xd->block;
  const int have_top = (block_idx >> 2) || xd->up_available;
  const int have_left = (block_idx & 3)  || xd->left_available;
  uint8_t *ptr = *(x->base_dst) + x->dst;
  int stride = x->dst_stride;
  int tx = have_left ? 4 : 0;
  int ty = have_top ? 4 : 0;
  if (!have_left && !have_top)
    return B_DC_PRED;
  return vp9_find_dominant_direction(ptr, stride, 4, tx, ty);
}
#endif

void vp9_intra4x4_predict(MACROBLOCKD *xd,
                          BLOCKD *x,
                          int b_mode,
                          uint8_t *predictor,
                          int ps) {
  int i, r, c;
  const int block_idx = x - xd->block;
  const int have_top = (block_idx >> 2) || xd->up_available;
  const int have_left = (block_idx & 3)  || xd->left_available;
  const int have_right = (block_idx & 3) != 3 || xd->right_available;
  uint8_t left[4], above[8], top_left;
  /*
   * 127 127 127 .. 127 127 127 127 127 127
   * 129  A   B  ..  Y   Z
   * 129  C   D  ..  W   X
   * 129  E   F  ..  U   V
   * 129  G   H  ..  S   T   T   T   T   T
   *  ..
   */

  if (have_left) {
    uint8_t *left_ptr = *(x->base_dst) + x->dst - 1;
    const int stride = x->dst_stride;

    left[0] = left_ptr[0 * stride];
    left[1] = left_ptr[1 * stride];
    left[2] = left_ptr[2 * stride];
    left[3] = left_ptr[3 * stride];
  } else {
    left[0] = left[1] = left[2] = left[3] = 129;
  }

  if (have_top) {
    uint8_t *above_ptr = *(x->base_dst) + x->dst - x->dst_stride;
    top_left = have_left ? above_ptr[-1] : 127;

    above[0] = above_ptr[0];
    above[1] = above_ptr[1];
    above[2] = above_ptr[2];
    above[3] = above_ptr[3];
    if (((block_idx & 3) != 3) ||
        (have_right && block_idx == 3 &&
         ((xd->mb_index != 3 && xd->sb_index != 3) ||
          ((xd->mb_index & 1) == 0 && xd->sb_index == 3)))) {
      above[4] = above_ptr[4];
      above[5] = above_ptr[5];
      above[6] = above_ptr[6];
      above[7] = above_ptr[7];
    } else if (have_right) {
      uint8_t *above_right = above_ptr + 4;

      if (xd->sb_index == 3 && (xd->mb_index & 1))
        above_right -= 32 * x->dst_stride;
      if (xd->mb_index == 3)
        above_right -= 16 * x->dst_stride;
      above_right -= (block_idx & ~3) * x->dst_stride;

      /* use a more distant above-right (from closest available top-right
       * corner), but with a "localized DC" (similar'ish to TM-pred):
       *
       *  A   B   C   D   E   F   G   H
       *  I   J   K   L
       *  M   N   O   P
       *  Q   R   S   T
       *  U   V   W   X   x1  x2  x3  x4
       *
       * Where:
       * x1 = clip_pixel(E + X - D)
       * x2 = clip_pixel(F + X - D)
       * x3 = clip_pixel(G + X - D)
       * x4 = clip_pixel(H + X - D)
       *
       * This is applied anytime when we use a "distant" above-right edge
       * that is not immediately top-right to the block that we're going
       * to do intra prediction for.
       */
      above[4] = clip_pixel(above_right[0] + above_ptr[3] - above_right[-1]);
      above[5] = clip_pixel(above_right[1] + above_ptr[3] - above_right[-1]);
      above[6] = clip_pixel(above_right[2] + above_ptr[3] - above_right[-1]);
      above[7] = clip_pixel(above_right[3] + above_ptr[3] - above_right[-1]);
    } else {
      // extend edge
      above[4] = above[5] = above[6] = above[7] = above[3];
    }
  } else {
    above[0] = above[1] = above[2] = above[3] = 127;
    above[4] = above[5] = above[6] = above[7] = 127;
    top_left = 127;
  }

#if CONFIG_NEWBINTRAMODES
  if (b_mode == B_CONTEXT_PRED)
    b_mode = x->bmi.as_mode.context;
#endif

  switch (b_mode) {
    case B_DC_PRED: {
      int expected_dc = 128;
      if (have_top || have_left) {
        int average = 0;
        int count = 0;
        if (have_top) {
          for (i = 0; i < 4; i++)
            average += above[i];
          count += 4;
        }
        if (have_left) {
          for (i = 0; i < 4; i++)
            average += left[i];
          count += 4;
        }
        expected_dc = (average + (count >> 1)) / count;
      }
      for (r = 0; r < 4; r++) {
        for (c = 0; c < 4; c++)
          predictor[c] = expected_dc;
        predictor += ps;
      }
    }
    break;
    case B_TM_PRED: {
      /* prediction similar to true_motion prediction */
      for (r = 0; r < 4; r++) {
        for (c = 0; c < 4; c++)
          predictor[c] = clip_pixel(above[c] - top_left + left[r]);
        predictor += ps;
      }
    }
    break;
    case B_V_PRED:
      for (r = 0; r < 4; r++) {
        for (c = 0; c < 4; c++)
          predictor[c] = above[c];
        predictor += ps;
      }
      break;
    case B_H_PRED:
      for (r = 0; r < 4; r++) {
        for (c = 0; c < 4; c++)
          predictor[c] = left[r];
        predictor += ps;
      }
      break;
    case B_D45_PRED: {
      uint8_t *p = above;

      predictor[0 * ps + 0] = ROUND_POWER_OF_TWO(p[0] + p[1] * 2 + p[2], 2);
      predictor[0 * ps + 1] =
        predictor[1 * ps + 0] = ROUND_POWER_OF_TWO(p[1] + p[2] * 2 + p[3], 2);
      predictor[0 * ps + 2] =
        predictor[1 * ps + 1] =
          predictor[2 * ps + 0] = ROUND_POWER_OF_TWO(p[2] + p[3] * 2 + p[4], 2);
      predictor[0 * ps + 3] =
        predictor[1 * ps + 2] =
          predictor[2 * ps + 1] =
            predictor[3 * ps + 0] =
              ROUND_POWER_OF_TWO(p[3] + p[4] * 2 + p[5], 2);
      predictor[1 * ps + 3] =
        predictor[2 * ps + 2] =
          predictor[3 * ps + 1] = ROUND_POWER_OF_TWO(p[4] + p[5] * 2 + p[6], 2);
      predictor[2 * ps + 3] =
        predictor[3 * ps + 2] = ROUND_POWER_OF_TWO(p[5] + p[6] * 2 + p[7], 2);
      predictor[3 * ps + 3] = ROUND_POWER_OF_TWO(p[6] + p[7] * 2 + p[7], 2);

    }
    break;
    case B_D135_PRED: {
      uint8_t p[9] = { left[3], left[2], left[1], left[0],
                       top_left,
                       above[0], above[1], above[2], above[3] };

      predictor[3 * ps + 0] = ROUND_POWER_OF_TWO(p[0] + p[1] * 2 + p[2], 2);
      predictor[3 * ps + 1] =
        predictor[2 * ps + 0] = ROUND_POWER_OF_TWO(p[1] + p[2] * 2 + p[3], 2);
      predictor[3 * ps + 2] =
        predictor[2 * ps + 1] =
          predictor[1 * ps + 0] = ROUND_POWER_OF_TWO(p[2] + p[3] * 2 + p[4], 2);
      predictor[3 * ps + 3] =
        predictor[2 * ps + 2] =
          predictor[1 * ps + 1] =
            predictor[0 * ps + 0] =
              ROUND_POWER_OF_TWO(p[3] + p[4] * 2 + p[5], 2);
      predictor[2 * ps + 3] =
        predictor[1 * ps + 2] =
          predictor[0 * ps + 1] = ROUND_POWER_OF_TWO(p[4] + p[5] * 2 + p[6], 2);
      predictor[1 * ps + 3] =
        predictor[0 * ps + 2] = ROUND_POWER_OF_TWO(p[5] + p[6] * 2 + p[7], 2);
      predictor[0 * ps + 3] = ROUND_POWER_OF_TWO(p[6] + p[7] * 2 + p[8], 2);

    }
    break;
    case B_D117_PRED: {
      uint8_t p[9] = { left[3], left[2], left[1], left[0],
                       top_left,
                       above[0], above[1], above[2], above[3] };

      predictor[3 * ps + 0] = ROUND_POWER_OF_TWO(p[1] + p[2] * 2 + p[3], 2);
      predictor[2 * ps + 0] = ROUND_POWER_OF_TWO(p[2] + p[3] * 2 + p[4], 2);
      predictor[3 * ps + 1] =
        predictor[1 * ps + 0] = ROUND_POWER_OF_TWO(p[3] + p[4] * 2 + p[5], 2);
      predictor[2 * ps + 1] =
        predictor[0 * ps + 0] = ROUND_POWER_OF_TWO(p[4] + p[5], 1);
      predictor[3 * ps + 2] =
        predictor[1 * ps + 1] = ROUND_POWER_OF_TWO(p[4] + p[5] * 2 + p[6], 2);
      predictor[2 * ps + 2] =
        predictor[0 * ps + 1] = ROUND_POWER_OF_TWO(p[5] + p[6], 1);
      predictor[3 * ps + 3] =
        predictor[1 * ps + 2] = ROUND_POWER_OF_TWO(p[5] + p[6] * 2 + p[7], 2);
      predictor[2 * ps + 3] =
        predictor[0 * ps + 2] = ROUND_POWER_OF_TWO(p[6] + p[7], 1);
      predictor[1 * ps + 3] = ROUND_POWER_OF_TWO(p[6] + p[7] * 2 + p[8], 2);
      predictor[0 * ps + 3] = ROUND_POWER_OF_TWO(p[7] + p[8], 1);

    }
    break;
    case B_D63_PRED: {
      uint8_t *p = above;

      predictor[0 * ps + 0] = ROUND_POWER_OF_TWO(p[0] + p[1], 1);
      predictor[1 * ps + 0] = ROUND_POWER_OF_TWO(p[0] + p[1] * 2 + p[2], 2);
      predictor[2 * ps + 0] =
        predictor[0 * ps + 1] = ROUND_POWER_OF_TWO(p[1] + p[2], 1);
      predictor[1 * ps + 1] =
        predictor[3 * ps + 0] = ROUND_POWER_OF_TWO(p[1] + p[2] * 2 + p[3], 2);
      predictor[2 * ps + 1] =
        predictor[0 * ps + 2] = ROUND_POWER_OF_TWO(p[2] + p[3], 1);
      predictor[3 * ps + 1] =
        predictor[1 * ps + 2] = ROUND_POWER_OF_TWO(p[2] + p[3] * 2 + p[4], 2);
      predictor[0 * ps + 3] =
        predictor[2 * ps + 2] = ROUND_POWER_OF_TWO(p[3] + p[4], 1);
      predictor[1 * ps + 3] =
        predictor[3 * ps + 2] = ROUND_POWER_OF_TWO(p[3] + p[4] * 2 + p[5], 2);
      predictor[2 * ps + 3] = ROUND_POWER_OF_TWO(p[4] + p[5] * 2 + p[6], 2);
      predictor[3 * ps + 3] = ROUND_POWER_OF_TWO(p[5] + p[6] * 2 + p[7], 2);
    }
    break;
    case B_D153_PRED: {
      uint8_t p[9] = { left[3], left[2], left[1], left[0],
                       top_left,
                       above[0], above[1], above[2], above[3] };

      predictor[3 * ps + 0] = ROUND_POWER_OF_TWO(p[0] + p[1], 1);
      predictor[3 * ps + 1] = ROUND_POWER_OF_TWO(p[0] + p[1] * 2 + p[2], 2);
      predictor[2 * ps + 0] =
        predictor[3 * ps + 2] = ROUND_POWER_OF_TWO(p[1] + p[2], 1);
      predictor[2 * ps + 1] =
        predictor[3 * ps + 3] = ROUND_POWER_OF_TWO(p[1] + p[2] * 2 + p[3], 2);
      predictor[2 * ps + 2] =
        predictor[1 * ps + 0] = ROUND_POWER_OF_TWO(p[2] + p[3], 1);
      predictor[2 * ps + 3] =
        predictor[1 * ps + 1] = ROUND_POWER_OF_TWO(p[2] + p[3] * 2 + p[4], 2);
      predictor[1 * ps + 2] =
        predictor[0 * ps + 0] = ROUND_POWER_OF_TWO(p[3] + p[4], 1);
      predictor[1 * ps + 3] =
        predictor[0 * ps + 1] = ROUND_POWER_OF_TWO(p[3] + p[4] * 2 + p[5], 2);
      predictor[0 * ps + 2] = ROUND_POWER_OF_TWO(p[4] + p[5] * 2 + p[6], 2);
      predictor[0 * ps + 3] = ROUND_POWER_OF_TWO(p[5] + p[6] * 2 + p[7], 2);
    }
    break;
    case B_D27_PRED: {
      uint8_t *p = left;
      predictor[0 * ps + 0] = ROUND_POWER_OF_TWO(p[0] + p[1], 1);
      predictor[0 * ps + 1] = ROUND_POWER_OF_TWO(p[0] + p[1] * 2 + p[2], 2);
      predictor[0 * ps + 2] =
        predictor[1 * ps + 0] = ROUND_POWER_OF_TWO(p[1] + p[2], 1);
      predictor[0 * ps + 3] =
        predictor[1 * ps + 1] = ROUND_POWER_OF_TWO(p[1] + p[2] * 2 + p[3], 2);
      predictor[1 * ps + 2] =
        predictor[2 * ps + 0] = ROUND_POWER_OF_TWO(p[2] + p[3], 1);
      predictor[1 * ps + 3] =
        predictor[2 * ps + 1] = ROUND_POWER_OF_TWO(p[2] + p[3] * 2 + p[3], 2);
      predictor[2 * ps + 2] =
        predictor[2 * ps + 3] =
          predictor[3 * ps + 0] =
            predictor[3 * ps + 1] =
              predictor[3 * ps + 2] =
                predictor[3 * ps + 3] = p[3];
    }
    break;

#if CONFIG_NEWBINTRAMODES
    case B_CONTEXT_PRED:
    break;
    /*
    case B_CORNER_PRED:
    corner_predictor(predictor, 16, 4, above, left);
    break;
    */
#endif
  }
}
