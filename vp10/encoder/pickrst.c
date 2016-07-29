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
#include <float.h>
#include <limits.h>
#include <math.h>

#include "./vpx_scale_rtcd.h"

#include "vpx_dsp/psnr.h"
#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem.h"

#include "vp10/common/onyxc_int.h"
#include "vp10/common/quant_common.h"

#include "vp10/encoder/encoder.h"
#include "vp10/encoder/quantize.h"
#include "vp10/encoder/picklpf.h"
#include "vp10/encoder/pickrst.h"

static int64_t try_restoration_frame(const YV12_BUFFER_CONFIG *sd,
                                     VP10_COMP *const cpi,
                                     RestorationInfo *rsi,
                                     int partial_frame) {
  VP10_COMMON *const cm = &cpi->common;
  int64_t filt_err;
  vp10_loop_restoration_frame(cm->frame_to_show, cm,
                              rsi, 1, partial_frame);
#if CONFIG_VPX_HIGHBITDEPTH
  if (cm->use_highbitdepth) {
    filt_err = vpx_highbd_get_y_sse(sd, cm->frame_to_show);
  } else {
    filt_err = vpx_get_y_sse(sd, cm->frame_to_show);
  }
#else
  filt_err = vpx_get_y_sse(sd, cm->frame_to_show);
#endif  // CONFIG_VPX_HIGHBITDEPTH

  // Re-instate the unfiltered frame
  vpx_yv12_copy_y(&cpi->last_frame_db, cm->frame_to_show);
  return filt_err;
}

static int search_bilateral_level(const YV12_BUFFER_CONFIG *sd,
                                  VP10_COMP *cpi,
                                  int filter_level, int partial_frame,
                                  double *best_cost_ret) {
  VP10_COMMON *const cm = &cpi->common;
  int i, restoration_best;
  int64_t err;
  double best_cost;
  double cost;
  const int restoration_level_bits = vp10_restoration_level_bits(&cpi->common);
  const int restoration_levels = 1 << restoration_level_bits;
  MACROBLOCK *x = &cpi->td.mb;
  int bits;
  RestorationInfo rsi;

  //  Make a copy of the unfiltered / processed recon buffer
  vpx_yv12_copy_y(cm->frame_to_show, &cpi->last_frame_uf);
  vp10_loop_filter_frame(cm->frame_to_show, cm, &cpi->td.mb.e_mbd, filter_level,
                         1, partial_frame);
  vpx_yv12_copy_y(cm->frame_to_show, &cpi->last_frame_db);

  restoration_best = -1;
  rsi.restoration_type = RESTORE_NONE;
  err = try_restoration_frame(sd, cpi, &rsi, partial_frame);
  bits = 0;
  best_cost = RDCOST_DBL(x->rdmult, x->rddiv,
                         (bits << (VP9_PROB_COST_SHIFT - 4)), err);
  for (i = 0; i < restoration_levels; ++i) {
    rsi.restoration_type = RESTORE_BILATERAL;
    rsi.restoration_level = i;
    err = try_restoration_frame(sd, cpi, &rsi, partial_frame);
    // Normally the rate is rate in bits * 256 and dist is sum sq err * 64
    // when RDCOST is used.  However below we just scale both in the correct
    // ratios appropriately but not exactly by these values.
    bits = restoration_level_bits;
    cost = RDCOST_DBL(x->rdmult, x->rddiv,
                      (bits << (VP9_PROB_COST_SHIFT - 4)), err);
    if (cost < best_cost) {
      restoration_best = i;
      best_cost = cost;
    }
  }
  if (best_cost_ret) *best_cost_ret = best_cost;
  vpx_yv12_copy_y(&cpi->last_frame_uf, cm->frame_to_show);
  return restoration_best;
}

static int search_filter_bilateral_level(const YV12_BUFFER_CONFIG *sd,
                                         VP10_COMP *cpi,
                                         int partial_frame,
                                         int *restoration_level,
                                         double *best_cost_ret) {
  const VP10_COMMON *const cm = &cpi->common;
  const struct loopfilter *const lf = &cm->lf;
  const int min_filter_level = 0;
  const int max_filter_level = vp10_get_max_filter_level(cpi);
  int filt_direction = 0;
  int filt_best, restoration_best;
  double best_err;
  int i;
  int bilateral_lev;

  // Start the search at the previous frame filter level unless it is now out of
  // range.
  int filt_mid = clamp(lf->filter_level, min_filter_level, max_filter_level);
  int filter_step = filt_mid < 16 ? 4 : filt_mid / 4;
  double ss_err[MAX_LOOP_FILTER + 1];

  // Set each entry to -1
  for (i = 0; i <= MAX_LOOP_FILTER; ++i)
    ss_err[i] = -1.0;

  bilateral_lev = search_bilateral_level(sd, cpi, filt_mid,
                                         partial_frame, &best_err);
  filt_best = filt_mid;
  restoration_best = bilateral_lev;
  ss_err[filt_mid] = best_err;

  while (filter_step > 0) {
    const int filt_high = VPXMIN(filt_mid + filter_step, max_filter_level);
    const int filt_low = VPXMAX(filt_mid - filter_step, min_filter_level);

    // Bias against raising loop filter in favor of lowering it.
    double bias = (best_err / (1 << (15 - (filt_mid / 8)))) * filter_step;

    if ((cpi->oxcf.pass == 2) && (cpi->twopass.section_intra_rating < 20))
      bias = (bias * cpi->twopass.section_intra_rating) / 20;

    // yx, bias less for large block size
    if (cm->tx_mode != ONLY_4X4)
      bias /= 2;

    if (filt_direction <= 0 && filt_low != filt_mid) {
      // Get Low filter error score
      if (ss_err[filt_low] < 0) {
        bilateral_lev = search_bilateral_level(
            sd, cpi, filt_low, partial_frame, &ss_err[filt_low]);
      }
      // If value is close to the best so far then bias towards a lower loop
      // filter value.
      if (ss_err[filt_low] < (best_err + bias)) {
        // Was it actually better than the previous best?
        if (ss_err[filt_low] < best_err) {
          best_err = ss_err[filt_low];
        }
        filt_best = filt_low;
        restoration_best = bilateral_lev;
      }
    }

    // Now look at filt_high
    if (filt_direction >= 0 && filt_high != filt_mid) {
      if (ss_err[filt_high] < 0) {
        bilateral_lev = search_bilateral_level(
            sd, cpi, filt_high, partial_frame, &ss_err[filt_high]);
      }
      // If value is significantly better than previous best, bias added against
      // raising filter value
      if (ss_err[filt_high] < (best_err - bias)) {
        best_err = ss_err[filt_high];
        filt_best = filt_high;
        restoration_best = bilateral_lev;
      }
    }

    // Half the step distance if the best filter value was the same as last time
    if (filt_best == filt_mid) {
      filter_step /= 2;
      filt_direction = 0;
    } else {
      filt_direction = (filt_best < filt_mid) ? -1 : 1;
      filt_mid = filt_best;
    }
  }

  // Update best error
  best_err = ss_err[filt_best];

  *restoration_level = restoration_best;
  if (best_cost_ret) *best_cost_ret = best_err;
  return filt_best;
}

static double find_average(uint8_t *src, int width, int height, int stride) {
  uint64_t sum = 0;
  double avg = 0;
  int i, j;
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++)
      sum += src[i * stride + j];
  avg = (double)sum / (height * width);
  return avg;
}

static void compute_stats(uint8_t *dgd, uint8_t *src, int width, int height,
                          int dgd_stride, int src_stride,
                          double *M, double *H) {
  int i, j, k, l;
  double Y[RESTORATION_WIN2];
  const double avg = find_average(dgd, width, height, dgd_stride);

  memset(M, 0, sizeof(*M) * RESTORATION_WIN2);
  memset(H, 0, sizeof(*H) * RESTORATION_WIN2 * RESTORATION_WIN2);
  for (i = RESTORATION_HALFWIN; i < height - RESTORATION_HALFWIN; i++) {
    for (j = RESTORATION_HALFWIN; j < width - RESTORATION_HALFWIN; j++) {
      const double X = (double)src[i * src_stride + j] - avg;
      int idx = 0;
      for (k = -RESTORATION_HALFWIN; k <= RESTORATION_HALFWIN; k++) {
        for (l = -RESTORATION_HALFWIN; l <= RESTORATION_HALFWIN; l++) {
          Y[idx] = (double)dgd[(i + l) * dgd_stride + (j + k)] - avg;
          idx++;
        }
      }
      for (k = 0; k < RESTORATION_WIN2; ++k) {
        M[k] += Y[k] * X;
        H[k * RESTORATION_WIN2 + k] += Y[k] * Y[k];
        for (l = k + 1; l < RESTORATION_WIN2; ++l) {
          double value = Y[k] * Y[l];
          H[k * RESTORATION_WIN2 + l] += value;
          H[l * RESTORATION_WIN2 + k] += value;
        }
      }
    }
  }
}

#if CONFIG_VPX_HIGHBITDEPTH
static double find_average_highbd(uint16_t *src,
                                  int width, int height, int stride) {
  uint64_t sum = 0;
  double avg = 0;
  int i, j;
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++)
      sum += src[i * stride + j];
  avg = (double)sum / (height * width);
  return avg;
}

static void compute_stats_highbd(
    uint8_t *dgd8, uint8_t *src8, int width, int height,
    int dgd_stride, int src_stride, double *M, double *H) {
  int i, j, k, l;
  double Y[RESTORATION_WIN2];
  uint16_t *src = CONVERT_TO_SHORTPTR(src8);
  uint16_t *dgd = CONVERT_TO_SHORTPTR(dgd8);
  const double avg = find_average_highbd(dgd, width, height, dgd_stride);

  memset(M, 0, sizeof(*M) * RESTORATION_WIN2);
  memset(H, 0, sizeof(*H) * RESTORATION_WIN2 * RESTORATION_WIN2);
  for (i = RESTORATION_HALFWIN; i < height - RESTORATION_HALFWIN; i++) {
    for (j = RESTORATION_HALFWIN; j < width - RESTORATION_HALFWIN; j++) {
      const double X = (double)src[i * src_stride + j] - avg;
      int idx = 0;
      for (k = -RESTORATION_HALFWIN; k <= RESTORATION_HALFWIN; k++) {
        for (l = -RESTORATION_HALFWIN; l <= RESTORATION_HALFWIN; l++) {
          Y[idx] = (double)dgd[(i + l) * dgd_stride + (j + k)] - avg;
          idx++;
        }
      }
      for (k = 0; k < RESTORATION_WIN2; ++k) {
        M[k] += Y[k] * X;
        H[k * RESTORATION_WIN2 + k] += Y[k] * Y[k];
        for (l = k + 1; l < RESTORATION_WIN2; ++l) {
          double value = Y[k] * Y[l];
          H[k * RESTORATION_WIN2 + l] += value;
          H[l * RESTORATION_WIN2 + k] += value;
        }
      }
    }
  }
}
#endif  // CONFIG_VPX_HIGHBITDEPTH

// Solves Ax = b, where x and b are column vectors
static int linsolve(int n, double *A, int stride, double *b, double *x) {
  int i, j, k;
  double c;
  // Partial pivoting
  for (i = n - 1; i > 0; i--) {
    if (A[(i - 1) * stride] < A[i * stride]) {
      for (j = 0; j < n; j++) {
        c = A[i * stride + j];
        A[i * stride + j] = A[(i - 1) * stride + j];
        A[(i - 1) * stride + j] = c;
      }
      c = b[i];
      b[i] = b[i - 1];
      b[i - 1] = c;
    }
  }
  // Forward elimination
  for (k = 0; k < n - 1; k++) {
    for (i = k; i < n - 1; i++) {
      c = A[(i + 1) * stride + k] / A[k * stride + k];
      for (j = 0; j < n; j++)
        A[(i + 1) * stride + j] -= c * A[k * stride + j];
      b[i + 1] -= c * b[k];
    }
  }
  // Backward substitution
  for (i = n - 1; i >= 0; i--) {
    if (fabs(A[i * stride + i]) < 1e-10)
      return 0;
    c = 0;
    for (j = i + 1; j <= n - 1; j++)
      c += A[i * stride + j] * x[j];
    x[i] = (b[i] - c) / A[i * stride + i];
  }
  return 1;
}

static INLINE int wrap_index(int i) {
  return (i >= RESTORATION_HALFWIN1 ? RESTORATION_WIN - 1 - i : i);
}

// Fix vector b, update vector a
static void update_a_sep_sym(double **Mc, double **Hc, double *a, double *b) {
  int i, j;
  double S[RESTORATION_WIN];
  double A[RESTORATION_WIN], B[RESTORATION_WIN2];
  int w, w2;
  memset(A, 0, sizeof(A));
  memset(B, 0, sizeof(B));
  for (i = 0; i < RESTORATION_WIN; i ++) {
    int j;
    for (j = 0; j < RESTORATION_WIN; ++j) {
      const int jj = wrap_index(j);
      A[jj] += Mc[i][j] * b[i];
    }
  }
  for (i = 0; i < RESTORATION_WIN; i ++) {
    for (j = 0; j < RESTORATION_WIN; j ++) {
      int k, l;
      for (k = 0; k < RESTORATION_WIN; ++k)
        for (l = 0; l < RESTORATION_WIN; ++l) {
          const int kk = wrap_index(k);
          const int ll = wrap_index(l);
          B[ll * RESTORATION_HALFWIN1 + kk] +=
              Hc[j * RESTORATION_WIN + i][k * RESTORATION_WIN2 + l] *
              b[i] * b[j];
        }
    }
  }
  // Normalization enforcement in the system of equations itself
  w = RESTORATION_WIN;
  w2 = (w >> 1) + 1;
  for (i = 0; i < w2 - 1; ++i)
    A[i] -= A[w2 - 1] * 2 + B[i * w2 + w2 - 1]
              - 2 * B[(w2 - 1) * w2 + (w2 - 1)];
  for (i = 0; i < w2 - 1; ++i)
    for (j = 0; j < w2 - 1; ++j)
      B[i * w2 + j] -= 2 * (B[i * w2 + (w2 - 1)] + B[(w2 - 1) * w2 + j] -
                            2 * B[(w2 - 1) * w2 + (w2 - 1)]);
  if (linsolve(w2 - 1, B, w2, A, S)) {
    S[w2 - 1] = 1.0;
    for (i = w2; i < w; ++i) {
      S[i] = S[w - 1 - i];
      S[w2 - 1] -= 2 * S[i];
    }
    memcpy(a, S, w * sizeof(*a));
  }
}

// Fix vector a, update vector b
static void update_b_sep_sym(double **Mc, double **Hc, double *a, double *b) {
  int i, j;
  double S[RESTORATION_WIN];
  double A[RESTORATION_WIN], B[RESTORATION_WIN2];
  int w, w2;
  memset(A, 0, sizeof(A));
  memset(B, 0, sizeof(B));
  for (i = 0; i < RESTORATION_WIN; i ++) {
    int j;
    const int ii = wrap_index(i);
    for (j = 0; j < RESTORATION_WIN; j ++)
      A[ii] += Mc[i][j] * a[j];
  }

  for (i = 0; i < RESTORATION_WIN; i++) {
    for (j = 0; j < RESTORATION_WIN; j++) {
      const int ii = wrap_index(i);
      const int jj = wrap_index(j);
      int k, l;
      for (k = 0; k < RESTORATION_WIN; ++k)
        for (l = 0; l < RESTORATION_WIN; ++l)
          B[jj * RESTORATION_HALFWIN1 + ii] +=
              Hc[i * RESTORATION_WIN + j][k * RESTORATION_WIN2 + l] *
              a[k] * a[l];
    }
  }
  // Normalization enforcement in the system of equations itself
  w = RESTORATION_WIN;
  w2 = RESTORATION_HALFWIN1;
  for (i = 0; i < w2 - 1; ++i)
    A[i] -= A[w2 - 1] * 2 + B[i * w2 + w2 - 1]
              - 2 * B[(w2 - 1) * w2 + (w2 - 1)];
  for (i = 0; i < w2 - 1; ++i)
    for (j = 0; j < w2 - 1; ++j)
      B[i * w2 + j] -= 2 * (B[i * w2 + (w2 - 1)] + B[(w2 - 1) * w2 + j] -
                            2 * B[(w2 - 1) * w2 + (w2 - 1)]);
  if (linsolve(w2 - 1, B, w2, A, S)) {
    S[w2 - 1] = 1.0;
    for (i = w2; i < w; ++i) {
      S[i] = S[w - 1 - i];
      S[w2 - 1] -= 2 * S[i];
    }
    memcpy(b, S, w * sizeof(*b));
  }
}

static int wiener_decompose_sep_sym(double *M, double *H,
                                    double *a, double *b) {
  static const double init_filt[RESTORATION_WIN] = {
    0.035623, -0.127154,  0.211436,  0.760190,  0.211436, -0.127154,  0.035623,
  };
  int i, j, iter;
  double *Hc[RESTORATION_WIN2];
  double *Mc[RESTORATION_WIN];
  for (i = 0; i < RESTORATION_WIN; i++) {
    Mc[i] = M + i * RESTORATION_WIN;
    for (j = 0; j < RESTORATION_WIN; j++) {
      Hc[i * RESTORATION_WIN + j] =
          H + i * RESTORATION_WIN * RESTORATION_WIN2 + j * RESTORATION_WIN;
    }
  }
  memcpy(a, init_filt, sizeof(*a) * RESTORATION_WIN);
  memcpy(b, init_filt, sizeof(*b) * RESTORATION_WIN);

  iter = 1;
  while (iter < 10) {
    update_a_sep_sym(Mc, Hc, a, b);
    update_b_sep_sym(Mc, Hc, a, b);
    iter++;
  }
  return 1;
}

// Computes the function x'*A*x - x'*b for the learned filters, and compares
// against identity filters; Final score is defined as the difference between
// the function values
  static double compute_score(double *M, double *H, int *vfilt, int *hfilt) {
  double ab[RESTORATION_WIN * RESTORATION_WIN];
  int i, k, l;
  double P = 0, Q = 0;
  double iP = 0, iQ = 0;
  double Score, iScore;
  int w;
  double a[RESTORATION_WIN], b[RESTORATION_WIN];
  w = RESTORATION_WIN;
  a[RESTORATION_HALFWIN] = b[RESTORATION_HALFWIN] = 1.0;
  for (i = 0; i < RESTORATION_HALFWIN; ++i) {
    a[i] = a[RESTORATION_WIN - i - 1 ] =
        (double) vfilt[i] / RESTORATION_FILT_STEP;
    b[i] = b[RESTORATION_WIN - i - 1 ] =
        (double) hfilt[i] / RESTORATION_FILT_STEP;
    a[RESTORATION_HALFWIN] -= 2 * a[i];
    b[RESTORATION_HALFWIN] -= 2 * b[i];
  }
  for (k = 0; k < w; ++k) {
    for (l = 0; l < w; ++l) {
      ab[k * w + l] = a[l] * b[k];
    }
  }
  for (k = 0; k < w * w; ++k) {
    P += ab[k] * M[k];
    for (l = 0; l < w * w; ++l)
      Q += ab[k] * H[k * w * w + l] * ab[l];
  }
  Score = Q - 2 * P;

  iP = M[(w * w) >> 1];
  iQ = H[((w * w) >> 1) * w * w + ((w * w) >> 1)];
  iScore = iQ - 2 * iP;

  return Score - iScore;
}

#define CLIP(x, lo, hi) ((x) < (lo) ? (lo) : (x) > (hi) ? (hi) : (x))
#define RINT(x) ((x) < 0 ? (int)((x) - 0.5) : (int)((x) + 0.5))

static void quantize_sym_filter(double *f, int *fi) {
  int i;
  for (i = 0; i < RESTORATION_HALFWIN; ++i) {
    fi[i] = RINT(f[i] * RESTORATION_FILT_STEP);
  }
  // Specialize for 7-tap filter
  fi[0] = CLIP(fi[0], WIENER_FILT_TAP0_MINV, WIENER_FILT_TAP0_MAXV);
  fi[1] = CLIP(fi[1], WIENER_FILT_TAP1_MINV, WIENER_FILT_TAP1_MAXV);
  fi[2] = CLIP(fi[2], WIENER_FILT_TAP2_MINV, WIENER_FILT_TAP2_MAXV);
}

static int search_wiener_filter(const YV12_BUFFER_CONFIG *src,
                                VP10_COMP *cpi,
                                int filter_level,
                                int partial_frame,
                                int *vfilter, int *hfilter,
                                double *best_cost_ret) {
  VP10_COMMON *const cm = &cpi->common;
  RestorationInfo rsi;
  int64_t err;
  int bits;
  double cost_wiener, cost_norestore;
  MACROBLOCK *x = &cpi->td.mb;
  double M[RESTORATION_WIN2];
  double H[RESTORATION_WIN2 * RESTORATION_WIN2];
  double vfilterd[RESTORATION_WIN], hfilterd[RESTORATION_WIN];
  const YV12_BUFFER_CONFIG *dgd = cm->frame_to_show;
  const int width = cm->width;
  const int height = cm->height;
  const int src_stride = src->y_stride;
  const int dgd_stride = dgd->y_stride;
  double score;

  assert(width == dgd->y_crop_width);
  assert(height == dgd->y_crop_height);
  assert(width == src->y_crop_width);
  assert(height == src->y_crop_height);

  //  Make a copy of the unfiltered / processed recon buffer
  vpx_yv12_copy_y(cm->frame_to_show, &cpi->last_frame_uf);
  vp10_loop_filter_frame(cm->frame_to_show, cm, &cpi->td.mb.e_mbd, filter_level,
                         1, partial_frame);
  vpx_yv12_copy_y(cm->frame_to_show, &cpi->last_frame_db);

  rsi.restoration_type = RESTORE_NONE;
  err = try_restoration_frame(src, cpi, &rsi, partial_frame);
  bits = 0;
  cost_norestore = RDCOST_DBL(x->rdmult, x->rddiv,
                              (bits << (VP9_PROB_COST_SHIFT - 4)), err);

#if CONFIG_VPX_HIGHBITDEPTH
  if (cm->use_highbitdepth)
    compute_stats_highbd(dgd->y_buffer, src->y_buffer, width, height,
                         dgd_stride, src_stride, M, H);
  else
#endif  // CONFIG_VPX_HIGHBITDEPTH
    compute_stats(dgd->y_buffer, src->y_buffer, width, height,
                  dgd_stride, src_stride, M, H);

  if (!wiener_decompose_sep_sym(M, H, vfilterd, hfilterd)) {
    *best_cost_ret = DBL_MAX;
    return 0;
  }
  quantize_sym_filter(vfilterd, vfilter);
  quantize_sym_filter(hfilterd, hfilter);

  // Filter score computes the value of the function x'*A*x - x'*b for the
  // learned filter and compares it against identity filer. If there is no
  // reduction in the function, the filter is reverted back to identity
  score = compute_score(M, H, vfilter, hfilter);
  if (score > 0.0) {
    int i;
    for (i = 0; i < RESTORATION_HALFWIN; ++i)
      vfilter[i] = hfilter[i] = 0;
    rsi.restoration_type = RESTORE_NONE;
    if (best_cost_ret) *best_cost_ret = cost_norestore;
    vpx_yv12_copy_y(&cpi->last_frame_uf, cm->frame_to_show);
    return 0;
  }

  rsi.restoration_type = RESTORE_WIENER;
  memcpy(rsi.vfilter, vfilter, sizeof(rsi.vfilter));
  memcpy(rsi.hfilter, hfilter, sizeof(rsi.hfilter));
  err = try_restoration_frame(src, cpi, &rsi, partial_frame);
  bits = WIENER_FILT_BITS;
  cost_wiener = RDCOST_DBL(x->rdmult, x->rddiv,
                           (bits << (VP9_PROB_COST_SHIFT - 4)), err);

  vpx_yv12_copy_y(&cpi->last_frame_uf, cm->frame_to_show);

  if (cost_wiener < cost_norestore) {
    if (best_cost_ret) *best_cost_ret = cost_wiener;
    return 1;
  } else {
    if (best_cost_ret) *best_cost_ret = cost_norestore;
    return 0;
  }
}

void vp10_pick_filter_restoration(
    const YV12_BUFFER_CONFIG *sd, VP10_COMP *cpi, LPF_PICK_METHOD method) {
  VP10_COMMON *const cm = &cpi->common;
  struct loopfilter *const lf = &cm->lf;
  int wiener_success = 0;
  double cost_bilateral = DBL_MAX;
  double cost_wiener = DBL_MAX;
  double cost_norestore = DBL_MAX;

  lf->sharpness_level =
      cm->frame_type == KEY_FRAME ? 0 : cpi->oxcf.sharpness;

  if (method == LPF_PICK_MINIMAL_LPF && lf->filter_level) {
      lf->filter_level = 0;
      cm->rst_info.restoration_type = RESTORE_NONE;
  } else if (method >= LPF_PICK_FROM_Q) {
    const int min_filter_level = 0;
    const int max_filter_level = vp10_get_max_filter_level(cpi);
    const int q = vp10_ac_quant(cm->base_qindex, 0, cm->bit_depth);
    // These values were determined by linear fitting the result of the
    // searched level, filt_guess = q * 0.316206 + 3.87252
#if CONFIG_VPX_HIGHBITDEPTH
    int filt_guess;
    switch (cm->bit_depth) {
      case VPX_BITS_8:
        filt_guess = ROUND_POWER_OF_TWO(q * 20723 + 1015158, 18);
        break;
      case VPX_BITS_10:
        filt_guess = ROUND_POWER_OF_TWO(q * 20723 + 4060632, 20);
        break;
      case VPX_BITS_12:
        filt_guess = ROUND_POWER_OF_TWO(q * 20723 + 16242526, 22);
        break;
      default:
        assert(0 && "bit_depth should be VPX_BITS_8, VPX_BITS_10 "
                    "or VPX_BITS_12");
        return;
    }
#else
    int filt_guess = ROUND_POWER_OF_TWO(q * 20723 + 1015158, 18);
#endif  // CONFIG_VPX_HIGHBITDEPTH
    if (cm->frame_type == KEY_FRAME)
      filt_guess -= 4;
    lf->filter_level = clamp(filt_guess, min_filter_level, max_filter_level);
    cm->rst_info.restoration_level = search_bilateral_level(
        sd, cpi, lf->filter_level, method == LPF_PICK_FROM_SUBIMAGE,
        &cost_bilateral);
    wiener_success = search_wiener_filter(
        sd, cpi, lf->filter_level, method == LPF_PICK_FROM_SUBIMAGE,
        cm->rst_info.vfilter, cm->rst_info.hfilter, &cost_wiener);
    if (cost_bilateral < cost_wiener) {
      if (cm->rst_info.restoration_level != -1)
        cm->rst_info.restoration_type = RESTORE_BILATERAL;
      else
        cm->rst_info.restoration_type = RESTORE_NONE;
    } else {
      if (wiener_success)
        cm->rst_info.restoration_type = RESTORE_WIENER;
      else
        cm->rst_info.restoration_type = RESTORE_NONE;
    }
  } else {
    int blf_filter_level = -1;
    blf_filter_level = search_filter_bilateral_level(
        sd, cpi, method == LPF_PICK_FROM_SUBIMAGE,
        &cm->rst_info.restoration_level, &cost_bilateral);
    lf->filter_level = vp10_search_filter_level(
        sd, cpi, method == LPF_PICK_FROM_SUBIMAGE, &cost_norestore);
    wiener_success = search_wiener_filter(
        sd, cpi, lf->filter_level, method == LPF_PICK_FROM_SUBIMAGE,
        cm->rst_info.vfilter, cm->rst_info.hfilter, &cost_wiener);
    if (cost_bilateral < cost_wiener) {
      lf->filter_level = blf_filter_level;
      if (cm->rst_info.restoration_level != -1)
        cm->rst_info.restoration_type = RESTORE_BILATERAL;
      else
        cm->rst_info.restoration_type = RESTORE_NONE;
    } else {
      if (wiener_success)
        cm->rst_info.restoration_type = RESTORE_WIENER;
      else
        cm->rst_info.restoration_type = RESTORE_NONE;
    }
    // printf("[%d] Costs %g %g (%d) %g (%d)\n", cm->rst_info.restoration_type,
    //         cost_norestore, cost_bilateral, lf->filter_level, cost_wiener,
    //         wiener_success);
  }
}
