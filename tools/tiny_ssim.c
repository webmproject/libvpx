/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vpx/vpx_integer.h"

void vp8_ssim_parms_8x8_c(unsigned char *s, int sp, unsigned char *r, int rp,
                          uint32_t *sum_s, uint32_t *sum_r, uint32_t *sum_sq_s,
                          uint32_t *sum_sq_r, uint32_t *sum_sxr) {
  int i, j;
  for (i = 0; i < 8; i++, s += sp, r += rp) {
    for (j = 0; j < 8; j++) {
      *sum_s += s[j];
      *sum_r += r[j];
      *sum_sq_s += s[j] * s[j];
      *sum_sq_r += r[j] * r[j];
      *sum_sxr += s[j] * r[j];
    }
  }
}

static const int64_t cc1 = 26634;   // (64^2*(.01*255)^2
static const int64_t cc2 = 239708;  // (64^2*(.03*255)^2

static double similarity(uint32_t sum_s, uint32_t sum_r, uint32_t sum_sq_s,
                         uint32_t sum_sq_r, uint32_t sum_sxr, int count) {
  int64_t ssim_n, ssim_d;
  int64_t c1, c2;

  // scale the constants by number of pixels
  c1 = (cc1 * count * count) >> 12;
  c2 = (cc2 * count * count) >> 12;

  ssim_n = (2 * sum_s * sum_r + c1) *
           ((int64_t)2 * count * sum_sxr - (int64_t)2 * sum_s * sum_r + c2);

  ssim_d = (sum_s * sum_s + sum_r * sum_r + c1) *
           ((int64_t)count * sum_sq_s - (int64_t)sum_s * sum_s +
            (int64_t)count * sum_sq_r - (int64_t)sum_r * sum_r + c2);

  return ssim_n * 1.0 / ssim_d;
}

static double ssim_8x8(unsigned char *s, int sp, unsigned char *r, int rp) {
  uint32_t sum_s = 0, sum_r = 0, sum_sq_s = 0, sum_sq_r = 0, sum_sxr = 0;
  vp8_ssim_parms_8x8_c(s, sp, r, rp, &sum_s, &sum_r, &sum_sq_s, &sum_sq_r,
                       &sum_sxr);
  return similarity(sum_s, sum_r, sum_sq_s, sum_sq_r, sum_sxr, 64);
}

// We are using a 8x8 moving window with starting location of each 8x8 window
// on the 4x4 pixel grid. Such arrangement allows the windows to overlap
// block boundaries to penalize blocking artifacts.
double vp8_ssim2(unsigned char *img1, unsigned char *img2, int stride_img1,
                 int stride_img2, int width, int height) {
  int i, j;
  int samples = 0;
  double ssim_total = 0;

  // sample point start with each 4x4 location
  for (i = 0; i <= height - 8;
       i += 4, img1 += stride_img1 * 4, img2 += stride_img2 * 4) {
    for (j = 0; j <= width - 8; j += 4) {
      double v = ssim_8x8(img1 + j, stride_img1, img2 + j, stride_img2);
      ssim_total += v;
      samples++;
    }
  }
  ssim_total /= samples;
  return ssim_total;
}

static uint64_t calc_plane_error(uint8_t *orig, int orig_stride, uint8_t *recon,
                                 int recon_stride, unsigned int cols,
                                 unsigned int rows) {
  unsigned int row, col;
  uint64_t total_sse = 0;
  int diff;

  for (row = 0; row < rows; row++) {
    for (col = 0; col < cols; col++) {
      diff = orig[col] - recon[col];
      total_sse += diff * diff;
    }

    orig += orig_stride;
    recon += recon_stride;
  }

  return total_sse;
}

#define MAX_PSNR 100

double vp9_mse2psnr(double samples, double peak, double mse) {
  double psnr;

  if (mse > 0.0)
    psnr = 10.0 * log10(peak * peak * samples / mse);
  else
    psnr = MAX_PSNR;  // Limit to prevent / 0

  if (psnr > MAX_PSNR) psnr = MAX_PSNR;

  return psnr;
}

int main(int argc, char *argv[]) {
  FILE *f[2];
  uint8_t *buf[2];
  int w, h, tl_skip = 0, tl_skips_remaining = 0;
  double ssim = 0, psnravg = 0, psnrglb = 0;
  double psnryglb = 0, psnruglb = 0, psnrvglb = 0;
  double ssimyavg = 0, ssimuavg = 0, ssimvavg = 0;
  double psnryavg = 0, psnruavg = 0, psnrvavg = 0;
  double *ssimy = NULL, *ssimu = NULL, *ssimv = NULL;
  uint64_t *psnry = NULL, *psnru = NULL, *psnrv = NULL;
  size_t i, n_frames = 0, allocated_frames = 0;

  if (argc < 4) {
    fprintf(stderr, "Usage: %s file1.yuv file2.yuv WxH [tl_skip={0,1,3}]\n",
            argv[0]);
    return 1;
  }
  f[0] = strcmp(argv[1], "-") ? fopen(argv[1], "rb") : stdin;
  f[1] = strcmp(argv[2], "-") ? fopen(argv[2], "rb") : stdin;
  sscanf(argv[3], "%dx%d", &w, &h);
  // Number of frames to skip from file1.yuv for every frame used. Normal values
  // 0, 1 and 3 correspond to TL2, TL1 and TL0 respectively for a 3TL encoding
  // in mode 10. 7 would be reasonable for comparing TL0 of a 4-layer encoding.
  if (argc > 4) {
    sscanf(argv[4], "%d", &tl_skip);
  }
  if (!f[0] || !f[1]) {
    fprintf(stderr, "Could not open input files: %s\n", strerror(errno));
    return 1;
  }
  if (w <= 0 || h <= 0 || w & 1 || h & 1) {
    fprintf(stderr, "Invalid size %dx%d\n", w, h);
    return 1;
  }
  buf[0] = malloc(w * h * 3 / 2);
  buf[1] = malloc(w * h * 3 / 2);
  while (1) {
    size_t r1, r2;
    r1 = fread(buf[0], w * h * 3 / 2, 1, f[0]);
    if (r1) {
      // Reading parts of file1.yuv that were not used in temporal layer.
      if (tl_skips_remaining > 0) {
        --tl_skips_remaining;
        continue;
      }
      // Use frame, but skip |tl_skip| after it.
      tl_skips_remaining = tl_skip;
    }
    r2 = fread(buf[1], w * h * 3 / 2, 1, f[1]);
    if (r1 && r2 && r1 != r2) {
      fprintf(stderr, "Failed to read data: %s [%d/%d]\n", strerror(errno),
              (int)r1, (int)r2);
      return 1;
    } else if (r1 == 0 || r2 == 0) {
      break;
    }
#define psnr_and_ssim(ssim, psnr, buf0, buf1, w, h) \
  ssim = vp8_ssim2(buf0, buf1, w, w, w, h);         \
  psnr = calc_plane_error(buf0, w, buf1, w, w, h);
    if (n_frames == allocated_frames) {
      allocated_frames = allocated_frames == 0 ? 1024 : allocated_frames * 2;
      ssimy = realloc(ssimy, allocated_frames * sizeof(*ssimy));
      ssimu = realloc(ssimu, allocated_frames * sizeof(*ssimu));
      ssimv = realloc(ssimv, allocated_frames * sizeof(*ssimv));

      psnry = realloc(psnry, allocated_frames * sizeof(*psnry));
      psnru = realloc(psnru, allocated_frames * sizeof(*psnru));
      psnrv = realloc(psnrv, allocated_frames * sizeof(*psnrv));
    }
    psnr_and_ssim(ssimy[n_frames], psnry[n_frames], buf[0], buf[1], w, h);
    psnr_and_ssim(ssimu[n_frames], psnru[n_frames], buf[0] + w * h,
                  buf[1] + w * h, w / 2, h / 2);
    psnr_and_ssim(ssimv[n_frames], psnrv[n_frames], buf[0] + w * h * 5 / 4,
                  buf[1] + w * h * 5 / 4, w / 2, h / 2);
    n_frames++;
  }
  free(buf[0]);
  free(buf[1]);

  for (i = 0; i < n_frames; ++i) {
    ssimyavg += ssimy[i];
    ssimuavg += ssimu[i];
    ssimvavg += ssimv[i];

    psnryavg += vp9_mse2psnr(w * h * 4 / 4, 255.0, (double)psnry[i]);
    psnruavg += vp9_mse2psnr(w * h * 1 / 4, 255.0, (double)psnru[i]);
    psnrvavg += vp9_mse2psnr(w * h * 1 / 4, 255.0, (double)psnrv[i]);

    psnravg += vp9_mse2psnr(w * h * 6 / 4, 255.0,
                            (double)psnry[i] + psnru[i] + psnrv[i]);
    psnryglb += psnry[i];
    psnruglb += psnru[i];
    psnrvglb += psnrv[i];
  }

  ssim = 0.8 * ssimyavg + 0.1 * (ssimuavg + ssimvavg);
  ssim /= n_frames;
  ssimyavg /= n_frames;
  ssimuavg /= n_frames;
  ssimvavg /= n_frames;

  printf("VpxSSIM: %lf\n", 100 * pow(ssim, 8.0));
  printf("SSIM: %lf\n", ssim);
  printf("SSIM-Y: %lf\n", ssimyavg);
  printf("SSIM-U: %lf\n", ssimuavg);
  printf("SSIM-V: %lf\n", ssimvavg);
  puts("");

  psnravg /= n_frames;
  psnryavg /= n_frames;
  psnruavg /= n_frames;
  psnrvavg /= n_frames;

  printf("AvgPSNR: %lf\n", psnravg);
  printf("AvgPSNR-Y: %lf\n", psnryavg);
  printf("AvgPSNR-U: %lf\n", psnruavg);
  printf("AvgPSNR-V: %lf\n", psnrvavg);
  puts("");

  psnrglb = psnryglb + psnruglb + psnrvglb;
  psnrglb = vp9_mse2psnr((double)n_frames * w * h * 6 / 4, 255.0, psnrglb);
  psnryglb = vp9_mse2psnr((double)n_frames * w * h * 4 / 4, 255.0, psnryglb);
  psnruglb = vp9_mse2psnr((double)n_frames * w * h * 1 / 4, 255.0, psnruglb);
  psnrvglb = vp9_mse2psnr((double)n_frames * w * h * 1 / 4, 255.0, psnrvglb);

  printf("GlbPSNR: %lf\n", psnrglb);
  printf("GlbPSNR-Y: %lf\n", psnryglb);
  printf("GlbPSNR-U: %lf\n", psnruglb);
  printf("GlbPSNR-V: %lf\n", psnrvglb);
  puts("");

  printf("Nframes: %d\n", (int)n_frames);

  // TODO(pbos): Add command-line flag for printing per-frame SSIM/PSNR metrics.

  if (strcmp(argv[1], "-")) fclose(f[0]);
  if (strcmp(argv[2], "-")) fclose(f[1]);

  free(ssimy);
  free(ssimu);
  free(ssimv);

  free(psnry);
  free(psnru);
  free(psnrv);

  return 0;
}
