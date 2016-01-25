/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be
 *  found  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <memory.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "vp9_ransac.h"

#define MAX_PARAMDIM 9
#define MAX_MINPTS   4

#define MAX_DEGENERATE_ITER 10
#define MINPTS_MULTIPLIER   5

// svdcmp
// Adopted from Numerical Recipes in C

static const double TINY_NEAR_ZERO = 1.0E-12;

static inline double SIGN(double a, double b) {
  return ((b) >= 0 ? fabs(a) : -fabs(a));
}

static inline double PYTHAG(double a, double b) {
  double absa, absb, ct;
  absa = fabs(a);
  absb = fabs(b);

  if(absa > absb) {
    ct = absb / absa;
    return absa * sqrt(1.0 + ct * ct);
  } else {
    ct = absa / absb;
    return (absb == 0) ? 0 : absb * sqrt(1.0 + ct * ct);
  }
}

int IMIN(int a, int b) {
  return (((a) < (b)) ? (a) : (b));
}

int IMAX(int a, int b) {
  return (((a) < (b)) ? (b) : (a));
}

void MultiplyMat(double *m1, double *m2, double *res,
                 const int M1, const int N1, const int N2) {
  int timesInner = N1;
  int timesRows = M1;
  int timesCols = N2;
  double sum;

  int row, col, inner;
  for( row = 0; row < timesRows; ++row ) {
    for( col = 0; col < timesCols; ++col ) {
      sum = 0;
      for (inner = 0; inner < timesInner; ++inner )
        sum += m1[row * N1 + inner] * m2[inner * N2 + col];
      *(res++) = sum;
    }
  }
}

static int svdcmp_(double **u, int m, int n, double w[], double **v) {
  const int max_its = 30;
  int flag, i, its, j, jj, k, l, nm;
  double anorm, c, f, g, h, s, scale, x, y, z;
  double *rv1 = (double *)malloc(sizeof(double) * (n + 1));
  g = scale = anorm = 0.0;
  for (i = 0; i < n; i++) {
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;
    if (i < m) {
      for (k = i; k < m; k++) scale += fabs(u[k][i]);
      if (scale) {
        for (k = i; k < m; k++) {
          u[k][i] /= scale;
          s += u[k][i] * u[k][i];
        }
        f = u[i][i];
        g = -SIGN(sqrt(s), f);
        h = f * g - s;
        u[i][i] = f - g;
        for (j = l; j < n; j++) {
          for (s = 0.0, k = i; k < m; k++) s += u[k][i] * u[k][j];
          f = s / h;
          for (k = i; k < m; k++) u[k][j] += f * u[k][i];
        }
        for (k = i; k < m; k++) u[k][i] *= scale;
      }
    }
    w[i] = scale * g;
    g = s = scale = 0.0;
    if (i < m && i != n - 1) {
      for (k = l; k < n; k++)
        scale += fabs(u[i][k]);
      if (scale) {
        for (k = l; k < n; k++) {
          u[i][k] /= scale;
          s += u[i][k] * u[i][k];
        }
        f = u[i][l];
        g = -SIGN(sqrt(s),f);
        h = f * g - s;
        u[i][l] = f - g;
        for (k = l; k < n; k++) rv1[k] = u[i][k] / h;
        for (j = l; j < m; j++) {
          for (s = 0.0, k = l; k < n; k++) s += u[j][k] * u[i][k];
          for (k = l; k < n; k++) u[j][k] += s * rv1[k];
        }
        for (k = l; k < n; k++) u[i][k] *= scale;
      }
    }
    anorm = fmax(anorm, (fabs(w[i]) + fabs(rv1[i])));
  }

  for (i = n - 1; i >= 0; i--) {
    if (i < n - 1) {
      if (g) {
        for (j = l; j < n; j++) v[j][i] = (u[i][j] / u[i][l]) / g;
        for (j = l; j < n; j++) {
          for (s = 0.0, k = l; k < n; k++) s += u[i][k] * v[k][j];
          for (k = l; k < n; k++) v[k][j] += s * v[k][i];
        }
      }
      for (j = l; j < n; j++) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }

  for (i = IMIN(m, n) - 1; i >= 0; i--) {
    l = i + 1;
    g = w[i];
    for (j = l; j < n; j++) u[i][j] = 0.0;
    if (g) {
      g = 1.0 / g;
      for (j = l; j < n; j++) {
        for (s = 0.0, k = l; k < m; k++) s += u[k][i] * u[k][j];
        f = (s / u[i][i]) * g;
        for (k = i; k < m; k++) u[k][j] += f * u[k][i];
      }
      for (j = i; j < m; j++) u[j][i] *= g;
    } else {
      for (j = i; j < m; j++) u[j][i] = 0.0;
    }
    ++u[i][i];
  }
  for (k = n - 1; k >= 0; k--) {
    for (its = 0; its < max_its; its++) {
      flag = 1;
      for (l = k; l >= 0; l--) {
        nm = l - 1;
        if ((double)(fabs(rv1[l]) + anorm) == anorm || nm < 0) {
          flag = 0;
          break;
        }
        if ((double)(fabs(w[nm]) + anorm) == anorm) break;
      }
      if (flag) {
        c = 0.0;
        s = 1.0;
        for (i = l; i <= k; i++) {
          f = s * rv1[i];
          rv1[i] = c * rv1[i];
          if ((double)(fabs(f) + anorm) == anorm) break;
          g = w[i];
          h = PYTHAG(f, g);
          w[i] = h;
          h = 1.0 / h;
          c = g * h;
          s = -f * h;
          for (j = 0; j < m; j++) {
            y = u[j][nm];
            z = u[j][i];
            u[j][nm] = y * c + z * s;
            u[j][i] = z * c - y * s;
          }
        }
      }
      z = w[k];
      if (l == k) {
        if (z < 0.0) {
          w[k] = -z;
          for (j = 0; j < n; j++) v[j][k] = -v[j][k];
        }
        break;
      }
      if (its == max_its - 1) {
        return 1;
      }
      assert(k > 0);
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = PYTHAG(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
      c = s = 1.0;
      for (j = l; j <= nm; j++) {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = PYTHAG(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y *= c;
        for (jj = 0; jj < n; jj++) {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = x * c + z * s;
          v[jj][i] = z * c - x * s;
        }
        z = PYTHAG(f, h);
        w[j] = z;
        if (z) {
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }
        f = c * g + s * y;
        x = c * y - s * g;
        for (jj = 0; jj < m; jj++) {
          y = u[jj][j];
          z = u[jj][i];
          u[jj][j] = y * c + z * s;
          u[jj][i] = z * c - y * s;
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  free(rv1);
  return 0;
}

static int SVD(double *U, double *W, double *V, double *matx, int M, int N) {
  // Assumes allocation for U is MxN
  double **nrU, **nrV;
  int problem, i;

  nrU = (double **)malloc((M)*sizeof(double*));
  nrV = (double **)malloc((N)*sizeof(double*));
  problem = !(nrU && nrV);
  if (!problem) {
    problem = 0;
    for (i = 0; i < M; i++) {
      nrU[i] = &U[i * N];
    }
    for (i = 0; i < N; i++) {
      nrV[i] = &V[i * N];
    }
  }
  if (problem) {
    return 1;
  }

  /* copy from given matx into nrU */
  for (i = 0; i < M; i++) {
    memcpy(&(nrU[i][0]), matx + N * i, N * sizeof(*matx));
  }

  /* HERE IT IS: do SVD */
  if (svdcmp_(nrU, M, N, W, nrV)) {
    return 1;
  }

  /* free Numerical Recipes arrays */
  free(nrU);
  free(nrV);

  return 0;
}

int PseudoInverse(double *inv, double *matx, const int M, const int N) {
  double *U, *W, *V, ans;
  int i, j, k;
  U = (double *)malloc(M * N * sizeof(*matx));
  W = (double *)malloc(N * sizeof(*matx));
  V = (double *)malloc(N * N * sizeof(*matx));

  if (!(U && W && V)) {
    return 1;
  }
  if (SVD(U, W, V, matx, M, N)) {
    return 1;
  }
  for (i = 0; i < N; i++) {
    if (fabs(W[i]) < TINY_NEAR_ZERO) {
      return 1;
    }
  }

  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++) {
      ans = 0;
      for (k = 0; k < N; k++) {
        ans += V[k + N * i] * U[k + N * j] / W[k];
      }
      inv[j + M * i] = ans;
    }
  }
  free(U);
  free(W);
  free(V);
  return 0;
}

static double compute_error(projectPointsType projectPoints,
                            double *points1, int stride1,
                            double *points2, int stride2,
                            int npoints, double *H, int *mask) {
  int i, n = 0;
  double pt[2];
  double *mp1 = points1;
  double *mp2 = points2;
  double sqerr = 0.0;

  if (projectPoints == NULL) return -1.0;
  if (mask) {
    for (i = 0; i < npoints; ++i, mp1 += stride1, mp2 += stride2) {
      if (mask[i]) {
        projectPoints(H, mp1, pt, 1, stride1, stride2);
        sqerr += (pt[0] - mp2[0]) * (pt[0] - mp2[0]) +
                 (pt[1] - mp2[1]) * (pt[1] - mp2[1]);
        n++;
      }
    }
  } else {
    for (i = 0; i < npoints; ++i, mp1 += stride1, mp2 += stride2) {
      projectPoints(H, mp1, pt, 1, stride1, stride2);
      sqerr += (pt[0] - mp2[0]) * (pt[0] - mp2[0]) +
               (pt[1] - mp2[1]) * (pt[1] - mp2[1]);
      n++;
    }
  }
  return sqrt(sqerr / n);
}

////////////////////////////////////////////////////////////////////////////////
// ransac
typedef int  (*isDegenerateType)(double *p);
typedef void (*normalizeType)(double *p, int np, double *T);
typedef void (*denormalizeType)(double *H, double *T1, double *T2);
typedef int  (*findTransformationType)(int points,
                                       double *points1,
                                       double *points2,
                                       double *H);

static int get_rand_indices(int npoints, int minpts, int *indices) {
  int i, j;
  int ptr = rand() % npoints;
  if (minpts > npoints)
    return 0;
  indices[0] = ptr;
  ptr = (ptr == npoints - 1 ? 0 : ptr + 1);
  i = 1;
  while (i < minpts) {
    int index = rand() % npoints;
    while (index) {
      ptr = (ptr == npoints - 1 ? 0 : ptr + 1);
      for (j = 0; j < i; ++j) {
        if (indices[j] == ptr)
          break;
      }
      if (j == i)
        index--;
    }
    indices[i++] = ptr;
  }
  return 1;
}

int ransac_(double *matched_points,
            int npoints,
            int *number_of_inliers,
            int *best_inlier_mask,
            double *bestH,
            const int minpts,
            const int paramdim,
            isDegenerateType isDegenerate,
            normalizeType normalize,
            denormalizeType denormalize,
            findTransformationType findTransformation,
            projectPointsType projectPoints) {

  static const double INLIER_THRESHOLD_NORMALIZED = 0.1;
  static const double INLIER_THRESHOLD_UNNORMALIZED = 1.0;
  static const double PROBABILITY_REQUIRED = 0.9;
  static const double EPS = 1e-12;
  static const int MIN_TRIALS = 20;

  const double inlier_threshold = (normalize && denormalize ?
                                   INLIER_THRESHOLD_NORMALIZED :
                                   INLIER_THRESHOLD_UNNORMALIZED);
  int N = 10000, trial_count = 0;
  int i;
  int ret_val = 0;

  int max_inliers = 0;
  double best_variance = 0.0;
  double H[MAX_PARAMDIM];
  double points1[2 * MAX_MINPTS];
  double points2[2 * MAX_MINPTS];
  int indices[MAX_MINPTS];

  double *best_inlier_set1;
  double *best_inlier_set2;
  double *inlier_set1;
  double *inlier_set2;
  double *corners1;
  double *corners2;
  double *image1_coord;
  double *image2_coord;
  int *inlier_mask;

  double *cnp1, *cnp2;
  double T1[9], T2[9];

  // srand((unsigned)time(NULL)) ;
  // better to make this deterministic for a given sequence for ease of testing
  srand(npoints);

  *number_of_inliers = 0;
  if (npoints < minpts * MINPTS_MULTIPLIER) {
    printf("Cannot find motion with %d matches\n", npoints);
    return 1;
  }

  best_inlier_set1 = (double *)malloc(sizeof(double) * npoints * 2);
  best_inlier_set2 = (double *)malloc(sizeof(double) * npoints * 2);
  inlier_set1 = (double *)malloc(sizeof(double) * npoints * 2);
  inlier_set2 = (double *)malloc(sizeof(double) * npoints * 2);
  corners1 = (double *)malloc(sizeof(double) * npoints * 2);
  corners2 = (double *)malloc(sizeof(double) * npoints * 2);
  image1_coord = (double *)malloc(sizeof(double) * npoints * 2);
  image2_coord = (double *)malloc(sizeof(double) * npoints * 2);
  inlier_mask = (int*)malloc(sizeof(int) * npoints);

  for(cnp1 = corners1, cnp2 = corners2, i = 0; i < npoints; ++i) {
    *(cnp1++) = *(matched_points++);
    *(cnp1++) = *(matched_points++);
    *(cnp2++) = *(matched_points++);
    *(cnp2++) = *(matched_points++);
  }
  matched_points -= 4 * npoints;

  if (normalize && denormalize) {
    normalize(corners1, npoints, T1);
    normalize(corners2, npoints, T2);
  }

  while (N > trial_count) {
    int num_inliers = 0;
    double sum_distance = 0.0;
    double sum_distance_squared = 0.0;

    int degenerate = 1;
    int num_degenerate_iter = 0;
    while (degenerate) {
      num_degenerate_iter++;
      if (!get_rand_indices(npoints, minpts, indices)) {
        ret_val = 1;
        goto finish_ransac;
      }
      i = 0;
      while (i < minpts) {
        int index = indices[i];
        // add to list
        points1[i*2] = corners1[index*2];
        points1[i*2+1] = corners1[index*2+1];
        points2[i*2] = corners2[index*2];
        points2[i*2+1] = corners2[index*2+1];
        i++;
      }
      degenerate = isDegenerate(points1);
      if (num_degenerate_iter > MAX_DEGENERATE_ITER) {
        ret_val = 1;
        goto finish_ransac;
      }
    }

    if (findTransformation(minpts, points1, points2, H)) {
      trial_count++;
      continue;
    }

    projectPoints(H, corners1, image1_coord, npoints, 2, 2);

    for( i = 0; i < npoints; ++i ) {
      double dx = image1_coord[i*2] - corners2[i*2];
      double dy = image1_coord[i*2 + 1] - corners2[i*2 + 1];
      double distance = sqrt(dx*dx + dy*dy);

      inlier_mask[i] = distance < inlier_threshold;
      if (inlier_mask[i]) {
        inlier_set1[num_inliers*2] =            corners1[i*2];
        inlier_set1[num_inliers*2 + 1] =        corners1[i*2 + 1];
        inlier_set2[num_inliers*2] =            corners2[i*2];
        inlier_set2[num_inliers*2 + 1] =        corners2[i*2 + 1];
        num_inliers++;
        sum_distance += distance;
        sum_distance_squared += distance*distance;
      }
    }

    if (num_inliers >= max_inliers) {
      double mean_distance = sum_distance / ((double)num_inliers);
      double variance = sum_distance_squared / ((double)num_inliers - 1.0) -
          mean_distance * mean_distance *
          ((double)num_inliers) / ((double)num_inliers - 1.0);
      if ((num_inliers > max_inliers) ||
          (num_inliers==max_inliers && variance < best_variance)) {
        best_variance = variance;
        max_inliers = num_inliers;
        memcpy(bestH, H, paramdim * sizeof(double));
        memcpy(best_inlier_set1, inlier_set1, num_inliers*2 * sizeof(double));
        memcpy(best_inlier_set2, inlier_set2, num_inliers*2 * sizeof(double));
        memcpy(best_inlier_mask, inlier_mask, npoints * sizeof(int));

        if (num_inliers > 0) {
          double fracinliers =  (double)num_inliers/(double)npoints;
          double pNoOutliers = 1 - pow(fracinliers, minpts);
          int temp;
          pNoOutliers = fmax(EPS, pNoOutliers);
          pNoOutliers = fmin(1 - EPS, pNoOutliers);
          temp = (int)(log(1.0 - PROBABILITY_REQUIRED)/log(pNoOutliers));
          if (temp > 0 && temp < N) {
            N = IMAX(temp, MIN_TRIALS);
          }
        }
      }
    }
    trial_count++;
  }
  // printf("Number of trials = %d\n", trial_count);
  findTransformation(max_inliers, best_inlier_set1, best_inlier_set2, bestH);
  if (normalize && denormalize) {
    denormalize(bestH, T1, T2);
  }
  *number_of_inliers = max_inliers;
  /*
  printf("Error score (all) = %g\n",
         compute_error(projectPoints, matched_points, 4,
                       matched_points + 2, 4,
                       npoints, bestH, NULL));
  printf("Error score (inliers) = %g\n",
         compute_error(projectPoints, matched_points, 4,
                       matched_points + 2, 4,
                       npoints, bestH,
                       best_inlier_mask));
  */
finish_ransac:
  free(best_inlier_set1);
  free(best_inlier_set2);
  free(inlier_set1);
  free(inlier_set2);
  free(corners1);
  free(corners2);
  free(image1_coord);
  free(image2_coord);
  free(inlier_mask);
  return ret_val;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////

static void normalizeHomography(double *pts, int n, double *T) {
  // Assume the points are 2d coordinates with scale = 1
  double *p = pts;
  double mean[2] = {0, 0};
  double msqe = 0;
  double scale;
  int i;
  for (i = 0; i < n; ++i, p+=2) {
    mean[0] += p[0];
    mean[1] += p[1];
  }
  mean[0] /= n;
  mean[1] /= n;
  for (p = pts, i = 0; i < n; ++i, p+=2) {
    p[0] -= mean[0];
    p[1] -= mean[1];
    msqe += sqrt(p[0] * p[0] + p[1] * p[1]);
  }
  msqe /= n;
  scale = sqrt(2)/msqe;
  T[0] = scale;
  T[1] = 0;
  T[2] = -scale * mean[0];
  T[3] = 0;
  T[4] = scale;
  T[5] = -scale * mean[1];
  T[6] = 0;
  T[7] = 0;
  T[8] = 1;
  for (p = pts, i = 0; i < n; ++i, p+=2) {
    p[0] *= scale;
    p[1] *= scale;
  }
}

static void invnormalize_mat(double *T, double *iT) {
  double is = 1.0/T[0];
  double m0 = -T[2]*is;
  double m1 = -T[5]*is;
  iT[0] = is;
  iT[1] = 0;
  iT[2] = m0;
  iT[3] = 0;
  iT[4] = is;
  iT[5] = m1;
  iT[6] = 0;
  iT[7] = 0;
  iT[8] = 1;
}

static void denormalizeHomography(double *H, double *T1, double *T2) {
  double iT2[9];
  double H2[9];
  invnormalize_mat(T2, iT2);
  MultiplyMat(H, T1, H2, 3, 3, 3);
  MultiplyMat(iT2, H2, H, 3, 3, 3);
}

static void denormalizeAffine(double *H, double *T1, double *T2) {
  double Ha[MAX_PARAMDIM];
  Ha[0] = H[0];
  Ha[1] = H[1];
  Ha[2] = H[4];
  Ha[3] = H[2];
  Ha[4] = H[3];
  Ha[5] = H[5];
  Ha[6] = Ha[7] = 0;
  Ha[8] = 1;
  denormalizeHomography(Ha, T1, T2);
  H[0] = Ha[0];
  H[1] = Ha[1];
  H[2] = Ha[3];
  H[3] = Ha[4];
  H[4] = Ha[2];
  H[5] = Ha[5];
}

static void denormalizeRotZoom(double *H, double *T1, double *T2) {
  double Ha[MAX_PARAMDIM];
  Ha[0] = H[0];
  Ha[1] = H[1];
  Ha[2] = H[2];
  Ha[3] = -H[1];
  Ha[4] = H[0];
  Ha[5] = H[3];
  Ha[6] = Ha[7] = 0;
  Ha[8] = 1;
  denormalizeHomography(Ha, T1, T2);
  H[0] = Ha[0];
  H[1] = Ha[1];
  H[2] = Ha[2];
  H[3] = Ha[5];
}

static void denormalizeTranslation(double *H, double *T1, double *T2) {
  double Ha[MAX_PARAMDIM];
  Ha[0] = 1;
  Ha[1] = 0;
  Ha[2] = H[0];
  Ha[3] = 0;
  Ha[4] = 1;
  Ha[5] = H[1];
  Ha[6] = Ha[7] = 0;
  Ha[8] = 1;
  denormalizeHomography(Ha, T1, T2);
  H[0] = Ha[2];
  H[1] = Ha[5];
}

static int is_collinear3(double *p1, double *p2, double *p3) {
  static const double collinear_eps = 1e-3;
  const double v = (p2[0] - p1[0]) * (p3[1] - p1[1]) -
                   (p2[1] - p1[1]) * (p3[0] - p1[0]);
  return fabs(v) < collinear_eps;
}

static int isDegenerateTranslation(double *p) {
  return (p[0] - p[2]) * (p[0] - p[2]) + (p[1] - p[3]) * (p[1] - p[3]) <= 2;
}

static int isDegenerateAffine(double *p) {
  return is_collinear3(p, p + 2, p + 4);
}

static int isDegenerateHomography(double *p) {
  return is_collinear3(p, p + 2, p + 4) ||
         is_collinear3(p, p + 2, p + 6) ||
         is_collinear3(p, p + 4, p + 6) ||
         is_collinear3(p + 2, p + 4, p + 6);
}

int findTranslation(const int np, double *pts1, double *pts2, double *mat) {
  int i;
  double sx, sy, dx, dy;
  double sumx, sumy;

  double T1[9], T2[9];
  normalizeHomography(pts1, np, T1);
  normalizeHomography(pts2, np, T2);

  sumx = 0;
  sumy = 0;
  for (i = 0; i < np; ++i) {
    dx = *(pts2++);
    dy = *(pts2++);
    sx = *(pts1++);
    sy = *(pts1++);

    sumx += dx - sx;
    sumy += dy - sy;
  }
  mat[0] = sumx / np;
  mat[1] = sumy / np;
  denormalizeTranslation(mat, T1, T2);
  return 0;
}

int findRotZoom(const int np, double *pts1, double *pts2, double *mat) {
  const int np2 = np * 2;
  double *a = (double *)malloc(sizeof(double) * np2 * 9);
  double *b = a + np2 * 4;
  double *temp = b + np2;
  int i;
  double sx, sy, dx, dy;

  double T1[9], T2[9];
  normalizeHomography(pts1, np, T1);
  normalizeHomography(pts2, np, T2);

  for (i = 0; i < np; ++i) {
    dx = *(pts2++);
    dy = *(pts2++);
    sx = *(pts1++);
    sy = *(pts1++);

    a[i * 2 * 4 + 0] = sx;
    a[i * 2 * 4 + 1] = sy;
    a[i * 2 * 4 + 2] = 1;
    a[i * 2 * 4 + 3] = 0;
    a[(i * 2 + 1) * 4 + 0] = sy;
    a[(i * 2 + 1) * 4 + 1] = -sx;
    a[(i * 2 + 1) * 4 + 2] = 0;
    a[(i * 2 + 1) * 4 + 3] = 1;

    b[2 * i] = dx;
    b[2 * i + 1] = dy;
  }
  if (PseudoInverse(temp, a, np2, 4)){
    free(a);
    return 1;
  }
  MultiplyMat(temp, b, mat, 4, np2, 1);
  denormalizeRotZoom(mat, T1, T2);
  free(a);
  return 0;
}

int findAffine(const int np, double *pts1, double *pts2, double *mat) {
  const int np2 = np * 2;
  double *a = (double *)malloc(sizeof(double) * np2 * 13);
  double *b = a + np2 * 6;
  double *temp = b + np2;
  int i;
  double sx, sy, dx, dy;

  double T1[9], T2[9];
  normalizeHomography(pts1, np, T1);
  normalizeHomography(pts2, np, T2);

  for (i = 0; i < np; ++i) {
    dx = *(pts2++);
    dy = *(pts2++);
    sx = *(pts1++);
    sy = *(pts1++);

    a[i * 2 * 6 + 0] = sx;
    a[i * 2 * 6 + 1] = sy;
    a[i * 2 * 6 + 2] = 0;
    a[i * 2 * 6 + 3] = 0;
    a[i * 2 * 6 + 4] = 1;
    a[i * 2 * 6 + 5] = 0;
    a[(i * 2 + 1) * 6 + 0] = 0;
    a[(i * 2 + 1) * 6 + 1] = 0;
    a[(i * 2 + 1) * 6 + 2] = sx;
    a[(i * 2 + 1) * 6 + 3] = sy;
    a[(i * 2 + 1) * 6 + 4] = 0;
    a[(i * 2 + 1) * 6 + 5] = 1;

    b[2 * i] = dx;
    b[2 * i + 1] = dy;
  }
  if (PseudoInverse(temp, a, np2, 6)){
    free(a);
    return 1;
  }
  MultiplyMat(temp, b, mat, 6, np2, 1);
  denormalizeAffine(mat, T1, T2);
  free(a);
  return 0;
}

int findHomography(const int np, double *pts1, double *pts2, double *mat) {
  // Implemented from Peter Kovesi's normalized implementation
  const int np3 = np * 3;
  double *a = (double *)malloc(sizeof(double) * np3 * 18);
  double *U = a + np3 * 9;
  double S[9], V[9 * 9];
  int i, mini;
  double sx, sy, dx, dy;
  double T1[9], T2[9];

  normalizeHomography(pts1, np, T1);
  normalizeHomography(pts2, np, T2);

  for (i = 0; i < np; ++i) {
    dx = *(pts2++);
    dy = *(pts2++);
    sx = *(pts1++);
    sy = *(pts1++);

    a[i * 3 * 9 + 0] = a[i * 3 * 9 + 1] = a[i * 3 * 9 + 2] = 0;
    a[i * 3 * 9 + 3] = -sx;
    a[i * 3 * 9 + 4] = -sy;
    a[i * 3 * 9 + 5] = -1;
    a[i * 3 * 9 + 6] = dy * sx;
    a[i * 3 * 9 + 7] = dy * sy;
    a[i * 3 * 9 + 8] = dy;

    a[(i * 3 + 1) * 9 + 0] = sx;
    a[(i * 3 + 1) * 9 + 1] = sy;
    a[(i * 3 + 1) * 9 + 2] = 1;
    a[(i * 3 + 1) * 9 + 3] = a[(i * 3 + 1) * 9 + 4] =
        a[(i * 3 + 1) * 9 + 5] = 0;
    a[(i * 3 + 1) * 9 + 6] = -dx * sx;
    a[(i * 3 + 1) * 9 + 7] = -dx * sy;
    a[(i * 3 + 1) * 9 + 8] = -dx;

    a[(i * 3 + 2) * 9 + 0] = -dy * sx;
    a[(i * 3 + 2) * 9 + 1] = -dy * sy;
    a[(i * 3 + 2) * 9 + 2] = -dy;
    a[(i * 3 + 2) * 9 + 3] = dx * sx;
    a[(i * 3 + 2) * 9 + 4] = dx * sy;
    a[(i * 3 + 2) * 9 + 5] = dx;
    a[(i * 3 + 2) * 9 + 6] = a[(i * 3 + 2) * 9 + 7] =
        a[(i * 3 + 2) * 9 + 8] = 0;
  }

  if (SVD(U, S, V, a, np3, 9)) {
    free(a);
    return 1;
  } else {
    double minS = 1e12;
    mini = -1;
    for (i = 0; i < 9; ++i) {
      if (S[i] < minS) {
        minS = S[i];
        mini = i;
      }
    }
  }

  for (i = 0; i < 9; i++) mat[i] = V[i * 9 + mini];
  denormalizeHomography(mat, T1, T2);
  free(a);
  if (mat[8] == 0.0) {
    return 1;
  }
  return 0;
}

int findHomographyScale1(const int np, double *pts1, double *pts2,
                         double *mat) {
  // This implementation assumes h33 = 1, but does not seem to give good results
  const int np2 = np * 2;
  double *a = (double *)malloc(sizeof(double) * np2 * 17);
  double *b = a + np2 * 8;
  double *temp = b + np2;
  int i, j;
  double sx, sy, dx, dy;
  double T1[9], T2[9];

  normalizeHomography(pts1, np, T1);
  normalizeHomography(pts2, np, T2);

  for (i = 0, j = np; i < np; ++i, ++j) {
    dx = *(pts2++);
    dy = *(pts2++);
    sx = *(pts1++);
    sy = *(pts1++);
    a[i * 8 + 0] = a[j * 8 + 3] = sx;
    a[i * 8 + 1] = a[j * 8 + 4] = sy;
    a[i * 8 + 2] = a[j * 8 + 5] = 1;
    a[i * 8 + 3] = a[i * 8 + 4] = a[i * 8 + 5] =
        a[j * 8 + 0] = a[j * 8 + 1] = a[j * 8 + 2] = 0;
    a[i * 8 + 6] = -dx * sx;
    a[i * 8 + 7] = -dx * sy;
    a[j * 8 + 6] = -dy * sx;
    a[j * 8 + 7] = -dy * sy;
    b[i] = dx;
    b[j] = dy;
  }

  if (PseudoInverse(temp, a, np2, 8)) {
    free(a);
    return 1;
  }
  MultiplyMat(temp, b, &*mat, 8, np2, 1);
  mat[8] = 1;

  denormalizeHomography(mat, T1, T2);
  free(a);
  return 0;
}

int ransacTranslation(double *matched_points, int npoints,
                      int *number_of_inliers,
                      int *best_inlier_mask, double *bestH) {
  return ransac_(matched_points,
                 npoints,
                 number_of_inliers,
                 best_inlier_mask,
                 bestH,
                 3,
                 2,
                 isDegenerateTranslation,
                 NULL,  // normalizeHomography,
                 NULL,  // denormalizeRotZoom,
                 findTranslation,
                 projectPointsTranslation);
}

int ransacRotZoom(double *matched_points, int npoints,
                  int *number_of_inliers,
                  int *best_inlier_mask, double *bestH) {
  return ransac_(matched_points,
                 npoints,
                 number_of_inliers,
                 best_inlier_mask,
                 bestH,
                 3,
                 4,
                 isDegenerateAffine,
                 NULL,  // normalizeHomography,
                 NULL,  // denormalizeRotZoom,
                 findRotZoom,
                 projectPointsRotZoom);
}

int ransacAffine(double *matched_points, int npoints,
                 int *number_of_inliers,
                 int *best_inlier_mask, double *bestH) {
  return ransac_(matched_points,
                 npoints,
                 number_of_inliers,
                 best_inlier_mask,
                 bestH,
                 3,
                 6,
                 isDegenerateAffine,
                 NULL, // normalizeHomography,
                 NULL, // denormalizeAffine,
                 findAffine,
                 projectPointsAffine);
}

int ransacHomography(double *matched_points, int npoints,
                     int *number_of_inliers,
                     int *best_inlier_mask, double *bestH) {
  int result =
      ransac_(matched_points,
              npoints,
              number_of_inliers,
              best_inlier_mask,
              bestH,
              4,
              8,
              isDegenerateHomography,
              NULL, // normalizeHomography,
              NULL, // denormalizeHomography,
              findHomography,
              projectPointsHomography);
  if (!result) {
    // normalize so that H33 = 1
    int i;
    double m = 1.0 / bestH[8];
    for (i = 0; i < 8; ++i)
      bestH[i] *= m;
    bestH[8] = 1.0;
  }
  return result;
}
