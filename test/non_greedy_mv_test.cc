/*
 *  Copyright (c) 2019 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>
#include "third_party/googletest/src/include/gtest/gtest.h"
#include "vp9/encoder/vp9_non_greedy_mv.h"

namespace {

static void read_in_mf(const char *filename, int *rows_ptr, int *cols_ptr,
                       MV **buffer_ptr) {
  FILE *input = fopen(filename, "rb");
  int row, col;
  int idx;

  ASSERT_NE(input, nullptr) << "Cannot open file: " << filename << std::endl;

  fscanf(input, "%d,%d\n", rows_ptr, cols_ptr);

  *buffer_ptr = (MV *)malloc((*rows_ptr) * (*cols_ptr) * sizeof(MV));

  for (idx = 0; idx < (*rows_ptr) * (*cols_ptr); ++idx) {
    fscanf(input, "%d,%d;", &row, &col);
    (*buffer_ptr)[idx].row = row;
    (*buffer_ptr)[idx].col = col;
  }
  fclose(input);
}

static void read_in_local_var(const char *filename, int *rows_ptr,
                              int *cols_ptr,
                              int (**M_ptr)[MF_LOCAL_STRUCTURE_SIZE]) {
  FILE *input = fopen(filename, "rb");
  int M00, M01, M10, M11;
  int idx;
  int int_type;

  ASSERT_NE(input, nullptr) << "Cannot open file: " << filename << std::endl;

  fscanf(input, "%d,%d\n", rows_ptr, cols_ptr);

  *M_ptr = (int(*)[MF_LOCAL_STRUCTURE_SIZE])malloc(
      (*rows_ptr) * (*cols_ptr) * MF_LOCAL_STRUCTURE_SIZE * sizeof(int_type));

  for (idx = 0; idx < (*rows_ptr) * (*cols_ptr); ++idx) {
    fscanf(input, "%d,%d,%d,%d;", &M00, &M01, &M10, &M11);
    (*M_ptr)[idx][0] = M00;
    (*M_ptr)[idx][1] = M01;
    (*M_ptr)[idx][2] = M10;
    (*M_ptr)[idx][3] = M11;
  }
  fclose(input);
}

static void compare_mf(const MV *mf1, const MV *mf2, int rows, int cols,
                       float *mean_ptr, float *std_ptr) {
  float float_type;
  float *diffs = (float *)malloc(rows * cols * sizeof(float_type));
  int idx;
  float accu = 0.0f;
  for (idx = 0; idx < rows * cols; ++idx) {
    MV mv1 = mf1[idx];
    MV mv2 = mf2[idx];
    float row_diff2 = (float)((mv1.row - mv2.row) * (mv1.row - mv2.row));
    float col_diff2 = (float)((mv1.col - mv2.col) * (mv1.col - mv2.col));
    diffs[idx] = sqrt(row_diff2 + col_diff2);
    accu += diffs[idx];
  }
  *mean_ptr = accu / rows / cols;
  *std_ptr = 0;
  for (idx = 0; idx < rows * cols; ++idx) {
    *std_ptr += (diffs[idx] - (*mean_ptr)) * (diffs[idx] - (*mean_ptr));
  }
  *std_ptr = sqrt(*std_ptr / rows / cols);
  free(diffs);
}

TEST(non_greedy_mv, smooth_mf) {
  const char *search_mf_file = "non_greedy_mv_test_files/exhaust_32x32.txt";
  const char *local_var_file = "non_greedy_mv_test_files/localVar_32x32.txt";
  const char *estimation_file = "non_greedy_mv_test_files/estimation_32x32.txt";
  const char *ground_truth_file =
      "non_greedy_mv_test_files/ground_truth_32x32.txt";
  BLOCK_SIZE bsize = BLOCK_32X32;
  MV *search_mf = NULL;
  MV *smooth_mf = NULL;
  MV *estimation = NULL;
  MV *ground_truth = NULL;
  int(*local_var)[MF_LOCAL_STRUCTURE_SIZE] = NULL;
  int rows = 0, cols = 0;

  int alpha = 100, max_iter = 100;

  read_in_mf(search_mf_file, &rows, &cols, &search_mf);
  read_in_local_var(local_var_file, &rows, &cols, &local_var);
  read_in_mf(estimation_file, &rows, &cols, &estimation);
  read_in_mf(ground_truth_file, &rows, &cols, &ground_truth);

  float sm_mean, sm_std;
  float est_mean, est_std;

  smooth_mf = (MV *)malloc(rows * cols * sizeof(MV));
  vp9_get_smooth_motion_field(search_mf, local_var, rows, cols, bsize, alpha,
                              max_iter, smooth_mf);

  compare_mf(smooth_mf, ground_truth, rows, cols, &sm_mean, &sm_std);
  compare_mf(smooth_mf, estimation, rows, cols, &est_mean, &est_std);

  EXPECT_LE(sm_mean, 3);
  EXPECT_LE(est_mean, 2);

  free(search_mf);
  free(local_var);
  free(estimation);
  free(ground_truth);
  free(smooth_mf);
}
}  // namespace
