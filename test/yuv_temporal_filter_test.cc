/*
 *  Copyright (c) 2019 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./vp9_rtcd.h"
#include "test/acm_random.h"
#include "test/buffer.h"
#include "test/register_state_check.h"
#include "vpx_ports/vpx_timer.h"

namespace {

using ::libvpx_test::ACMRandom;
using ::libvpx_test::Buffer;

typedef void (*YUVTemporalFilterFunc)(
    const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre,
    int y_pre_stride, const uint8_t *u_src, const uint8_t *v_src,
    int uv_src_stride, const uint8_t *u_pre, const uint8_t *v_pre,
    int uv_pre_stride, unsigned int block_width, unsigned int block_height,
    int ss_x, int ss_y, int strength, const int *const blk_fw, int use_32x32,
    uint32_t *y_accumulator, uint16_t *y_count, uint32_t *u_accumulator,
    uint16_t *u_count, uint32_t *v_accumulator, uint16_t *v_count);

int GetFilterWeight(unsigned int row, unsigned int col,
                    unsigned int block_height, unsigned int block_width,
                    const int *const blk_fw, int use_32x32) {
  if (use_32x32) {
    return blk_fw[0];
  }

  return blk_fw[2 * (row >= block_height / 2) + (col >= block_width / 2)];
}

int GetModIndex(int sum_dist, int index, int rounding, int strength,
                int filter_weight) {
  int mod = (sum_dist * 3) / index;
  mod += rounding;
  mod >>= strength;

  mod = VPXMIN(16, mod);

  mod = 16 - mod;
  mod *= filter_weight;

  return mod;
}

void ApplyReferenceFilter(
    const Buffer<uint8_t> &y_src, const Buffer<uint8_t> &y_pre,
    const Buffer<uint8_t> &u_src, const Buffer<uint8_t> &v_src,
    const Buffer<uint8_t> &u_pre, const Buffer<uint8_t> &v_pre,
    unsigned int block_width, unsigned int block_height, int ss_x, int ss_y,
    int strength, const int *const blk_fw, int use_32x32,
    Buffer<uint32_t> *y_accumulator, Buffer<uint16_t> *y_count,
    Buffer<uint32_t> *u_accumulator, Buffer<uint16_t> *u_count,
    Buffer<uint32_t> *v_accumulator, Buffer<uint16_t> *v_count) {
  // blk_fw means block_filter_weight
  // Set up buffer to store squared_diffs
  Buffer<int> y_dif = Buffer<int>(block_width, block_height, 0);
  const int uv_block_width = block_width >> ss_x;
  const int uv_block_height = block_height >> ss_y;
  Buffer<int> u_dif = Buffer<int>(uv_block_width, uv_block_height, 0);
  Buffer<int> v_dif = Buffer<int>(uv_block_width, uv_block_height, 0);
  ASSERT_TRUE(y_dif.Init());
  ASSERT_TRUE(u_dif.Init());
  ASSERT_TRUE(v_dif.Init());
  y_dif.Set(0);
  u_dif.Set(0);
  v_dif.Set(0);

  // How many bits to we want to round
  ASSERT_GE(strength, 0);
  ASSERT_LE(strength, 6);
  int rounding = 0;
  if (strength > 0) {
    rounding = 1 << (strength - 1);
  }

  // Check that the buffers are valid
  ASSERT_TRUE(y_src.TopLeftPixel() != NULL);
  ASSERT_TRUE(y_pre.TopLeftPixel() != NULL);
  ASSERT_TRUE(y_dif.TopLeftPixel() != NULL);
  ASSERT_TRUE(u_src.TopLeftPixel() != NULL);
  ASSERT_TRUE(u_pre.TopLeftPixel() != NULL);
  ASSERT_TRUE(u_dif.TopLeftPixel() != NULL);
  ASSERT_TRUE(v_src.TopLeftPixel() != NULL);
  ASSERT_TRUE(v_pre.TopLeftPixel() != NULL);
  ASSERT_TRUE(v_dif.TopLeftPixel() != NULL);

  // Get the square diffs
  for (int row = 0; row < static_cast<int>(block_height); row++) {
    for (int col = 0; col < static_cast<int>(block_width); col++) {
      const int diff = y_src.TopLeftPixel()[row * y_src.stride() + col] -
                       y_pre.TopLeftPixel()[row * y_pre.stride() + col];
      y_dif.TopLeftPixel()[row * y_dif.stride() + col] = diff * diff;
    }
  }

  for (int row = 0; row < uv_block_height; row++) {
    for (int col = 0; col < uv_block_width; col++) {
      const int u_diff = u_src.TopLeftPixel()[row * u_src.stride() + col] -
                         u_pre.TopLeftPixel()[row * u_pre.stride() + col];
      const int v_diff = v_src.TopLeftPixel()[row * v_src.stride() + col] -
                         v_pre.TopLeftPixel()[row * v_pre.stride() + col];
      u_dif.TopLeftPixel()[row * u_dif.stride() + col] = u_diff * u_diff;
      v_dif.TopLeftPixel()[row * v_dif.stride() + col] = v_diff * v_diff;
    }
  }

  // Apply the filter
  for (int row = 0; row < static_cast<int>(block_height); row++) {
    for (int col = 0; col < static_cast<int>(block_width); col++) {
      const int uv_r = row >> ss_y;
      const int uv_c = col >> ss_x;
      const int filter_weight = GetFilterWeight(row, col, block_height,
                                                block_width, blk_fw, use_32x32);

      // First we get the modifier for the current y pixel
      const int y_pixel = y_pre.TopLeftPixel()[row * y_pre.stride() + col];
      int y_num_used = 0;
      int y_mod = 0;

      // Sum the neighboring 3x3 y pixels
      for (int row_step = -1; row_step <= 1; row_step++) {
        for (int col_step = -1; col_step <= 1; col_step++) {
          const int sub_row = row + row_step;
          const int sub_col = col + col_step;

          if (sub_row >= 0 && sub_row < static_cast<int>(block_height) &&
              sub_col >= 0 && sub_col < static_cast<int>(block_width)) {
            y_mod += y_dif.TopLeftPixel()[sub_row * y_dif.stride() + sub_col];
            y_num_used++;
          }
        }
      }

      ASSERT_GE(y_num_used, 0);

      // Sum the corresponding uv pixels to the current y modifier
      // Note we are rounding down instead of rounding to the nearest pixel.
      y_mod += u_dif.TopLeftPixel()[uv_r * uv_block_width + uv_c];
      y_mod += v_dif.TopLeftPixel()[uv_r * uv_block_width + uv_c];

      y_num_used += 2;

      // Set the modifier
      y_mod = GetModIndex(y_mod, y_num_used, rounding, strength, filter_weight);

      // Accumulate the result
      y_count->TopLeftPixel()[row * y_count->stride() + col] += y_mod;
      y_accumulator->TopLeftPixel()[row * y_accumulator->stride() + col] +=
          y_mod * y_pixel;

      // Get the modifier for chroma components
      if (!(row & ss_y) && !(col & ss_x)) {
        const int u_pixel = u_pre.TopLeftPixel()[uv_r * u_pre.stride() + uv_c];
        const int v_pixel = v_pre.TopLeftPixel()[uv_r * v_pre.stride() + uv_c];

        int uv_num_used = 0;
        int u_mod = 0, v_mod = 0;

        // Sum the neighboring 3x3 chromal pixels to the chroma modifier
        for (int row_step = -1; row_step <= 1; row_step++) {
          for (int col_step = -1; col_step <= 1; col_step++) {
            const int sub_row = uv_r + row_step;
            const int sub_col = uv_c + col_step;

            if (sub_row >= 0 && sub_row < uv_block_height && sub_col >= 0 &&
                sub_col < uv_block_width) {
              u_mod += u_dif.TopLeftPixel()[sub_row * uv_block_width + sub_col];
              v_mod += v_dif.TopLeftPixel()[sub_row * uv_block_width + sub_col];
              uv_num_used++;
            }
          }
        }

        ASSERT_GT(uv_num_used, 0);

        // Sum all the luma pixels associated with the current luma pixel
        for (int row_step = 0; row_step < 1 + ss_y; row_step++) {
          for (int col_step = 0; col_step < 1 + ss_x; col_step++) {
            const int sub_row = (uv_r << ss_y) + row_step;
            const int sub_col = (uv_c << ss_x) + col_step;
            const int y_diff =
                y_dif.TopLeftPixel()[sub_row * y_dif.stride() + sub_col];

            u_mod += y_diff;
            v_mod += y_diff;
            uv_num_used++;
          }
        }

        // Set the modifier
        u_mod =
            GetModIndex(u_mod, uv_num_used, rounding, strength, filter_weight);
        v_mod =
            GetModIndex(v_mod, uv_num_used, rounding, strength, filter_weight);

        // Accumulate the result
        u_count->TopLeftPixel()[uv_r * u_count->stride() + uv_c] += u_mod;
        u_accumulator->TopLeftPixel()[uv_r * u_accumulator->stride() + uv_c] +=
            u_mod * u_pixel;
        v_count->TopLeftPixel()[uv_r * u_count->stride() + uv_c] += v_mod;
        v_accumulator->TopLeftPixel()[uv_r * v_accumulator->stride() + uv_c] +=
            v_mod * v_pixel;
      }
    }
  }
}

class YUVTemporalFilterTest
    : public ::testing::TestWithParam<YUVTemporalFilterFunc> {
 public:
  virtual void SetUp() {
    filter_func_ = GetParam();
    rnd_.Reset(ACMRandom::DeterministicSeed());
  }

 protected:
  YUVTemporalFilterFunc filter_func_;
  ACMRandom rnd_;
};

TEST_P(YUVTemporalFilterTest, Use_32X32) {
  const int width = 32, height = 32;
  Buffer<uint8_t> y_src = Buffer<uint8_t>(width, height, 8);
  Buffer<uint8_t> y_pre = Buffer<uint8_t>(width, height, 0);
  Buffer<uint16_t> y_count_ref = Buffer<uint16_t>(width, height, 0);
  Buffer<uint32_t> y_accum_ref = Buffer<uint32_t>(width, height, 0);
  Buffer<uint16_t> y_count_tst = Buffer<uint16_t>(width, height, 0);
  Buffer<uint32_t> y_accum_tst = Buffer<uint32_t>(width, height, 0);
  ASSERT_TRUE(y_src.Init());
  ASSERT_TRUE(y_pre.Init());
  ASSERT_TRUE(y_count_ref.Init());
  ASSERT_TRUE(y_accum_ref.Init());
  ASSERT_TRUE(y_count_tst.Init());
  ASSERT_TRUE(y_accum_tst.Init());

  const int use_32x32 = 1;

  for (int ss_x = 0; ss_x <= 1; ss_x++) {
    for (int ss_y = 0; ss_y <= 1; ss_y++) {
      for (int filter_strength = 0; filter_strength <= 6;
           filter_strength += 2) {
        for (int filter_weight = 0; filter_weight <= 2; filter_weight++) {
          const int uv_width = width >> ss_x, uv_height = height >> ss_y;
          Buffer<uint8_t> u_src = Buffer<uint8_t>(uv_width, uv_height, 8);
          Buffer<uint8_t> u_pre = Buffer<uint8_t>(uv_width, uv_height, 0);
          Buffer<uint16_t> u_count_ref =
              Buffer<uint16_t>(uv_width, uv_height, 0);
          Buffer<uint32_t> u_accum_ref =
              Buffer<uint32_t>(uv_width, uv_height, 0);
          Buffer<uint16_t> u_count_tst =
              Buffer<uint16_t>(uv_width, uv_height, 0);
          Buffer<uint32_t> u_accum_tst =
              Buffer<uint32_t>(uv_width, uv_height, 0);
          ASSERT_TRUE(u_src.Init());
          ASSERT_TRUE(u_pre.Init());
          ASSERT_TRUE(u_count_ref.Init());
          ASSERT_TRUE(u_accum_ref.Init());
          ASSERT_TRUE(u_count_tst.Init());
          ASSERT_TRUE(u_accum_tst.Init());
          Buffer<uint8_t> v_src = Buffer<uint8_t>(uv_width, uv_height, 8);
          Buffer<uint8_t> v_pre = Buffer<uint8_t>(uv_width, uv_height, 0);
          Buffer<uint16_t> v_count_ref =
              Buffer<uint16_t>(uv_width, uv_height, 0);
          Buffer<uint32_t> v_accum_ref =
              Buffer<uint32_t>(uv_width, uv_height, 0);
          Buffer<uint16_t> v_count_tst =
              Buffer<uint16_t>(uv_width, uv_height, 0);
          Buffer<uint32_t> v_accum_tst =
              Buffer<uint32_t>(uv_width, uv_height, 0);
          ASSERT_TRUE(v_src.Init());
          ASSERT_TRUE(v_pre.Init());
          ASSERT_TRUE(v_count_ref.Init());
          ASSERT_TRUE(v_accum_ref.Init());
          ASSERT_TRUE(v_count_tst.Init());
          ASSERT_TRUE(v_accum_tst.Init());

          // The difference between the buffers must be small to pass the
          // threshold to apply the filter.
          y_src.Set(&rnd_, 0, 7);
          y_pre.Set(&rnd_, 0, 7);
          u_src.Set(&rnd_, 0, 7);
          u_pre.Set(&rnd_, 0, 7);
          v_src.Set(&rnd_, 0, 7);
          v_pre.Set(&rnd_, 0, 7);

          y_accum_ref.Set(rnd_.Rand8());
          y_accum_tst.CopyFrom(y_accum_ref);
          y_count_ref.Set(rnd_.Rand8());
          y_count_tst.CopyFrom(y_count_ref);
          u_accum_ref.Set(rnd_.Rand8());
          u_accum_tst.CopyFrom(u_accum_ref);
          u_count_ref.Set(rnd_.Rand8());
          u_count_tst.CopyFrom(u_count_ref);
          v_accum_ref.Set(rnd_.Rand8());
          v_accum_tst.CopyFrom(v_accum_ref);
          v_count_ref.Set(rnd_.Rand8());
          v_count_tst.CopyFrom(v_count_ref);

          ApplyReferenceFilter(y_src, y_pre, u_src, v_src, u_pre, v_pre, width,
                               height, ss_x, ss_y, filter_strength,
                               &filter_weight, use_32x32, &y_accum_ref,
                               &y_count_ref, &u_accum_ref, &u_count_ref,
                               &v_accum_ref, &v_count_ref);
          ASM_REGISTER_STATE_CHECK(filter_func_(
              y_src.TopLeftPixel(), y_src.stride(), y_pre.TopLeftPixel(),
              y_pre.stride(), u_src.TopLeftPixel(), v_src.TopLeftPixel(),
              u_src.stride(), u_pre.TopLeftPixel(), v_pre.TopLeftPixel(),
              u_pre.stride(), width, height, ss_x, ss_y, filter_strength,
              &filter_weight, use_32x32, y_accum_tst.TopLeftPixel(),
              y_count_tst.TopLeftPixel(), u_accum_tst.TopLeftPixel(),
              u_count_tst.TopLeftPixel(), v_accum_tst.TopLeftPixel(),
              v_count_tst.TopLeftPixel()));

          EXPECT_TRUE(y_accum_tst.CheckValues(y_accum_ref));
          EXPECT_TRUE(y_count_tst.CheckValues(y_count_ref));
          EXPECT_TRUE(u_accum_tst.CheckValues(u_accum_ref));
          EXPECT_TRUE(u_count_tst.CheckValues(u_count_ref));
          EXPECT_TRUE(v_accum_tst.CheckValues(v_accum_ref));
          EXPECT_TRUE(v_count_tst.CheckValues(v_count_ref));

          if (HasFailure()) {
            printf("SS_X: %d, SS_Y: %d, Weight: %d, Strength: %d\n", ss_x, ss_y,
                   filter_weight, filter_strength);
            y_accum_tst.PrintDifference(y_accum_ref);
            y_count_tst.PrintDifference(y_count_ref);
            u_accum_tst.PrintDifference(u_accum_ref);
            u_count_tst.PrintDifference(u_count_ref);
            v_accum_tst.PrintDifference(v_accum_ref);
            v_count_tst.PrintDifference(v_count_ref);
            return;
          }
        }
      }
    }
  }
}

TEST_P(YUVTemporalFilterTest, Use_16X16) {
  const int width = 32, height = 32;
  Buffer<uint8_t> y_src = Buffer<uint8_t>(width, height, 8);
  Buffer<uint8_t> y_pre = Buffer<uint8_t>(width, height, 0);
  Buffer<uint16_t> y_count_ref = Buffer<uint16_t>(width, height, 0);
  Buffer<uint32_t> y_accum_ref = Buffer<uint32_t>(width, height, 0);
  Buffer<uint16_t> y_count_tst = Buffer<uint16_t>(width, height, 0);
  Buffer<uint32_t> y_accum_tst = Buffer<uint32_t>(width, height, 0);
  ASSERT_TRUE(y_src.Init());
  ASSERT_TRUE(y_pre.Init());
  ASSERT_TRUE(y_count_ref.Init());
  ASSERT_TRUE(y_accum_ref.Init());
  ASSERT_TRUE(y_count_tst.Init());
  ASSERT_TRUE(y_accum_tst.Init());

  const int use_32x32 = 0;

  for (int ss_x = 0; ss_x <= 1; ss_x++) {
    for (int ss_y = 0; ss_y <= 1; ss_y++) {
      for (int filter_idx = 0; filter_idx < 3 * 3 * 3 * 3; filter_idx++) {
        // Set up the filter
        int filter_weight[4];
        int filter_idx_cp = filter_idx;
        for (int idx = 0; idx < 4; idx++) {
          filter_weight[idx] = filter_idx_cp % 3;
          filter_idx_cp /= 3;
        }

        // Test each parameter
        for (int filter_strength = 0; filter_strength <= 6;
             filter_strength += 2) {
          const int uv_width = width >> ss_x, uv_height = height >> ss_y;
          Buffer<uint8_t> u_src = Buffer<uint8_t>(uv_width, uv_height, 8);
          Buffer<uint8_t> u_pre = Buffer<uint8_t>(uv_width, uv_height, 0);
          Buffer<uint16_t> u_count_ref =
              Buffer<uint16_t>(uv_width, uv_height, 0);
          Buffer<uint32_t> u_accum_ref =
              Buffer<uint32_t>(uv_width, uv_height, 0);
          Buffer<uint16_t> u_count_tst =
              Buffer<uint16_t>(uv_width, uv_height, 0);
          Buffer<uint32_t> u_accum_tst =
              Buffer<uint32_t>(uv_width, uv_height, 0);
          ASSERT_TRUE(u_src.Init());
          ASSERT_TRUE(u_pre.Init());
          ASSERT_TRUE(u_count_ref.Init());
          ASSERT_TRUE(u_accum_ref.Init());
          ASSERT_TRUE(u_count_tst.Init());
          ASSERT_TRUE(u_accum_tst.Init());
          Buffer<uint8_t> v_src = Buffer<uint8_t>(uv_width, uv_height, 8);
          Buffer<uint8_t> v_pre = Buffer<uint8_t>(uv_width, uv_height, 0);
          Buffer<uint16_t> v_count_ref =
              Buffer<uint16_t>(uv_width, uv_height, 0);
          Buffer<uint32_t> v_accum_ref =
              Buffer<uint32_t>(uv_width, uv_height, 0);
          Buffer<uint16_t> v_count_tst =
              Buffer<uint16_t>(uv_width, uv_height, 0);
          Buffer<uint32_t> v_accum_tst =
              Buffer<uint32_t>(uv_width, uv_height, 0);
          ASSERT_TRUE(v_src.Init());
          ASSERT_TRUE(v_pre.Init());
          ASSERT_TRUE(v_count_ref.Init());
          ASSERT_TRUE(v_accum_ref.Init());
          ASSERT_TRUE(v_count_tst.Init());
          ASSERT_TRUE(v_accum_tst.Init());

          // The difference between the buffers must be small to pass the
          // threshold to apply the filter.
          y_src.Set(&rnd_, 0, 7);
          y_pre.Set(&rnd_, 0, 7);
          u_src.Set(&rnd_, 0, 7);
          u_pre.Set(&rnd_, 0, 7);
          v_src.Set(&rnd_, 0, 7);
          v_pre.Set(&rnd_, 0, 7);

          y_accum_ref.Set(rnd_.Rand8());
          y_accum_tst.CopyFrom(y_accum_ref);
          y_count_ref.Set(rnd_.Rand8());
          y_count_tst.CopyFrom(y_count_ref);
          u_accum_ref.Set(rnd_.Rand8());
          u_accum_tst.CopyFrom(u_accum_ref);
          u_count_ref.Set(rnd_.Rand8());
          u_count_tst.CopyFrom(u_count_ref);
          v_accum_ref.Set(rnd_.Rand8());
          v_accum_tst.CopyFrom(v_accum_ref);
          v_count_ref.Set(rnd_.Rand8());
          v_count_tst.CopyFrom(v_count_ref);

          ApplyReferenceFilter(y_src, y_pre, u_src, v_src, u_pre, v_pre, width,
                               height, ss_x, ss_y, filter_strength,
                               filter_weight, use_32x32, &y_accum_ref,
                               &y_count_ref, &u_accum_ref, &u_count_ref,
                               &v_accum_ref, &v_count_ref);
          ASM_REGISTER_STATE_CHECK(filter_func_(
              y_src.TopLeftPixel(), y_src.stride(), y_pre.TopLeftPixel(),
              y_pre.stride(), u_src.TopLeftPixel(), v_src.TopLeftPixel(),
              u_src.stride(), u_pre.TopLeftPixel(), v_pre.TopLeftPixel(),
              u_pre.stride(), width, height, ss_x, ss_y, filter_strength,
              filter_weight, use_32x32, y_accum_tst.TopLeftPixel(),
              y_count_tst.TopLeftPixel(), u_accum_tst.TopLeftPixel(),
              u_count_tst.TopLeftPixel(), v_accum_tst.TopLeftPixel(),
              v_count_tst.TopLeftPixel()));

          EXPECT_TRUE(y_accum_tst.CheckValues(y_accum_ref));
          EXPECT_TRUE(y_count_tst.CheckValues(y_count_ref));
          EXPECT_TRUE(u_accum_tst.CheckValues(u_accum_ref));
          EXPECT_TRUE(u_count_tst.CheckValues(u_count_ref));
          EXPECT_TRUE(v_accum_tst.CheckValues(v_accum_ref));
          EXPECT_TRUE(v_count_tst.CheckValues(v_count_ref));

          if (HasFailure()) {
            printf("SS_X: %d, SS_Y: %d, Weight Idx: %d, Strength: %d\n", ss_x,
                   ss_y, filter_idx, filter_strength);
            y_accum_tst.PrintDifference(y_accum_ref);
            y_count_tst.PrintDifference(y_count_ref);
            u_accum_tst.PrintDifference(u_accum_ref);
            u_count_tst.PrintDifference(u_count_ref);
            v_accum_tst.PrintDifference(v_accum_ref);
            v_count_tst.PrintDifference(v_count_ref);
            return;
          }
        }
      }
    }
  }
}

TEST_P(YUVTemporalFilterTest, DISABLED_Speed) {
  const int width = 32, height = 32;
  Buffer<uint8_t> y_src = Buffer<uint8_t>(width, height, 8);
  Buffer<uint8_t> y_pre = Buffer<uint8_t>(width, height, 0);
  Buffer<uint16_t> y_count = Buffer<uint16_t>(width, height, 0);
  Buffer<uint32_t> y_accum = Buffer<uint32_t>(width, height, 0);
  ASSERT_TRUE(y_src.Init());
  ASSERT_TRUE(y_pre.Init());
  ASSERT_TRUE(y_count.Init());
  ASSERT_TRUE(y_accum.Init());

  for (int use_32x32 = 0; use_32x32 <= 1; use_32x32++) {
    const int num_filter_weights = use_32x32 ? 3 : 3 * 3 * 3 * 3;
    for (int ss_x = 0; ss_x <= 1; ss_x++) {
      for (int ss_y = 0; ss_y <= 1; ss_y++) {
        for (int filter_idx = 0; filter_idx < num_filter_weights;
             filter_idx++) {
          // Set up the filter
          int filter_weight[4];
          int filter_idx_cp = filter_idx;
          for (int idx = 0; idx < 4; idx++) {
            filter_weight[idx] = filter_idx_cp % 3;
            filter_idx_cp /= 3;
          }

          // Test each parameter
          for (int filter_strength = 0; filter_strength <= 6;
               filter_strength += 2) {
            const int uv_width = width >> ss_x, uv_height = height >> ss_y;
            Buffer<uint8_t> u_src = Buffer<uint8_t>(uv_width, uv_height, 8);
            Buffer<uint8_t> u_pre = Buffer<uint8_t>(uv_width, uv_height, 0);
            Buffer<uint16_t> u_count = Buffer<uint16_t>(uv_width, uv_height, 0);
            Buffer<uint32_t> u_accum = Buffer<uint32_t>(uv_width, uv_height, 0);
            ASSERT_TRUE(u_src.Init());
            ASSERT_TRUE(u_pre.Init());
            ASSERT_TRUE(u_count.Init());
            ASSERT_TRUE(u_accum.Init());
            Buffer<uint8_t> v_src = Buffer<uint8_t>(uv_width, uv_height, 8);
            Buffer<uint8_t> v_pre = Buffer<uint8_t>(uv_width, uv_height, 0);
            Buffer<uint16_t> v_count = Buffer<uint16_t>(uv_width, uv_height, 0);
            Buffer<uint32_t> v_accum = Buffer<uint32_t>(uv_width, uv_height, 0);
            ASSERT_TRUE(v_src.Init());
            ASSERT_TRUE(v_pre.Init());
            ASSERT_TRUE(v_count.Init());
            ASSERT_TRUE(v_accum.Init());

            y_src.Set(&rnd_, 0, 7);
            y_pre.Set(&rnd_, 0, 7);
            u_src.Set(&rnd_, 0, 7);
            u_pre.Set(&rnd_, 0, 7);
            v_src.Set(&rnd_, 0, 7);
            v_pre.Set(&rnd_, 0, 7);

            y_accum.Set(0);
            y_count.Set(0);
            u_accum.Set(0);
            u_count.Set(0);
            v_accum.Set(0);
            v_count.Set(0);

            vpx_usec_timer timer;
            vpx_usec_timer_start(&timer);
            for (int num_calls = 0; num_calls < 1000; num_calls++) {
              filter_func_(
                  y_src.TopLeftPixel(), y_src.stride(), y_pre.TopLeftPixel(),
                  y_pre.stride(), u_src.TopLeftPixel(), v_src.TopLeftPixel(),
                  u_src.stride(), u_pre.TopLeftPixel(), v_pre.TopLeftPixel(),
                  u_pre.stride(), width, height, ss_x, ss_y, filter_strength,
                  filter_weight, use_32x32, y_accum.TopLeftPixel(),
                  y_count.TopLeftPixel(), u_accum.TopLeftPixel(),
                  u_count.TopLeftPixel(), v_accum.TopLeftPixel(),
                  v_count.TopLeftPixel());
            }

            vpx_usec_timer_mark(&timer);
            const int elapsed_time =
                static_cast<int>(vpx_usec_timer_elapsed(&timer));

            printf(
                "Use 32X32: %d, SS_X: %d, SS_Y: %d, Weight Idx: %d, Strength: "
                "%d, Time: %5d\n",
                use_32x32, ss_x, ss_y, filter_idx, filter_strength,
                elapsed_time);
          }
        }
      }
    }
  }
}

INSTANTIATE_TEST_CASE_P(C, YUVTemporalFilterTest,
                        ::testing::Values(&vp9_apply_temporal_filter));
}  // namespace
