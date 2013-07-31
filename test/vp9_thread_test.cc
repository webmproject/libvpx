/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/decoder/vp9_thread.h"

#include "third_party/googletest/src/include/gtest/gtest.h"

namespace {

class VP9WorkerThreadTest : public ::testing::Test {
 protected:
  virtual ~VP9WorkerThreadTest() {}
  virtual void SetUp() {
    vp9_worker_init(&worker_);
  }

  virtual void TearDown() {
    vp9_worker_end(&worker_);
  }

  VP9Worker worker_;
};

int ThreadHook(void* data, void* return_value) {
  int* const hook_data = reinterpret_cast<int*>(data);
  *hook_data = 5;
  return *reinterpret_cast<int*>(return_value);
}

TEST_F(VP9WorkerThreadTest, HookSuccess) {
  EXPECT_TRUE(vp9_worker_sync(&worker_));  // should be a no-op.

  for (int i = 0; i < 2; ++i) {
    EXPECT_TRUE(vp9_worker_reset(&worker_));

    int hook_data = 0;
    int return_value = 1;  // return successfully from the hook
    worker_.hook = ThreadHook;
    worker_.data1 = &hook_data;
    worker_.data2 = &return_value;

    vp9_worker_launch(&worker_);
    EXPECT_TRUE(vp9_worker_sync(&worker_));
    EXPECT_FALSE(worker_.had_error);
    EXPECT_EQ(5, hook_data);

    EXPECT_TRUE(vp9_worker_sync(&worker_));  // should be a no-op.
  }
}

TEST_F(VP9WorkerThreadTest, HookFailure) {
  EXPECT_TRUE(vp9_worker_reset(&worker_));

  int hook_data = 0;
  int return_value = 0;  // return failure from the hook
  worker_.hook = ThreadHook;
  worker_.data1 = &hook_data;
  worker_.data2 = &return_value;

  vp9_worker_launch(&worker_);
  EXPECT_FALSE(vp9_worker_sync(&worker_));
  EXPECT_TRUE(worker_.had_error);

  // Ensure _reset() clears the error and _launch() can be called again.
  return_value = 1;
  EXPECT_TRUE(vp9_worker_reset(&worker_));
  EXPECT_FALSE(worker_.had_error);
  vp9_worker_launch(&worker_);
  EXPECT_TRUE(vp9_worker_sync(&worker_));
  EXPECT_FALSE(worker_.had_error);
}

}  // namespace
