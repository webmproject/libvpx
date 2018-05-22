/*
 *  Copyright (c) 2018 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef TEST_BENCH_H_
#define TEST_BENCH_H_

// Number of iterations used to compute median run time.
#define VPX_BENCH_ROBUST_ITER 15

class AbstractBench {
 public:
  void runNTimes(int n);
  void printMedian(const char *title);

 protected:
  // Implement this method and put the code to benchmark in it.
  virtual void run() = 0;

 private:
  int times[VPX_BENCH_ROBUST_ITER];
};

#endif  // TEST_BENCH_H_
